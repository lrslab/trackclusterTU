//! Annotation, gene-context, and bacterial TU semantics.

use std::cmp::Reverse;
use std::collections::{BTreeMap, BinaryHeap, HashSet};
use std::fs::File;
use std::io::Write;
use std::path::Path;

use crate::io::delimited::{DelimitedWriter, Delimiter};
pub(super) use crate::io::gff::GeneRecord;
use crate::model::{Interval, Strand};
use crate::tools::domain_records::{read_bed6_as_genes, Bed6Validation};
use crate::tu::{ReadRecord, Tu, TuEndpointStats};

#[derive(Clone, Copy, Debug)]
pub(super) struct GeneOverlapPolicy {
    pub(super) min_overlap_bp: u32,
    pub(super) min_tu_fraction: f64,
    pub(super) min_gene_fraction: f64,
}

#[derive(Clone, Debug)]
pub(super) struct TuGeneContexts {
    pub(super) same_strand: Vec<Vec<(usize, u32)>>,
    pub(super) antisense: Vec<Vec<(usize, u32)>>,
}

pub(super) fn read_annotation_bed6(path: &Path) -> anyhow::Result<Vec<GeneRecord>> {
    read_bed6_as_genes(
        path,
        Bed6Validation {
            id_label: "annotation",
            require_nonempty_id: false,
            require_unique_ids: false,
            unknown_strand_context: Some("TU gene-context analysis requires '+' or '-'"),
            invalid_interval_label: "gene",
        },
    )
}

fn add_active(
    end: u32,
    idx: usize,
    active: &mut HashSet<usize>,
    ends: &mut BinaryHeap<Reverse<(u32, usize)>>,
) {
    active.insert(idx);
    ends.push(Reverse((end, idx)));
}

fn expire_active(
    current_start: u32,
    active: &mut HashSet<usize>,
    ends: &mut BinaryHeap<Reverse<(u32, usize)>>,
) {
    while let Some(Reverse((end, idx))) = ends.peek().copied() {
        if end <= current_start {
            ends.pop();
            active.remove(&idx);
        } else {
            break;
        }
    }
}

fn overlap_passes_policy(
    tu: &Tu,
    gene: &GeneRecord,
    overlap_bp: u32,
    policy: GeneOverlapPolicy,
) -> bool {
    overlap_bp > 0
        && overlap_bp >= policy.min_overlap_bp
        && (overlap_bp as f64 / tu.interval.len() as f64) >= policy.min_tu_fraction
        && (overlap_bp as f64 / gene.interval.len() as f64) >= policy.min_gene_fraction
}

fn opposite_strand(strand: Strand) -> Strand {
    match strand {
        Strand::Plus => Strand::Minus,
        Strand::Minus => Strand::Plus,
        Strand::Unknown => Strand::Unknown,
    }
}

fn tu_gene_overlaps_for_relation(
    tus: &[Tu],
    genes: &[GeneRecord],
    policy: GeneOverlapPolicy,
    antisense: bool,
) -> Vec<Vec<(usize, u32)>> {
    let mut overlaps: Vec<Vec<(usize, u32)>> = vec![Vec::new(); tus.len()];

    let mut tus_by_key: BTreeMap<(String, Strand), Vec<usize>> = BTreeMap::new();
    for (idx, tu) in tus.iter().enumerate() {
        tus_by_key
            .entry((tu.contig.clone(), tu.strand))
            .or_default()
            .push(idx);
    }

    let mut genes_by_key: BTreeMap<(String, Strand), Vec<usize>> = BTreeMap::new();
    for (idx, gene) in genes.iter().enumerate() {
        let lookup_strand = if antisense {
            opposite_strand(gene.strand)
        } else {
            gene.strand
        };
        genes_by_key
            .entry((gene.contig.clone(), lookup_strand))
            .or_default()
            .push(idx);
    }

    for (key, mut tu_indices) in tus_by_key {
        let Some(gene_indices_ref) = genes_by_key.get(&key) else {
            continue;
        };
        let mut gene_indices = gene_indices_ref.clone();

        tu_indices.sort_by(|&a, &b| {
            tus[a]
                .interval
                .start()
                .cmp(&tus[b].interval.start())
                .then_with(|| tus[a].interval.end().cmp(&tus[b].interval.end()))
                .then_with(|| tus[a].id.cmp(&tus[b].id))
        });

        gene_indices.sort_by(|&a, &b| {
            genes[a]
                .interval
                .start()
                .cmp(&genes[b].interval.start())
                .then_with(|| genes[a].interval.end().cmp(&genes[b].interval.end()))
                .then_with(|| genes[a].id.cmp(&genes[b].id))
        });

        let mut ti: usize = 0;
        let mut gi: usize = 0;

        let mut active_tus: HashSet<usize> = HashSet::new();
        let mut active_genes: HashSet<usize> = HashSet::new();

        let mut tu_ends: BinaryHeap<Reverse<(u32, usize)>> = BinaryHeap::new();
        let mut gene_ends: BinaryHeap<Reverse<(u32, usize)>> = BinaryHeap::new();

        loop {
            let next_tu_start = tu_indices
                .get(ti)
                .map(|&idx| tus[idx].interval.start().get())
                .unwrap_or(u32::MAX);
            let next_gene_start = gene_indices
                .get(gi)
                .map(|&idx| genes[idx].interval.start().get())
                .unwrap_or(u32::MAX);

            if next_tu_start == u32::MAX && next_gene_start == u32::MAX {
                break;
            }

            let current_start = next_tu_start.min(next_gene_start);
            expire_active(current_start, &mut active_tus, &mut tu_ends);
            expire_active(current_start, &mut active_genes, &mut gene_ends);

            if next_tu_start <= next_gene_start {
                let tu_idx = tu_indices[ti];
                ti += 1;

                let tu = &tus[tu_idx];
                if tu.interval.is_empty() {
                    continue;
                }

                add_active(
                    tu.interval.end().get(),
                    tu_idx,
                    &mut active_tus,
                    &mut tu_ends,
                );

                for &gene_idx in &active_genes {
                    let gene = &genes[gene_idx];
                    let overlap = tu.interval.overlap_len(gene.interval);
                    if overlap_passes_policy(tu, gene, overlap, policy) {
                        overlaps[tu_idx].push((gene_idx, overlap));
                    }
                }
            } else {
                let gene_idx = gene_indices[gi];
                gi += 1;

                let gene = &genes[gene_idx];
                if gene.interval.is_empty() {
                    continue;
                }

                add_active(
                    gene.interval.end().get(),
                    gene_idx,
                    &mut active_genes,
                    &mut gene_ends,
                );

                for &tu_idx in &active_tus {
                    let tu = &tus[tu_idx];
                    let overlap = tu.interval.overlap_len(gene.interval);
                    if overlap_passes_policy(tu, gene, overlap, policy) {
                        overlaps[tu_idx].push((gene_idx, overlap));
                    }
                }
            }
        }
    }

    for (tu_idx, overlaps_for_tu) in overlaps.iter_mut().enumerate() {
        overlaps_for_tu.sort_by(|(gene_a, _), (gene_b, _)| {
            let genomic = genes[*gene_a]
                .interval
                .start()
                .cmp(&genes[*gene_b].interval.start())
                .then_with(|| {
                    genes[*gene_a]
                        .interval
                        .end()
                        .cmp(&genes[*gene_b].interval.end())
                });
            let directional = match tus[tu_idx].strand {
                Strand::Plus | Strand::Unknown => genomic,
                Strand::Minus => genomic.reverse(),
            };
            directional.then_with(|| genes[*gene_a].id.cmp(&genes[*gene_b].id))
        });
        overlaps_for_tu.dedup_by(|(gene_a, _), (gene_b, _)| gene_a == gene_b);
    }

    overlaps
}

pub(super) fn tu_gene_contexts(
    tus: &[Tu],
    genes: &[GeneRecord],
    policy: GeneOverlapPolicy,
) -> TuGeneContexts {
    TuGeneContexts {
        same_strand: tu_gene_overlaps_for_relation(tus, genes, policy, false),
        antisense: tu_gene_overlaps_for_relation(tus, genes, policy, true),
    }
}

#[derive(Clone, Copy, Debug, Default)]
pub(super) struct TuClassification {
    alternative_start: bool,
    alternative_termination: bool,
    readthrough: bool,
    antisense: bool,
    intergenic: bool,
    probable_processing_product: bool,
}

impl TuClassification {
    fn labels(self) -> Vec<&'static str> {
        let mut labels = Vec::new();
        if self.alternative_start {
            labels.push("alternative_start");
        }
        if self.alternative_termination {
            labels.push("alternative_termination");
        }
        if self.readthrough {
            labels.push("readthrough");
        }
        if self.antisense {
            labels.push("antisense");
        }
        if self.intergenic {
            labels.push("intergenic");
        }
        if self.probable_processing_product {
            labels.push("probable_processing_product");
        }
        if labels.is_empty() {
            labels.push("canonical");
        }
        labels
    }
}

fn tu_five_prime_coord(tu: &Tu) -> u32 {
    match tu.strand {
        Strand::Plus | Strand::Unknown => tu.interval.start().get(),
        Strand::Minus => tu.interval.end().get(),
    }
}

fn tu_three_prime_coord(tu: &Tu) -> u32 {
    match tu.strand {
        Strand::Plus | Strand::Unknown => tu.interval.end().get(),
        Strand::Minus => tu.interval.start().get(),
    }
}

pub(super) fn classify_tus(tus: &[Tu], contexts: &TuGeneContexts) -> Vec<TuClassification> {
    let mut classifications = vec![TuClassification::default(); tus.len()];
    for (tu_idx, classification) in classifications.iter_mut().enumerate() {
        classification.readthrough = contexts.same_strand[tu_idx].len() >= 2;
        classification.antisense = !contexts.antisense[tu_idx].is_empty();
        classification.intergenic =
            contexts.same_strand[tu_idx].is_empty() && contexts.antisense[tu_idx].is_empty();
    }

    let mut comparable_tus: BTreeMap<(String, Strand, Vec<usize>), Vec<usize>> = BTreeMap::new();
    for (tu_idx, tu) in tus.iter().enumerate() {
        let signature: Vec<usize> = contexts.same_strand[tu_idx]
            .iter()
            .map(|(gene_idx, _)| *gene_idx)
            .collect();
        if !signature.is_empty() {
            comparable_tus
                .entry((tu.contig.clone(), tu.strand, signature))
                .or_default()
                .push(tu_idx);
        }
    }

    for indices in comparable_tus.into_values() {
        for left_position in 0..indices.len() {
            for right_position in (left_position + 1)..indices.len() {
                let left_idx = indices[left_position];
                let right_idx = indices[right_position];
                let left = &tus[left_idx];
                let right = &tus[right_idx];

                let left_five_prime = tu_five_prime_coord(left);
                let right_five_prime = tu_five_prime_coord(right);
                let left_three_prime = tu_three_prime_coord(left);
                let right_three_prime = tu_three_prime_coord(right);
                if left_three_prime == right_three_prime && left_five_prime != right_five_prime {
                    classifications[left_idx].alternative_start = true;
                    classifications[right_idx].alternative_start = true;
                }
                if left_five_prime == right_five_prime && left_three_prime != right_three_prime {
                    classifications[left_idx].alternative_termination = true;
                    classifications[right_idx].alternative_termination = true;
                }

                let left_inside_right = right.interval.start() <= left.interval.start()
                    && right.interval.end() >= left.interval.end()
                    && right.interval != left.interval;
                let right_inside_left = left.interval.start() <= right.interval.start()
                    && left.interval.end() >= right.interval.end()
                    && right.interval != left.interval;
                if left_five_prime != right_five_prime && left_three_prime != right_three_prime {
                    if left_inside_right {
                        classifications[left_idx].probable_processing_product = true;
                    }
                    if right_inside_left {
                        classifications[right_idx].probable_processing_product = true;
                    }
                }
            }
        }
    }

    classifications
}

fn escaped_gene_signature(hits: &[(usize, u32)], genes: &[GeneRecord]) -> String {
    if hits.is_empty() {
        return ".".to_owned();
    }
    hits.iter()
        .map(|(gene_idx, _)| percent_encode(genes[*gene_idx].id.as_bytes()))
        .collect::<Vec<_>>()
        .join(",")
}

fn percent_encode(bytes: &[u8]) -> String {
    const HEX: &[u8; 16] = b"0123456789ABCDEF";
    let mut escaped = String::with_capacity(bytes.len());
    for &byte in bytes {
        if byte.is_ascii_alphanumeric()
            || matches!(
                byte,
                b'.' | b'_' | b':' | b'^' | b'*' | b'$' | b'@' | b'!' | b'+' | b'?' | b'-' | b'|'
            )
        {
            escaped.push(byte as char);
        } else {
            escaped.push('%');
            escaped.push(HEX[(byte >> 4) as usize] as char);
            escaped.push(HEX[(byte & 0x0f) as usize] as char);
        }
    }
    escaped
}

pub(super) fn write_tu_semantics(
    path: &Path,
    tus: &[Tu],
    genes: &[GeneRecord],
    contexts: &TuGeneContexts,
    classifications: &[TuClassification],
    policy: GeneOverlapPolicy,
) -> anyhow::Result<()> {
    let metadata = vec![
        "#trackclustertu_tu_semantics_schema=v1".to_owned(),
        "#gene_context_order=transcription_direction".to_owned(),
        format!("#gene_min_overlap_bp={}", policy.min_overlap_bp),
        format!("#gene_min_tu_fraction={}", policy.min_tu_fraction),
        format!("#gene_min_gene_fraction={}", policy.min_gene_fraction),
        "#alternative_start=same_nonempty_gene_context_and_same_3prime_but_different_5prime".to_owned(),
        "#alternative_termination=same_nonempty_gene_context_and_same_5prime_but_different_3prime".to_owned(),
        "#readthrough=two_or_more_qualifying_same_strand_genes".to_owned(),
        "#antisense=one_or_more_qualifying_opposite_strand_genes".to_owned(),
        "#intergenic=no_qualifying_gene_on_either_strand".to_owned(),
        "#probable_processing_product=strictly_contained_same_nonempty_gene_context_with_both_boundaries_internal".to_owned(),
        "#classification_order=alternative_start,alternative_termination,readthrough,antisense,intergenic,probable_processing_product".to_owned(),
    ];
    let mut writer = DelimitedWriter::create(path, Delimiter::Tab, &metadata)?;
    writer.write_record([
        "tu_id",
        "contig",
        "strand",
        "start",
        "end",
        "gene_context_signature",
        "same_strand_gene_count",
        "antisense_gene_context",
        "antisense_gene_count",
        "classifications",
    ])?;
    for (tu_idx, tu) in tus.iter().enumerate() {
        writer.write_record([
            tu.id.clone(),
            tu.contig.clone(),
            tu.strand.as_char().to_string(),
            tu.interval.start().get().to_string(),
            tu.interval.end().get().to_string(),
            escaped_gene_signature(&contexts.same_strand[tu_idx], genes),
            contexts.same_strand[tu_idx].len().to_string(),
            escaped_gene_signature(&contexts.antisense[tu_idx], genes),
            contexts.antisense[tu_idx].len().to_string(),
            classifications[tu_idx].labels().join(","),
        ])?;
    }
    writer.flush()?;
    Ok(())
}

pub(super) fn write_tu_gene_relationships(
    path: &Path,
    tus: &[Tu],
    genes: &[GeneRecord],
    contexts: &TuGeneContexts,
    policy: GeneOverlapPolicy,
) -> anyhow::Result<()> {
    let metadata = vec![
        "#trackclustertu_tu_gene_schema=v2".to_owned(),
        "#gene_context_order=transcription_direction".to_owned(),
        format!("#gene_min_overlap_bp={}", policy.min_overlap_bp),
        format!("#gene_min_tu_fraction={}", policy.min_tu_fraction),
        format!("#gene_min_gene_fraction={}", policy.min_gene_fraction),
        "#contig\tstrand\ttu_id\ttu_start\ttu_end\trelation\tcontext_order\tgene_id\tgene_start\tgene_end\toverlap_bp\ttu_fraction\tgene_fraction".to_owned(),
    ];
    let mut writer = DelimitedWriter::create(path, Delimiter::Tab, &metadata)?;
    for (tu_idx, tu) in tus.iter().enumerate() {
        for (relation, hits) in [
            ("same_strand", &contexts.same_strand[tu_idx]),
            ("antisense", &contexts.antisense[tu_idx]),
        ] {
            for (order, &(gene_idx, overlap_bp)) in hits.iter().enumerate() {
                let gene = &genes[gene_idx];
                writer.write_record([
                    tu.contig.clone(),
                    tu.strand.as_char().to_string(),
                    tu.id.clone(),
                    tu.interval.start().get().to_string(),
                    tu.interval.end().get().to_string(),
                    relation.to_owned(),
                    (order + 1).to_string(),
                    gene.id.clone(),
                    gene.interval.start().get().to_string(),
                    gene.interval.end().get().to_string(),
                    overlap_bp.to_string(),
                    format!("{:.6}", overlap_bp as f64 / tu.interval.len() as f64),
                    format!("{:.6}", overlap_bp as f64 / gene.interval.len() as f64),
                ])?;
            }
        }
    }
    writer.flush()?;
    Ok(())
}

pub(super) fn write_tu_gff3(
    path: &Path,
    tus: &[Tu],
    genes: &[GeneRecord],
    contexts: &TuGeneContexts,
    classifications: &[TuClassification],
    endpoint_stats: &[TuEndpointStats],
) -> anyhow::Result<()> {
    let mut writer = std::io::BufWriter::new(File::create(path)?);
    writeln!(writer, "##gff-version 3")?;
    writeln!(writer, "#trackclustertu_tu_gff3_schema=v1")?;
    for (tu_idx, tu) in tus.iter().enumerate() {
        let escaped_tu_id = percent_encode(tu.id.as_bytes());
        let context = escaped_gene_signature(&contexts.same_strand[tu_idx], genes);
        let classifications = classifications[tu_idx].labels().join(",");
        writeln!(
            writer,
            "{}\ttrackclustertu\ttranscript\t{}\t{}\t.\t{}\t.\tID={};Name={};gene_context={};classification={};support={}",
            percent_encode(tu.contig.as_bytes()),
            tu.interval.start().get() as u64 + 1,
            tu.interval.end().get(),
            tu.strand.as_char(),
            escaped_tu_id,
            escaped_tu_id,
            context,
            classifications,
            endpoint_stats[tu_idx].support,
        )?;

        let mut relation_number = 0usize;
        for (relation, hits) in [
            ("same_strand", &contexts.same_strand[tu_idx]),
            ("antisense", &contexts.antisense[tu_idx]),
        ] {
            for &(gene_idx, overlap_bp) in hits {
                relation_number += 1;
                let gene = &genes[gene_idx];
                let overlap = tu
                    .interval
                    .intersection(gene.interval)
                    .expect("reported gene relationship must overlap");
                let relationship_id = percent_encode(
                    format!("{}.gene_overlap.{}", tu.id, relation_number).as_bytes(),
                );
                writeln!(
                    writer,
                    "{}\ttrackclustertu\tgene_overlap\t{}\t{}\t.\t{}\t.\tID={};Parent={};gene_id={};relation={};overlap_bp={}",
                    percent_encode(tu.contig.as_bytes()),
                    overlap.start().get() as u64 + 1,
                    overlap.end().get(),
                    tu.strand.as_char(),
                    relationship_id,
                    escaped_tu_id,
                    percent_encode(gene.id.as_bytes()),
                    relation,
                    overlap_bp,
                )?;
            }
        }
    }
    writer.flush()?;
    Ok(())
}

pub(super) fn format_name2(
    read_indices: &[usize],
    reads: &[ReadRecord],
    read_count: usize,
) -> String {
    if read_indices.is_empty() {
        return format!(",|{read_count}");
    }

    let mut name2 = String::new();
    for (i, &read_idx) in read_indices.iter().enumerate() {
        if i > 0 {
            name2.push(',');
        }
        name2.push_str(&reads[read_idx].id);
    }
    name2.push_str(",|");
    name2.push_str(&read_count.to_string());
    name2
}

pub(super) fn merge_overlapping_intervals(mut intervals: Vec<Interval>) -> Vec<Interval> {
    intervals.sort_by(|left, right| {
        left.start()
            .cmp(&right.start())
            .then_with(|| left.end().cmp(&right.end()))
    });

    let mut merged: Vec<Interval> = Vec::with_capacity(intervals.len());
    for interval in intervals {
        let Some(previous) = merged.last_mut() else {
            merged.push(interval);
            continue;
        };
        if interval.start() < previous.end() {
            if interval.end() > previous.end() {
                *previous = Interval::new(previous.start(), interval.end())
                    .expect("merged interval preserves coordinate order");
            }
        } else {
            merged.push(interval);
        }
    }
    merged
}
