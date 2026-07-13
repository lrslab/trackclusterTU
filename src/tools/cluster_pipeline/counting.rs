//! Membership parsing and TU/sample/group/gene counting.

use std::collections::{BTreeSet, HashSet};
use std::path::Path;

use super::annotation::GeneRecord;
use crate::io::delimited::{DelimitedWriter, Delimiter};
use crate::tools::membership::read_memberships;
use crate::tu::multi::{MultiSampleCounter, MultiSampleCounts, SampleManifestRecord};

#[derive(Clone, Debug)]
pub(super) struct CountMetrics {
    pub(super) unique: Vec<u64>,
    pub(super) full_length_evidence: Vec<u64>,
    pub(super) total: Vec<u64>,
    pub(super) fractional: Vec<f64>,
}

impl CountMetrics {
    pub(super) fn zeros(entity_count: usize) -> Self {
        Self {
            unique: vec![0; entity_count],
            full_length_evidence: vec![0; entity_count],
            total: vec![0; entity_count],
            fractional: vec![0.0; entity_count],
        }
    }
}

pub(super) fn write_count_metrics_csv(
    path: &Path,
    id_header: &str,
    ids: &[String],
    counts: &CountMetrics,
) -> anyhow::Result<()> {
    let mut writer = DelimitedWriter::create(path, Delimiter::Comma, &[])?;
    // `count` remains an alias of `total_count` for readers of the historical two-column CSV.
    writer.write_record([
        id_header,
        "count",
        "unique_count",
        "full_length_evidence_count",
        "total_count",
        "fractional_count",
    ])?;
    for (entity_idx, id) in ids.iter().enumerate() {
        writer.write_record([
            id.clone(),
            counts.total[entity_idx].to_string(),
            counts.unique[entity_idx].to_string(),
            counts.full_length_evidence[entity_idx].to_string(),
            counts.total[entity_idx].to_string(),
            counts.fractional[entity_idx].to_string(),
        ])?;
    }
    writer.flush()?;
    Ok(())
}

pub(super) const GENE_COUNT_SEMANTICS: &str = "nonexclusive_same_strand_tu_overlap;each_TU_assignment_is_counted_once_for_every_qualifying_gene;gene_totals_may_exceed_read_and_TU_totals;antisense_relationships_are_not_counted";

pub(super) fn write_gene_count_metrics_csv(
    path: &Path,
    ids: &[String],
    counts: &CountMetrics,
) -> anyhow::Result<()> {
    let mut writer = DelimitedWriter::create(
        path,
        Delimiter::Comma,
        &[format!("#count_semantics={GENE_COUNT_SEMANTICS}")],
    )?;
    writer.write_record([
        "gene_id",
        "count",
        "unique_count",
        "full_length_evidence_count",
        "total_count",
        "fractional_count",
    ])?;
    for (entity_idx, id) in ids.iter().enumerate() {
        writer.write_record([
            id.clone(),
            counts.total[entity_idx].to_string(),
            counts.unique[entity_idx].to_string(),
            counts.full_length_evidence[entity_idx].to_string(),
            counts.total[entity_idx].to_string(),
            counts.fractional[entity_idx].to_string(),
        ])?;
    }
    writer.flush()?;
    Ok(())
}

#[allow(clippy::too_many_arguments)]
pub(super) fn write_entity_metric_matrix(
    path: &Path,
    count_semantics: Option<&str>,
    id_header: &str,
    ids: &[String],
    column_names: &[String],
    unique_counts: &[Vec<u64>],
    full_length_evidence_counts: &[Vec<u64>],
    total_counts: &[Vec<u64>],
    fractional_counts: &[Vec<f64>],
) -> anyhow::Result<()> {
    let metadata = count_semantics
        .map(|value| vec![format!("#count_semantics={value}")])
        .unwrap_or_default();
    let mut writer = DelimitedWriter::create(path, Delimiter::Tab, &metadata)?;
    let mut header = vec![id_header.to_owned()];
    for column_name in column_names {
        for metric in [
            "unique_count",
            "full_length_evidence_count",
            "total_count",
            "fractional_count",
        ] {
            header.push(format!("{column_name}.{metric}"));
        }
    }
    writer.write_record(header)?;

    for (row_idx, id) in ids.iter().enumerate() {
        let mut row = vec![id.clone()];
        for column_idx in 0..column_names.len() {
            row.push(unique_counts[row_idx][column_idx].to_string());
            row.push(full_length_evidence_counts[row_idx][column_idx].to_string());
            row.push(total_counts[row_idx][column_idx].to_string());
            row.push(fractional_counts[row_idx][column_idx].to_string());
        }
        writer.write_record(row)?;
    }
    writer.flush()?;
    Ok(())
}

#[allow(clippy::too_many_arguments)]
pub(super) fn write_entity_metric_long_counts(
    path: &Path,
    id_header: &str,
    column_header: &str,
    ids: &[String],
    column_names: &[String],
    unique_counts: &[Vec<u64>],
    full_length_evidence_counts: &[Vec<u64>],
    total_counts: &[Vec<u64>],
    fractional_counts: &[Vec<f64>],
) -> anyhow::Result<()> {
    let mut writer = DelimitedWriter::create(path, Delimiter::Tab, &[])?;
    writer.write_record([
        id_header,
        column_header,
        "unique_count",
        "full_length_evidence_count",
        "total_count",
        "fractional_count",
    ])?;
    for (row_idx, id) in ids.iter().enumerate() {
        for (column_idx, column_name) in column_names.iter().enumerate() {
            let unique = unique_counts[row_idx][column_idx];
            let full_length = full_length_evidence_counts[row_idx][column_idx];
            let total = total_counts[row_idx][column_idx];
            let fractional = fractional_counts[row_idx][column_idx];
            if unique > 0 || full_length > 0 || total > 0 || fractional > 0.0 {
                writer.write_record([
                    id.clone(),
                    column_name.clone(),
                    unique.to_string(),
                    full_length.to_string(),
                    total.to_string(),
                    fractional.to_string(),
                ])?;
            }
        }
    }
    writer.flush()?;
    Ok(())
}

pub(super) fn total_counts_from_multi(counts: &MultiSampleCounts) -> CountMetrics {
    CountMetrics {
        unique: counts
            .unique_counts()
            .iter()
            .map(|row| row.iter().sum())
            .collect(),
        full_length_evidence: counts
            .full_length_evidence_counts()
            .iter()
            .map(|row| row.iter().sum())
            .collect(),
        total: counts.counts().iter().map(|row| row.iter().sum()).collect(),
        fractional: counts
            .fractional_counts()
            .iter()
            .map(|row| row.iter().sum())
            .collect(),
    }
}

pub(super) fn filter_multi_counts(
    counts: MultiSampleCounts,
    min_tu_count: Option<u64>,
) -> MultiSampleCounts {
    let Some(min_tu_count) = min_tu_count else {
        return counts;
    };
    counts.filter_min_total_count(min_tu_count)
}

pub(super) fn build_gene_total_counts(
    genes: &[GeneRecord],
    overlaps: &[Vec<(usize, u32)>],
    tu_counts: &CountMetrics,
) -> CountMetrics {
    let mut gene_counts = CountMetrics::zeros(genes.len());
    for (tu_idx, gene_hits) in overlaps.iter().enumerate() {
        let mut seen = HashSet::new();
        for &(gene_idx, _) in gene_hits {
            if seen.insert(gene_idx) {
                gene_counts.unique[gene_idx] += tu_counts.unique[tu_idx];
                gene_counts.full_length_evidence[gene_idx] +=
                    tu_counts.full_length_evidence[tu_idx];
                gene_counts.total[gene_idx] += tu_counts.total[tu_idx];
                gene_counts.fractional[gene_idx] += tu_counts.fractional[tu_idx];
            }
        }
    }
    gene_counts
}

pub(super) struct MatrixCountMetrics {
    pub(super) unique: Vec<Vec<u64>>,
    pub(super) full_length_evidence: Vec<Vec<u64>>,
    pub(super) total: Vec<Vec<u64>>,
    pub(super) fractional: Vec<Vec<f64>>,
}

impl MatrixCountMetrics {
    fn zeros(rows: usize, columns: usize) -> Self {
        Self {
            unique: vec![vec![0; columns]; rows],
            full_length_evidence: vec![vec![0; columns]; rows],
            total: vec![vec![0; columns]; rows],
            fractional: vec![vec![0.0; columns]; rows],
        }
    }
}

pub(super) fn build_gene_multi_counts(
    genes: &[GeneRecord],
    overlaps: &[Vec<(usize, u32)>],
    tu_counts: &MultiSampleCounts,
) -> (MatrixCountMetrics, MatrixCountMetrics) {
    let mut gene_sample_counts = MatrixCountMetrics::zeros(genes.len(), tu_counts.samples().len());
    let mut gene_group_counts = MatrixCountMetrics::zeros(genes.len(), tu_counts.groups().len());

    for (tu_idx, gene_hits) in overlaps.iter().enumerate() {
        let mut seen = HashSet::new();
        for &(gene_idx, _) in gene_hits {
            if !seen.insert(gene_idx) {
                continue;
            }
            for sample_idx in 0..tu_counts.samples().len() {
                gene_sample_counts.unique[gene_idx][sample_idx] +=
                    tu_counts.unique_counts()[tu_idx][sample_idx];
                gene_sample_counts.full_length_evidence[gene_idx][sample_idx] +=
                    tu_counts.full_length_evidence_counts()[tu_idx][sample_idx];
                gene_sample_counts.total[gene_idx][sample_idx] +=
                    tu_counts.counts()[tu_idx][sample_idx];
                gene_sample_counts.fractional[gene_idx][sample_idx] +=
                    tu_counts.fractional_counts()[tu_idx][sample_idx];
            }
            for group_idx in 0..tu_counts.groups().len() {
                gene_group_counts.unique[gene_idx][group_idx] +=
                    tu_counts.group_unique_counts()[tu_idx][group_idx];
                gene_group_counts.full_length_evidence[gene_idx][group_idx] +=
                    tu_counts.group_full_length_evidence_counts()[tu_idx][group_idx];
                gene_group_counts.total[gene_idx][group_idx] +=
                    tu_counts.group_counts()[tu_idx][group_idx];
                gene_group_counts.fractional[gene_idx][group_idx] +=
                    tu_counts.group_fractional_counts()[tu_idx][group_idx];
            }
        }
    }

    (gene_sample_counts, gene_group_counts)
}

pub(super) fn build_multi_counts_from_membership(
    samples: &[SampleManifestRecord],
    membership_path: &Path,
    min_tu_count: Option<u64>,
) -> anyhow::Result<MultiSampleCounts> {
    let memberships = read_memberships(membership_path)?;
    let tu_ids = memberships
        .iter()
        .flat_map(|membership| membership.candidate_tu_ids())
        .map(str::to_owned)
        .collect::<BTreeSet<_>>()
        .into_iter()
        .collect();
    let mut counter = MultiSampleCounter::new(samples, tu_ids);
    for membership in &memberships {
        for contribution in membership.contributions() {
            counter.add_assignment_metrics(
                &membership.read_id,
                contribution.tu_id,
                contribution.unique,
                contribution.full_length_evidence,
                contribution.total,
                contribution.fractional_weight,
            )?;
        }
    }
    Ok(filter_multi_counts(counter.finish(), min_tu_count))
}
