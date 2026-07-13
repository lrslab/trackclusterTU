//! Transaction-safe output records and versioned schema writers.

use std::fs::File;
use std::io::Write;
use std::path::Path;

use crate::io::delimited::{DelimitedWriter, Delimiter};
use crate::model::{Interval, Strand};
use crate::tools::membership::{
    membership_metadata, write_membership, MembershipRecord, MembershipStatus, V2MembershipFields,
    WeightedTu,
};
use crate::tu::{
    AssignmentStatus, ReadAssignment, ReadRecord, Tu, TuClusteringResult, TuEndpointStats,
};

use super::config::TuIdStyle;

#[derive(Clone, Debug)]
pub(super) struct PreparedOutputs {
    pub(super) tus: Vec<Tu>,
    pub(super) endpoint_stats: Vec<TuEndpointStats>,
    pub(super) id_mappings: Vec<TuIdMapping>,
}

#[derive(Clone, Debug)]
pub(super) struct TuIdMapping {
    emitted_tu_id: String,
    stable_tu_id: String,
    sequential_tu_id: String,
    contig: String,
    strand: Strand,
    interval: Interval,
}

fn format_tu_id(index: usize, width: usize) -> String {
    format!("TU{:0width$}", index + 1, width = width)
}

fn stable_tu_id(tu: &Tu) -> String {
    const HEX: &[u8; 16] = b"0123456789abcdef";
    let mut encoded_contig = String::with_capacity(tu.contig.len() * 2);
    for byte in tu.contig.as_bytes() {
        encoded_contig.push(HEX[(byte >> 4) as usize] as char);
        encoded_contig.push(HEX[(byte & 0x0f) as usize] as char);
    }
    let strand = match tu.strand {
        Strand::Plus => 'p',
        Strand::Minus => 'm',
        Strand::Unknown => 'u',
    };
    format!(
        "TUg_{}_{}_{}_{}",
        encoded_contig,
        tu.interval.start().get(),
        tu.interval.end().get(),
        strand
    )
}

pub(super) fn write_pooled_reads_bed6(path: &Path, reads: &[ReadRecord]) -> anyhow::Result<()> {
    let records: Vec<crate::io::bed::Bed6Record> = reads
        .iter()
        .map(|read| crate::io::bed::Bed6Record {
            chrom: read.contig.clone(),
            start: read.interval.start(),
            end: read.interval.end(),
            name: read.id.clone(),
            score: 0,
            strand: read.strand,
            extra_fields: Vec::new(),
        })
        .collect();
    crate::io::bed::write_bed6(path, records.iter())?;
    Ok(())
}

pub(super) fn write_tu_endpoint_stats(
    path: &Path,
    prepared: &PreparedOutputs,
    reads: &[ReadRecord],
) -> anyhow::Result<()> {
    let mut writer = DelimitedWriter::create(
        path,
        Delimiter::Tab,
        &["#trackclustertu_tu_endpoint_stats_schema=v1".to_owned()],
    )?;
    writer.write_record([
        "tu_id",
        "support",
        "five_prime_consensus",
        "three_prime_consensus",
        "five_prime_min",
        "five_prime_max",
        "five_prime_spread_bp",
        "three_prime_min",
        "three_prime_max",
        "three_prime_spread_bp",
        "example_read_id",
    ])?;
    for (tu, stats) in prepared.tus.iter().zip(prepared.endpoint_stats.iter()) {
        writer.write_record([
            tu.id.clone(),
            stats.support.to_string(),
            stats.five_prime_consensus.to_string(),
            stats.three_prime_consensus.to_string(),
            stats.five_prime_min.to_string(),
            stats.five_prime_max.to_string(),
            stats.five_prime_spread_bp().to_string(),
            stats.three_prime_min.to_string(),
            stats.three_prime_max.to_string(),
            stats.three_prime_spread_bp().to_string(),
            reads[tu.rep_read_index].id.clone(),
        ])?;
    }
    writer.flush()?;
    Ok(())
}

pub(super) fn write_tu_id_map(path: &Path, prepared: &PreparedOutputs) -> anyhow::Result<()> {
    let mut writer = DelimitedWriter::create(
        path,
        Delimiter::Tab,
        &["#trackclustertu_tu_id_map_schema=v1".to_owned()],
    )?;
    writer.write_record([
        "emitted_tu_id",
        "stable_tu_id",
        "sequential_tu_id",
        "contig",
        "strand",
        "start",
        "end",
    ])?;
    for mapping in &prepared.id_mappings {
        writer.write_record([
            mapping.emitted_tu_id.clone(),
            mapping.stable_tu_id.clone(),
            mapping.sequential_tu_id.clone(),
            mapping.contig.clone(),
            mapping.strand.as_char().to_string(),
            mapping.interval.start().get().to_string(),
            mapping.interval.end().get().to_string(),
        ])?;
    }
    writer.flush()?;
    Ok(())
}

pub(super) fn hard_assignment_tu_index(assignment: &ReadAssignment) -> Option<usize> {
    match assignment.status {
        AssignmentStatus::Unique | AssignmentStatus::Partial => {
            assignment.best.map(|candidate| candidate.tu_index)
        }
        AssignmentStatus::Ambiguous | AssignmentStatus::Unassigned => None,
    }
}

pub(super) fn create_membership_writer(
    path: &Path,
    ambiguity_margin: f64,
    fractional_assignment: bool,
) -> anyhow::Result<DelimitedWriter<std::io::BufWriter<File>>> {
    Ok(DelimitedWriter::create(
        path,
        Delimiter::Tab,
        &membership_metadata(ambiguity_margin, fractional_assignment),
    )?)
}

pub(super) fn write_membership_row<W: Write>(
    writer: &mut DelimitedWriter<W>,
    read: &ReadRecord,
    assignment: &ReadAssignment,
    tus: &[Tu],
    full_length_evidence: bool,
) -> anyhow::Result<()> {
    let best = assignment.best;
    let second = assignment.second_best;
    let hard_tu = hard_assignment_tu_index(assignment);
    let primary_tu_id = hard_tu
        .map(|tu_index| tus[tu_index].id.as_str())
        .unwrap_or(".");
    let best_tu_id = best
        .map(|candidate| tus[candidate.tu_index].id.as_str())
        .unwrap_or(".");
    let second_tu_id = second
        .map(|candidate| tus[candidate.tu_index].id.as_str())
        .unwrap_or(".");
    let best_score1 = best
        .map(|candidate| format!("{:.6}", candidate.score1))
        .unwrap_or_else(|| ".".to_owned());
    let best_score2 = best
        .map(|candidate| format!("{:.6}", candidate.score2))
        .unwrap_or_else(|| ".".to_owned());
    let second_score1 = second
        .map(|candidate| format!("{:.6}", candidate.score1))
        .unwrap_or_else(|| ".".to_owned());
    let second_score2 = second
        .map(|candidate| format!("{:.6}", candidate.score2))
        .unwrap_or_else(|| ".".to_owned());
    let best_five_prime_delta = best
        .map(|candidate| candidate.five_prime_delta_bp.to_string())
        .unwrap_or_else(|| ".".to_owned());
    let best_three_prime_delta = best
        .map(|candidate| candidate.three_prime_delta_bp.to_string())
        .unwrap_or_else(|| ".".to_owned());
    let second_five_prime_delta = second
        .map(|candidate| candidate.five_prime_delta_bp.to_string())
        .unwrap_or_else(|| ".".to_owned());
    let second_three_prime_delta = second
        .map(|candidate| candidate.three_prime_delta_bp.to_string())
        .unwrap_or_else(|| ".".to_owned());
    let best_assignment_score = best
        .map(|candidate| format!("{:.6}", candidate.assignment_score))
        .unwrap_or_else(|| ".".to_owned());
    let second_assignment_score = second
        .map(|candidate| format!("{:.6}", candidate.assignment_score))
        .unwrap_or_else(|| ".".to_owned());
    let score_margin = assignment
        .score_margin
        .map(|margin| format!("{margin:.6}"))
        .unwrap_or_else(|| ".".to_owned());
    let status = match assignment.status {
        AssignmentStatus::Unique => MembershipStatus::Unique,
        AssignmentStatus::Ambiguous => MembershipStatus::Ambiguous,
        AssignmentStatus::Partial => MembershipStatus::Partial,
        AssignmentStatus::Unassigned => MembershipStatus::Unassigned,
    };
    let primary_weight = match status {
        MembershipStatus::Unique | MembershipStatus::Partial => 1.0,
        MembershipStatus::Ambiguous => best
            .and_then(|best| {
                assignment
                    .fractional_assignments
                    .iter()
                    .find(|fraction| fraction.tu_index == best.tu_index)
                    .map(|fraction| fraction.weight)
            })
            .unwrap_or(0.0),
        MembershipStatus::Unassigned => 0.0,
        MembershipStatus::Legacy => unreachable!("new clustering emits v2 membership rows"),
    };
    let record = MembershipRecord {
        read_id: read.id.clone(),
        hard_tu_id: (primary_tu_id != ".").then(|| primary_tu_id.to_owned()),
        score1: best_score1,
        score2: best_score2,
        status,
        fractional_assignments: assignment
            .fractional_assignments
            .iter()
            .map(|fraction| WeightedTu {
                tu_id: tus[fraction.tu_index].id.clone(),
                weight: fraction.weight,
            })
            .collect(),
        v2: Some(V2MembershipFields {
            complete: true,
            best_candidate_tu_id: (best_tu_id != ".").then(|| best_tu_id.to_owned()),
            second_candidate_tu_id: (second_tu_id != ".").then(|| second_tu_id.to_owned()),
            second_score1,
            second_score2,
            best_five_prime_delta_bp: best_five_prime_delta,
            best_three_prime_delta_bp: best_three_prime_delta,
            second_five_prime_delta_bp: second_five_prime_delta,
            second_three_prime_delta_bp: second_three_prime_delta,
            best_assignment_score,
            second_assignment_score,
            score_margin,
            primary_weight: primary_weight.to_string(),
            full_length_evidence,
        }),
        line_number: 0,
    };
    write_membership(writer, &record)?;
    Ok(())
}

pub(super) fn prepare_outputs(
    result: &TuClusteringResult,
    min_tu_count: Option<u64>,
    id_style: TuIdStyle,
) -> PreparedOutputs {
    let mut tu_counts: Vec<u64> = vec![0; result.tus().len()];
    for &tu_idx in result.read_to_tu() {
        tu_counts[tu_idx] += 1;
    }

    let mut kept_old_tu_indices: Vec<usize> = Vec::new();
    for (tu_idx, &count) in tu_counts.iter().enumerate() {
        if min_tu_count.is_none_or(|min| count >= min) {
            kept_old_tu_indices.push(tu_idx);
        }
    }

    let new_width = 6usize.max(kept_old_tu_indices.len().to_string().len());
    let mut tus: Vec<Tu> = Vec::with_capacity(kept_old_tu_indices.len());
    let mut endpoint_stats: Vec<TuEndpointStats> = Vec::with_capacity(kept_old_tu_indices.len());
    let mut id_mappings: Vec<TuIdMapping> = Vec::with_capacity(kept_old_tu_indices.len());

    for (new_idx, &old_idx) in kept_old_tu_indices.iter().enumerate() {
        let tu = &result.tus()[old_idx];
        let sequential_tu_id = format_tu_id(new_idx, new_width);
        let mut emitted_tu = Tu {
            id: String::new(),
            contig: tu.contig.clone(),
            strand: tu.strand,
            interval: tu.interval,
            rep_read_index: tu.rep_read_index,
        };
        let stable_tu_id = stable_tu_id(&emitted_tu);
        emitted_tu.id = match id_style {
            TuIdStyle::Stable => stable_tu_id.clone(),
            TuIdStyle::Sequential => sequential_tu_id.clone(),
        };
        id_mappings.push(TuIdMapping {
            emitted_tu_id: emitted_tu.id.clone(),
            stable_tu_id,
            sequential_tu_id,
            contig: emitted_tu.contig.clone(),
            strand: emitted_tu.strand,
            interval: emitted_tu.interval,
        });
        tus.push(emitted_tu);
        endpoint_stats.push(result.endpoint_stats()[old_idx]);
    }

    PreparedOutputs {
        tus,
        endpoint_stats,
        id_mappings,
    }
}
