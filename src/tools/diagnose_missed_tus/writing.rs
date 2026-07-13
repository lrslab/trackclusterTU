use std::collections::HashMap;
use std::path::Path;

use anyhow::{Context, Result};

use crate::io::bed::Bed6Record;
use crate::io::delimited::{DelimitedWriter, Delimiter};
use crate::score::{score1_interval, score2_interval};
use crate::tools::membership::{
    write_membership, MembershipRecord, MembershipStatus, V2MembershipFields, WeightedTu,
    V2_COLUMNS,
};
use crate::tu::ReadRecord;

use super::detection::{
    five_prime_coord, five_prime_coord_interval, three_prime_coord, three_prime_coord_interval,
    CandidateMode,
};
use super::parsing::ExistingTu;
use super::rescue::{assignment_target, AssignmentOrigin, FinalTu, MembershipSchema, RescuePlan};

const CANDIDATE_REPORT_COLUMNS: [&str; 18] = [
    "contig",
    "strand",
    "candidate_start",
    "candidate_end",
    "three_prime_coord",
    "five_prime_coord",
    "family_support",
    "mode_support",
    "mode_fraction",
    "reason",
    "nearest_tu_id",
    "nearest_tu_start",
    "nearest_tu_end",
    "nearest_tu_three_prime_delta_bp",
    "nearest_tu_five_prime_delta_bp",
    "nearest_tu_score1",
    "nearest_tu_score2",
    "gene_hint",
];

fn cmp_candidate_for_report(a: &CandidateMode, b: &CandidateMode) -> std::cmp::Ordering {
    a.contig
        .cmp(&b.contig)
        .then_with(|| a.strand.cmp(&b.strand))
        .then_with(|| a.interval.start().cmp(&b.interval.start()))
        .then_with(|| a.interval.end().cmp(&b.interval.end()))
        .then_with(|| b.mode_support.cmp(&a.mode_support))
}

fn assign_report_ids_in_order(candidates: &mut [CandidateMode], prefix: &str) {
    let width = 4usize.max(candidates.len().to_string().len());
    for (idx, candidate) in candidates.iter_mut().enumerate() {
        candidate.report_id = format!("{prefix}{:0width$}", idx + 1, width = width);
    }
}

pub(super) fn assign_candidate_ids(candidates: &mut [CandidateMode], prefix: &str) {
    candidates.sort_by(cmp_candidate_for_report);
    assign_report_ids_in_order(candidates, prefix);
}

pub(super) fn assign_rescued_candidate_ids(
    candidates: &mut Vec<CandidateMode>,
    rescue_ids: &mut Vec<String>,
    prefix: &str,
) -> Result<()> {
    anyhow::ensure!(
        candidates.len() == rescue_ids.len(),
        "internal rescue candidate and TU ID vectors are not aligned"
    );
    let mut paired: Vec<(CandidateMode, String)> = std::mem::take(candidates)
        .into_iter()
        .zip(std::mem::take(rescue_ids))
        .collect();
    paired.sort_by(|(a, _), (b, _)| cmp_candidate_for_report(a, b));
    for (candidate, rescue_id) in paired {
        candidates.push(candidate);
        rescue_ids.push(rescue_id);
    }
    assign_report_ids_in_order(candidates, prefix);
    Ok(())
}

fn candidate_report_fields(candidate: &CandidateMode) -> Vec<String> {
    let nearest = candidate.nearest_tu.as_ref();
    vec![
        candidate.contig.clone(),
        candidate.strand.as_char().to_string(),
        candidate.interval.start().get().to_string(),
        candidate.interval.end().get().to_string(),
        candidate.three_prime_coord.to_string(),
        candidate.five_prime_coord.to_string(),
        candidate.family_support.to_string(),
        candidate.mode_support.to_string(),
        format!("{:.6}", candidate.mode_fraction),
        candidate.reason.clone(),
        nearest
            .map(|tu| tu.id.clone())
            .unwrap_or_else(|| ".".to_owned()),
        nearest
            .map(|tu| tu.interval.start().get().to_string())
            .unwrap_or_else(|| ".".to_owned()),
        nearest
            .map(|tu| tu.interval.end().get().to_string())
            .unwrap_or_else(|| ".".to_owned()),
        nearest
            .map(|tu| tu.three_prime_delta_bp.to_string())
            .unwrap_or_else(|| ".".to_owned()),
        nearest
            .map(|tu| tu.five_prime_delta_bp.to_string())
            .unwrap_or_else(|| ".".to_owned()),
        nearest
            .map(|tu| format!("{:.6}", tu.score1))
            .unwrap_or_else(|| ".".to_owned()),
        nearest
            .map(|tu| format!("{:.6}", tu.score2))
            .unwrap_or_else(|| ".".to_owned()),
        candidate.gene_hint.clone(),
    ]
}

pub(super) fn write_tsv_report(path: &Path, candidates: &[CandidateMode]) -> Result<()> {
    let mut writer = DelimitedWriter::create(path, Delimiter::Tab, &[])?;
    let mut header = vec!["candidate_id"];
    header.extend(CANDIDATE_REPORT_COLUMNS);
    writer.write_record(header)?;

    for candidate in candidates {
        let mut row = vec![candidate.report_id.clone()];
        row.extend(candidate_report_fields(candidate));
        writer.write_record(row)?;
    }

    writer.flush()?;
    Ok(())
}

pub(super) fn write_rescued_tsv_report(
    path: &Path,
    candidates: &[CandidateMode],
    rescue_ids: &[String],
    final_tus: &[FinalTu],
) -> Result<()> {
    anyhow::ensure!(
        candidates.len() == rescue_ids.len(),
        "internal rescue candidate and TU ID vectors are not aligned"
    );
    let hard_counts: HashMap<&str, u64> = final_tus
        .iter()
        .map(|tu| (tu.id.as_str(), tu.count))
        .collect();
    let mut writer = DelimitedWriter::create(path, Delimiter::Tab, &[])?;
    let mut header = vec!["candidate_id", "rescued_tu_id", "rescued_hard_count"];
    header.extend(CANDIDATE_REPORT_COLUMNS);
    writer.write_record(header)?;

    for (candidate, rescue_id) in candidates.iter().zip(rescue_ids) {
        let rescued_hard_count = hard_counts.get(rescue_id.as_str()).with_context(|| {
            format!("rescued TU {rescue_id:?} is missing from the final TU count set")
        })?;
        let mut row = vec![
            candidate.report_id.clone(),
            rescue_id.clone(),
            rescued_hard_count.to_string(),
        ];
        row.extend(candidate_report_fields(candidate));
        writer.write_record(row)?;
    }

    writer.flush()?;
    Ok(())
}

fn write_bed_report_with_rescue_ids(
    path: &Path,
    candidates: &[CandidateMode],
    rescue_ids: Option<&[String]>,
) -> Result<()> {
    if let Some(rescue_ids) = rescue_ids {
        anyhow::ensure!(
            candidates.len() == rescue_ids.len(),
            "internal rescue candidate and TU ID vectors are not aligned"
        );
    }
    let records: Vec<Bed6Record> = candidates
        .iter()
        .enumerate()
        .map(|(index, candidate)| Bed6Record {
            chrom: candidate.contig.clone(),
            start: candidate.interval.start(),
            end: candidate.interval.end(),
            name: match rescue_ids {
                Some(rescue_ids) => format!(
                    "{}|{}|{}|{}",
                    candidate.report_id, rescue_ids[index], candidate.reason, candidate.gene_hint
                ),
                None => format!(
                    "{}|{}|{}",
                    candidate.report_id, candidate.reason, candidate.gene_hint
                ),
            },
            score: candidate.mode_support.min(u32::MAX as usize) as u32,
            strand: candidate.strand,
            extra_fields: Vec::new(),
        })
        .collect();
    crate::io::bed::write_bed6(path, records.iter())
        .with_context(|| format!("failed to write BED report {path:?}"))?;
    Ok(())
}

pub(super) fn write_bed_report(path: &Path, candidates: &[CandidateMode]) -> Result<()> {
    write_bed_report_with_rescue_ids(path, candidates, None)
}

pub(super) fn write_rescued_bed_report(
    path: &Path,
    candidates: &[CandidateMode],
    rescue_ids: &[String],
) -> Result<()> {
    write_bed_report_with_rescue_ids(path, candidates, Some(rescue_ids))
}

pub(super) fn write_rescued_tu_bed(path: &Path, final_tus: &[FinalTu]) -> Result<()> {
    let records: Vec<Bed6Record> = final_tus
        .iter()
        .map(|tu| Bed6Record {
            chrom: tu.contig.clone(),
            start: tu.interval.start(),
            end: tu.interval.end(),
            name: tu.id.clone(),
            score: 0,
            strand: tu.strand,
            extra_fields: Vec::new(),
        })
        .collect();
    crate::io::bed::write_bed6(path, records.iter())
        .with_context(|| format!("failed to write rescued TU BED {path:?}"))?;
    Ok(())
}

pub(super) fn write_rescued_membership(
    path: &Path,
    reads: &[ReadRecord],
    existing_memberships: &[Option<MembershipRecord>],
    membership_schema: MembershipSchema,
    source_metadata: &[String],
    existing_tus: &[ExistingTu],
    plan: &RescuePlan,
) -> Result<()> {
    anyhow::ensure!(
        plan.hard_assignments.len() == reads.len() && existing_memberships.len() == reads.len(),
        "internal rescue membership vectors are not aligned with the retained reads"
    );
    let mut read_indices: Vec<usize> = (0..reads.len()).collect();
    read_indices.sort_unstable_by(|&a, &b| cmp_read(&reads[a], &reads[b]).then_with(|| a.cmp(&b)));

    let mut records = Vec::with_capacity(reads.len());
    for read_idx in read_indices {
        let original = existing_memberships[read_idx].as_ref();
        let record = match plan.hard_assignments[read_idx] {
            Some(AssignmentOrigin::Rescued(_)) => {
                let assignment = plan.hard_assignments[read_idx]
                    .expect("rescued assignment matched immediately above");
                let (tu_id, interval) = assignment_target(
                    assignment,
                    existing_tus,
                    &plan.candidates,
                    &plan.rescue_ids,
                )?;
                rescued_membership_record(
                    &reads[read_idx],
                    tu_id,
                    interval,
                    original,
                    membership_schema,
                )
            }
            Some(AssignmentOrigin::Existing(_)) | None => {
                let Some(original) = original else {
                    continue;
                };
                original.clone()
            }
        };
        records.push(record);
    }

    let metadata = rescued_membership_metadata(&records, membership_schema, source_metadata);
    let mut writer = DelimitedWriter::create_flexible(path, Delimiter::Tab, &metadata)?;
    for record in records {
        write_membership(&mut writer, &record)?;
    }
    writer.flush()?;
    Ok(())
}

fn rescued_membership_metadata(
    records: &[MembershipRecord],
    source_schema: MembershipSchema,
    source_metadata: &[String],
) -> Vec<String> {
    const LEGACY_COLUMNS: [&str; 4] = ["read_id", "tu_id", "score1", "score2"];
    const COMPATIBLE_V2_COLUMNS: [&str; 6] = [
        "read_id",
        "tu_id",
        "score1",
        "score2",
        "schema_version",
        "assignment_status",
    ];

    let mut contains_legacy = false;
    let mut contains_compatible_v2 = false;
    let mut contains_canonical_v2 = false;
    for record in records {
        match record.v2.as_ref() {
            None => contains_legacy = true,
            Some(v2) if v2.complete => contains_canonical_v2 = true,
            Some(_) => contains_compatible_v2 = true,
        }
    }
    let output_is_v2 = contains_compatible_v2
        || contains_canonical_v2
        || (records.is_empty() && source_schema != MembershipSchema::Legacy);

    let preserved_metadata = source_metadata
        .iter()
        .filter(|line| !is_structural_membership_metadata(line) && !line.starts_with("#rescue_"));
    if !output_is_v2 {
        return preserved_metadata.cloned().collect();
    }

    let mut metadata = vec!["#trackclustertu_membership_schema=v2".to_owned()];
    let distinct_widths = usize::from(contains_legacy)
        + usize::from(contains_compatible_v2)
        + usize::from(contains_canonical_v2);
    if distinct_widths <= 1 {
        let columns = if contains_legacy {
            &LEGACY_COLUMNS[..]
        } else if contains_compatible_v2 {
            &COMPATIBLE_V2_COLUMNS[..]
        } else {
            // Canonical v2 is also the declared shape for an empty v2 table.
            &V2_COLUMNS[..]
        };
        metadata.push(format!("#columns={}", columns.join("\t")));
    } else {
        let mut widths = Vec::with_capacity(3);
        if contains_legacy {
            widths.push(LEGACY_COLUMNS.len().to_string());
            metadata.push(format!("#legacy_columns={}", LEGACY_COLUMNS.join("\t")));
        }
        if contains_compatible_v2 {
            widths.push(COMPATIBLE_V2_COLUMNS.len().to_string());
            metadata.push(format!(
                "#compatible_v2_columns={}",
                COMPATIBLE_V2_COLUMNS.join("\t")
            ));
        }
        if contains_canonical_v2 {
            widths.push(V2_COLUMNS.len().to_string());
            metadata.push(format!("#canonical_v2_columns={}", V2_COLUMNS.join("\t")));
        }
        metadata.insert(1, format!("#row_widths={}", widths.join(",")));
    }
    if contains_legacy {
        metadata.push("#contains_legacy_rows=true".to_owned());
    }
    if contains_compatible_v2 {
        metadata.push("#contains_compatible_v2_rows=true".to_owned());
    }
    metadata.extend(preserved_metadata.cloned());
    metadata.extend([
        "#rescue_unmodified_rows=preserved_from_input".to_owned(),
        "#rescue_promoted_rows=canonical_unique_hard_assignment".to_owned(),
        "#rescue_count_semantics=hard_assignments".to_owned(),
    ]);
    metadata
}

fn is_structural_membership_metadata(line: &str) -> bool {
    [
        "#trackclustertu_membership_schema=",
        "#columns=",
        "#row_widths=",
        "#legacy_columns=",
        "#compatible_v2_columns=",
        "#canonical_v2_columns=",
        "#contains_legacy_rows=",
        "#contains_compatible_v2_rows=",
    ]
    .iter()
    .any(|prefix| line.starts_with(prefix))
}

fn rescued_membership_record(
    read: &ReadRecord,
    tu_id: &str,
    interval: crate::model::Interval,
    original: Option<&MembershipRecord>,
    table_schema: MembershipSchema,
) -> MembershipRecord {
    let score1 = score1_interval(read.interval, interval);
    let score2 = score2_interval(read.interval, interval);
    let use_v2 = original.is_some_and(|record| record.v2.is_some())
        || (original.is_none() && table_schema != MembershipSchema::Legacy);
    if !use_v2 {
        return MembershipRecord::legacy(
            read.id.clone(),
            tu_id.to_owned(),
            format!("{score1:.6}"),
            format!("{score2:.6}"),
        );
    }

    let full_length_evidence = original
        .and_then(|record| record.v2.as_ref())
        .is_some_and(|v2| v2.full_length_evidence);
    MembershipRecord {
        read_id: read.id.clone(),
        hard_tu_id: Some(tu_id.to_owned()),
        score1: format!("{score1:.6}"),
        score2: format!("{score2:.6}"),
        status: MembershipStatus::Unique,
        fractional_assignments: vec![WeightedTu {
            tu_id: tu_id.to_owned(),
            weight: 1.0,
        }],
        v2: Some(V2MembershipFields {
            complete: true,
            best_candidate_tu_id: Some(tu_id.to_owned()),
            second_candidate_tu_id: None,
            second_score1: ".".to_owned(),
            second_score2: ".".to_owned(),
            best_five_prime_delta_bp: five_prime_coord(read)
                .abs_diff(five_prime_coord_interval(interval, read.strand))
                .to_string(),
            best_three_prime_delta_bp: three_prime_coord(read)
                .abs_diff(three_prime_coord_interval(interval, read.strand))
                .to_string(),
            second_five_prime_delta_bp: ".".to_owned(),
            second_three_prime_delta_bp: ".".to_owned(),
            best_assignment_score: format!("{score2:.6}"),
            second_assignment_score: ".".to_owned(),
            score_margin: ".".to_owned(),
            primary_weight: "1".to_owned(),
            full_length_evidence,
        }),
        line_number: 0,
    }
}

pub(super) fn write_rescued_counts_csv(path: &Path, final_tus: &[FinalTu]) -> Result<()> {
    let mut writer = DelimitedWriter::create(path, Delimiter::Comma, &[])?;
    writer.write_record(["tu_id", "count"])?;
    for tu in final_tus {
        writer.write_record([tu.id.clone(), tu.count.to_string()])?;
    }
    writer.flush()?;
    Ok(())
}

fn cmp_read(a: &ReadRecord, b: &ReadRecord) -> std::cmp::Ordering {
    a.contig
        .cmp(&b.contig)
        .then_with(|| a.strand.cmp(&b.strand))
        .then_with(|| a.interval.start().cmp(&b.interval.start()))
        .then_with(|| a.interval.end().cmp(&b.interval.end()))
        .then_with(|| a.id.cmp(&b.id))
}
