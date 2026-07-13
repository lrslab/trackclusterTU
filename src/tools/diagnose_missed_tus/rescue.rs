use std::cmp::Ordering;
use std::collections::{HashMap, HashSet};

use thiserror::Error;

use crate::model::Interval;
use crate::tu::ReadRecord;

use super::detection::{
    five_prime_coord, five_prime_coord_interval, three_prime_coord, three_prime_coord_interval,
    CandidateMode, DetectionParams,
};
use super::parsing::{ExistingMembershipRow, ExistingTu};

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(super) enum AssignmentOrigin {
    Existing(usize),
    Rescued(usize),
}

#[derive(Clone, Debug)]
pub(super) struct FinalTu {
    pub(super) id: String,
    pub(super) contig: String,
    pub(super) strand: crate::model::Strand,
    pub(super) interval: Interval,
    pub(super) count: u64,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
struct BoundaryQuality {
    three_prime_delta_bp: u32,
    five_prime_delta_bp: u32,
}

#[derive(Debug)]
pub(super) struct RescuePlan {
    pub(super) candidates: Vec<CandidateMode>,
    pub(super) hard_assignments: Vec<Option<AssignmentOrigin>>,
    pub(super) rescue_ids: Vec<String>,
    pub(super) final_tus: Vec<FinalTu>,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(super) enum MembershipSchema {
    Legacy,
    V2,
    Mixed,
}

#[derive(Debug)]
pub(super) struct ValidatedExistingMemberships {
    pub(super) read_to_existing: Vec<Option<usize>>,
    pub(super) rows_by_read: Vec<Option<ExistingMembershipRow>>,
    pub(super) schema: MembershipSchema,
}

#[derive(Clone, Debug, Error, PartialEq, Eq)]
pub(super) enum MembershipValidationError {
    #[error("{line}: membership read ID {read_id:?} does not exist in the input BED")]
    UnknownRead { line: usize, read_id: String },
    #[error(
        "{line}: membership TU ID {tu_id:?} referenced by read {read_id:?} does not exist in the existing TU BED"
    )]
    UnknownTu {
        line: usize,
        read_id: String,
        tu_id: String,
    },
}

#[derive(Clone, Debug, Error, PartialEq, Eq)]
pub(super) enum RescueError {
    #[error(
        "unable to allocate {count} unique rescue IDs with prefix {prefix:?}: numeric suffix exhausted"
    )]
    RescueIdSpaceExhausted { count: usize, prefix: String },
    #[error("internal error: existing TU assignment index {index} is invalid")]
    InvalidExistingIndex { index: usize },
    #[error("internal error: rescued candidate assignment index {index} is invalid")]
    InvalidCandidateIndex { index: usize },
    #[error("internal error: rescue ID for candidate index {index} is missing")]
    MissingRescueId { index: usize },
    #[error(
        "internal rescue consistency error: TU ID {tu_id:?} would occur more than once in the rescued TU BED"
    )]
    DuplicateFinalTu { tu_id: String },
    #[error(
        "internal rescue consistency error: membership TU ID {tu_id:?} is absent from the rescued TU BED"
    )]
    MissingFinalTu { tu_id: String },
    #[error(
        "internal rescue consistency error: TU hard counts sum to {tu_count}, but membership contains {hard_assignment_count} hard assignments"
    )]
    CountSumMismatch {
        tu_count: u64,
        hard_assignment_count: u64,
    },
    #[error(
        "internal rescue consistency error: TU {tu_id:?} has count {tu_count}, but {assignment_count} hard assignments"
    )]
    TuCountMismatch {
        tu_id: String,
        tu_count: u64,
        assignment_count: u64,
    },
}

pub(super) fn validate_existing_memberships(
    rows: Vec<ExistingMembershipRow>,
    source_metadata: &[String],
    input_read_ids: &HashSet<String>,
    reads: &[ReadRecord],
    existing_tus: &[ExistingTu],
) -> Result<ValidatedExistingMemberships, MembershipValidationError> {
    let read_lookup: HashMap<&str, usize> = reads
        .iter()
        .enumerate()
        .map(|(idx, read)| (read.id.as_str(), idx))
        .collect();
    let tu_lookup: HashMap<&str, usize> = existing_tus
        .iter()
        .enumerate()
        .map(|(idx, tu)| (tu.id.as_str(), idx))
        .collect();
    let mut read_to_existing = vec![None; reads.len()];
    let mut rows_by_read = vec![None; reads.len()];
    let mut saw_legacy = false;
    let mut saw_v2 = source_metadata.iter().any(|line| {
        line.strip_prefix("#trackclustertu_membership_schema=")
            .is_some_and(|value| value.trim() == "v2")
    });

    for row in rows {
        if row.v2.is_some() {
            saw_v2 = true;
        } else {
            saw_legacy = true;
        }
        if !input_read_ids.contains(&row.read_id) {
            return Err(MembershipValidationError::UnknownRead {
                line: row.line_number,
                read_id: row.read_id.clone(),
            });
        }
        for tu_id in row.candidate_tu_ids() {
            if !tu_lookup.contains_key(tu_id) {
                return Err(MembershipValidationError::UnknownTu {
                    line: row.line_number,
                    read_id: row.read_id.clone(),
                    tu_id: tu_id.to_owned(),
                });
            }
        }
        let hard_tu_idx = row.hard_tu_id.as_deref().map(|tu_id| tu_lookup[tu_id]);
        if let Some(&read_idx) = read_lookup.get(row.read_id.as_str()) {
            read_to_existing[read_idx] = hard_tu_idx;
            rows_by_read[read_idx] = Some(row);
        }
    }
    let schema = match (saw_legacy, saw_v2) {
        (false, true) => MembershipSchema::V2,
        (true, true) => MembershipSchema::Mixed,
        _ => MembershipSchema::Legacy,
    };
    Ok(ValidatedExistingMemberships {
        read_to_existing,
        rows_by_read,
        schema,
    })
}

pub(super) fn build_rescue_plan(
    reads: &[ReadRecord],
    existing_tus: &[ExistingTu],
    candidate_modes: Vec<CandidateMode>,
    read_to_existing: Vec<Option<usize>>,
    params: DetectionParams,
    rescue_prefix: &str,
) -> Result<RescuePlan, RescueError> {
    let candidates = deduplicate_candidates(candidate_modes, params);
    let mut best_candidate_for_read: Vec<Option<(usize, BoundaryQuality)>> =
        vec![None; reads.len()];
    for (candidate_idx, candidate) in candidates.iter().enumerate() {
        for &read_idx in &candidate.support_read_indices {
            let quality = boundary_quality_for_read(&reads[read_idx], candidate.interval);
            let replace = match best_candidate_for_read[read_idx] {
                None => true,
                Some((current_idx, current_quality)) => better_candidate_choice(
                    candidate_idx,
                    quality,
                    candidate.mode_support,
                    current_idx,
                    current_quality,
                    candidates[current_idx].mode_support,
                ),
            };
            if replace {
                best_candidate_for_read[read_idx] = Some((candidate_idx, quality));
            }
        }
    }

    let mut provisional_candidate_counts = vec![0usize; candidates.len()];
    for (candidate_idx, _) in best_candidate_for_read.iter().flatten() {
        provisional_candidate_counts[*candidate_idx] += 1;
    }
    let candidate_kept: Vec<bool> = provisional_candidate_counts
        .iter()
        .map(|&count| count >= params.min_mode_support)
        .collect();

    let mut final_read_to_candidate: Vec<Option<usize>> = vec![None; reads.len()];
    let mut final_candidate_counts = vec![0usize; candidates.len()];
    for (read_idx, read) in reads.iter().enumerate() {
        let Some((candidate_idx, candidate_quality)) = best_candidate_for_read[read_idx] else {
            continue;
        };
        if !candidate_kept[candidate_idx] {
            continue;
        }

        let current_tu = read_to_existing[read_idx].map(|tu_idx| &existing_tus[tu_idx]);
        let use_candidate = match current_tu {
            Some(tu) => {
                let current_quality = boundary_quality_for_read(read, tu.interval);
                better_boundary_quality(
                    candidate_quality,
                    candidates[candidate_idx].mode_support,
                    current_quality,
                    0,
                )
            }
            None => true,
        };

        if use_candidate {
            final_read_to_candidate[read_idx] = Some(candidate_idx);
            final_candidate_counts[candidate_idx] += 1;
        }
    }

    let final_candidate_kept: Vec<bool> = final_candidate_counts
        .iter()
        .map(|&count| count >= params.min_mode_support)
        .collect();
    let kept_candidate_index_map = build_kept_candidate_index_map(&final_candidate_kept);
    let mut candidates: Vec<CandidateMode> = candidates
        .into_iter()
        .enumerate()
        .filter_map(|(candidate_idx, candidate)| {
            final_candidate_kept[candidate_idx].then_some(candidate)
        })
        .collect();

    let reserved_tu_ids: HashSet<String> = existing_tus.iter().map(|tu| tu.id.clone()).collect();
    let rescue_ids = allocate_rescue_ids(rescue_prefix, candidates.len(), &reserved_tu_ids)?;

    let mut final_assignments: Vec<Option<AssignmentOrigin>> = read_to_existing
        .iter()
        .map(|tu_idx| tu_idx.map(AssignmentOrigin::Existing))
        .collect();
    for read_idx in 0..reads.len() {
        let Some(original_candidate_idx) = final_read_to_candidate[read_idx] else {
            continue;
        };
        let Some(new_candidate_idx) = kept_candidate_index_map[original_candidate_idx] else {
            continue;
        };
        final_assignments[read_idx] = Some(AssignmentOrigin::Rescued(new_candidate_idx));
    }

    let mut existing_final_counts: Vec<u64> = vec![0; existing_tus.len()];
    let mut rescue_final_counts: Vec<u64> = vec![0; candidates.len()];
    for assignment in &final_assignments {
        match assignment {
            Some(AssignmentOrigin::Existing(tu_idx)) => existing_final_counts[*tu_idx] += 1,
            Some(AssignmentOrigin::Rescued(candidate_idx)) => {
                rescue_final_counts[*candidate_idx] += 1;
            }
            None => {}
        }
    }

    let mut final_tus: Vec<FinalTu> = Vec::new();
    for (tu_idx, tu) in existing_tus.iter().enumerate() {
        let count = existing_final_counts[tu_idx];
        final_tus.push(FinalTu {
            id: tu.id.clone(),
            contig: tu.contig.clone(),
            strand: tu.strand,
            interval: tu.interval,
            count,
        });
    }
    for (candidate_idx, candidate) in candidates.iter().enumerate() {
        let count = rescue_final_counts[candidate_idx];
        if count > 0 {
            final_tus.push(FinalTu {
                id: rescue_ids[candidate_idx].clone(),
                contig: candidate.contig.clone(),
                strand: candidate.strand,
                interval: candidate.interval,
                count,
            });
        }
    }
    final_tus.sort_by(cmp_final_tu);

    validate_final_rescue_state(
        &final_assignments,
        &final_tus,
        existing_tus,
        &candidates,
        &rescue_ids,
    )?;

    Ok(RescuePlan {
        candidates: std::mem::take(&mut candidates),
        hard_assignments: final_assignments,
        rescue_ids,
        final_tus,
    })
}

fn deduplicate_candidates(
    mut candidates: Vec<CandidateMode>,
    params: DetectionParams,
) -> Vec<CandidateMode> {
    candidates.sort_by(|a, b| {
        b.mode_support
            .cmp(&a.mode_support)
            .then_with(|| b.mode_fraction.total_cmp(&a.mode_fraction))
            .then_with(|| a.contig.cmp(&b.contig))
            .then_with(|| a.strand.cmp(&b.strand))
            .then_with(|| a.interval.start().cmp(&b.interval.start()))
            .then_with(|| a.interval.end().cmp(&b.interval.end()))
    });

    let mut kept: Vec<CandidateMode> = Vec::new();
    'candidate: for candidate in candidates {
        for existing in &kept {
            if existing.contig != candidate.contig || existing.strand != candidate.strand {
                continue;
            }
            let existing_quality = boundary_quality_for_interval(
                candidate.interval,
                existing.interval,
                candidate.strand,
            );
            if existing_quality.three_prime_delta_bp <= params.three_prime_window_bp
                && existing_quality.five_prime_delta_bp <= params.five_prime_window_bp
            {
                continue 'candidate;
            }
        }
        kept.push(candidate);
    }

    kept.sort_by(|a, b| {
        a.contig
            .cmp(&b.contig)
            .then_with(|| a.strand.cmp(&b.strand))
            .then_with(|| a.interval.start().cmp(&b.interval.start()))
            .then_with(|| a.interval.end().cmp(&b.interval.end()))
            .then_with(|| b.mode_support.cmp(&a.mode_support))
    });
    kept
}

fn allocate_rescue_ids(
    prefix: &str,
    count: usize,
    reserved_tu_ids: &HashSet<String>,
) -> Result<Vec<String>, RescueError> {
    let width = 4usize.max(count.to_string().len());
    let mut used_ids = reserved_tu_ids.clone();
    let mut rescue_ids = Vec::with_capacity(count);
    let mut serial = 1usize;

    while rescue_ids.len() < count {
        let candidate_id = format!("{prefix}{serial:0width$}");
        if used_ids.insert(candidate_id.clone()) {
            rescue_ids.push(candidate_id);
        }
        if rescue_ids.len() < count {
            serial = serial
                .checked_add(1)
                .ok_or_else(|| RescueError::RescueIdSpaceExhausted {
                    count,
                    prefix: prefix.to_owned(),
                })?;
        }
    }

    Ok(rescue_ids)
}

fn boundary_quality_for_read(read: &ReadRecord, interval: Interval) -> BoundaryQuality {
    BoundaryQuality {
        three_prime_delta_bp: three_prime_coord(read)
            .abs_diff(three_prime_coord_interval(interval, read.strand)),
        five_prime_delta_bp: five_prime_coord(read)
            .abs_diff(five_prime_coord_interval(interval, read.strand)),
    }
}

fn boundary_quality_for_interval(
    a: Interval,
    b: Interval,
    strand: crate::model::Strand,
) -> BoundaryQuality {
    BoundaryQuality {
        three_prime_delta_bp: three_prime_coord_interval(a, strand)
            .abs_diff(three_prime_coord_interval(b, strand)),
        five_prime_delta_bp: five_prime_coord_interval(a, strand)
            .abs_diff(five_prime_coord_interval(b, strand)),
    }
}

fn better_boundary_quality(
    candidate_quality: BoundaryQuality,
    candidate_support: usize,
    current_quality: BoundaryQuality,
    current_support: usize,
) -> bool {
    (
        candidate_quality.three_prime_delta_bp,
        candidate_quality.five_prime_delta_bp,
    ) < (
        current_quality.three_prime_delta_bp,
        current_quality.five_prime_delta_bp,
    ) || ((
        candidate_quality.three_prime_delta_bp,
        candidate_quality.five_prime_delta_bp,
    ) == (
        current_quality.three_prime_delta_bp,
        current_quality.five_prime_delta_bp,
    ) && candidate_support > current_support)
}

fn better_candidate_choice(
    candidate_idx: usize,
    candidate_quality: BoundaryQuality,
    candidate_support: usize,
    current_idx: usize,
    current_quality: BoundaryQuality,
    current_support: usize,
) -> bool {
    better_boundary_quality(
        candidate_quality,
        candidate_support,
        current_quality,
        current_support,
    ) || (candidate_quality == current_quality
        && candidate_support == current_support
        && candidate_idx < current_idx)
}

fn build_kept_candidate_index_map(final_candidate_kept: &[bool]) -> Vec<Option<usize>> {
    let mut kept_index_map = vec![None; final_candidate_kept.len()];
    let mut next_idx = 0usize;
    for (idx, kept) in final_candidate_kept.iter().copied().enumerate() {
        if kept {
            kept_index_map[idx] = Some(next_idx);
            next_idx += 1;
        }
    }
    kept_index_map
}

fn cmp_final_tu(a: &FinalTu, b: &FinalTu) -> Ordering {
    a.contig
        .cmp(&b.contig)
        .then_with(|| a.strand.cmp(&b.strand))
        .then_with(|| a.interval.start().cmp(&b.interval.start()))
        .then_with(|| a.interval.end().cmp(&b.interval.end()))
        .then_with(|| a.id.cmp(&b.id))
}

pub(super) fn assignment_target<'a>(
    assignment: AssignmentOrigin,
    existing_tus: &'a [ExistingTu],
    candidates: &'a [CandidateMode],
    rescue_ids: &'a [String],
) -> Result<(&'a str, Interval), RescueError> {
    match assignment {
        AssignmentOrigin::Existing(tu_idx) => {
            let tu = existing_tus
                .get(tu_idx)
                .ok_or(RescueError::InvalidExistingIndex { index: tu_idx })?;
            Ok((tu.id.as_str(), tu.interval))
        }
        AssignmentOrigin::Rescued(candidate_idx) => {
            let candidate =
                candidates
                    .get(candidate_idx)
                    .ok_or(RescueError::InvalidCandidateIndex {
                        index: candidate_idx,
                    })?;
            let rescue_id = rescue_ids
                .get(candidate_idx)
                .ok_or(RescueError::MissingRescueId {
                    index: candidate_idx,
                })?;
            Ok((rescue_id.as_str(), candidate.interval))
        }
    }
}

fn validate_final_rescue_state(
    assignments: &[Option<AssignmentOrigin>],
    final_tus: &[FinalTu],
    existing_tus: &[ExistingTu],
    candidates: &[CandidateMode],
    rescue_ids: &[String],
) -> Result<(), RescueError> {
    let mut emitted_counts: HashMap<&str, u64> = HashMap::new();
    for tu in final_tus {
        if emitted_counts.insert(tu.id.as_str(), tu.count).is_some() {
            return Err(RescueError::DuplicateFinalTu {
                tu_id: tu.id.clone(),
            });
        }
    }

    let mut assignment_counts: HashMap<&str, u64> = HashMap::new();
    let mut hard_assignment_count = 0u64;
    for assignment in assignments.iter().flatten().copied() {
        let (tu_id, _) = assignment_target(assignment, existing_tus, candidates, rescue_ids)?;
        if !emitted_counts.contains_key(tu_id) {
            return Err(RescueError::MissingFinalTu {
                tu_id: tu_id.to_owned(),
            });
        }
        *assignment_counts.entry(tu_id).or_default() += 1;
        hard_assignment_count += 1;
    }

    let emitted_count_sum: u64 = final_tus.iter().map(|tu| tu.count).sum();
    if emitted_count_sum != hard_assignment_count {
        return Err(RescueError::CountSumMismatch {
            tu_count: emitted_count_sum,
            hard_assignment_count,
        });
    }
    for (tu_id, emitted_count) in emitted_counts {
        let assignment_count = assignment_counts.get(tu_id).copied().unwrap_or(0);
        if emitted_count != assignment_count {
            return Err(RescueError::TuCountMismatch {
                tu_id: tu_id.to_owned(),
                tu_count: emitted_count,
                assignment_count,
            });
        }
    }

    Ok(())
}
