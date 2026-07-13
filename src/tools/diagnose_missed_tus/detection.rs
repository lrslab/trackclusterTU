use std::cmp::Ordering;
use std::collections::{BTreeMap, HashMap};

use thiserror::Error;

use crate::endpoint::{detect_endpoint_modes, TieDirection};
use crate::model::{Coord, Interval, Strand};
use crate::score::{score1_interval, score2_interval};
use crate::tu::ReadRecord;

use super::parsing::{ExistingTu, GeneRecord};

#[derive(Clone, Copy, Debug)]
pub(super) struct DetectionParams {
    pub(super) three_prime_window_bp: u32,
    pub(super) max_three_prime_family_diameter_bp: u32,
    pub(super) five_prime_window_bp: u32,
    pub(super) min_family_support: usize,
    pub(super) min_mode_support: usize,
    pub(super) min_mode_fraction: f64,
    pub(super) max_candidates_per_family: usize,
}

#[derive(Clone, Copy, Debug, Error, PartialEq)]
pub(super) enum DetectionParamsError {
    #[error("--min-mode-fraction must be between 0 and 1, got {value}")]
    InvalidModeFraction { value: f64 },
    #[error("--max-candidates-per-family must be >= 1")]
    NoCandidatesPerFamily,
}

impl DetectionParams {
    #[allow(clippy::too_many_arguments)]
    pub(super) fn try_new(
        three_prime_window_bp: u32,
        max_three_prime_family_diameter_bp: Option<u32>,
        five_prime_window_bp: u32,
        min_family_support: usize,
        min_mode_support: usize,
        min_mode_fraction: f64,
        max_candidates_per_family: usize,
    ) -> Result<Self, DetectionParamsError> {
        if !(0.0..=1.0).contains(&min_mode_fraction) {
            return Err(DetectionParamsError::InvalidModeFraction {
                value: min_mode_fraction,
            });
        }
        if max_candidates_per_family == 0 {
            return Err(DetectionParamsError::NoCandidatesPerFamily);
        }
        Ok(Self {
            three_prime_window_bp,
            max_three_prime_family_diameter_bp: max_three_prime_family_diameter_bp
                .unwrap_or(three_prime_window_bp),
            five_prime_window_bp,
            min_family_support,
            min_mode_support,
            min_mode_fraction,
            max_candidates_per_family,
        })
    }
}

#[derive(Clone, Debug)]
pub(super) struct BestExistingTu {
    pub(super) id: String,
    pub(super) interval: Interval,
    pub(super) three_prime_delta_bp: u32,
    pub(super) five_prime_delta_bp: u32,
    pub(super) score1: f64,
    pub(super) score2: f64,
}

#[derive(Clone, Debug)]
pub(super) struct CandidateMode {
    pub(super) report_id: String,
    pub(super) contig: String,
    pub(super) strand: Strand,
    pub(super) interval: Interval,
    pub(super) three_prime_coord: u32,
    pub(super) five_prime_coord: u32,
    pub(super) family_support: usize,
    pub(super) mode_support: usize,
    pub(super) mode_fraction: f64,
    pub(super) reason: String,
    pub(super) nearest_tu: Option<BestExistingTu>,
    pub(super) gene_hint: String,
    pub(super) support_read_indices: Vec<usize>,
}

#[derive(Clone, Debug)]
struct ThreePrimeFamily {
    contig: String,
    strand: Strand,
    three_prime_coord: u32,
    read_indices: Vec<usize>,
}

#[derive(Clone, Debug, PartialEq, Eq)]
struct BoundaryPeak {
    five_prime_coord: u32,
    support: usize,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum BoundaryKind {
    ThreePrime,
    FivePrime,
}

pub(super) fn detect_candidate_modes(
    reads: &[ReadRecord],
    existing_tus: &[ExistingTu],
    genes: &[GeneRecord],
    params: DetectionParams,
) -> Vec<CandidateMode> {
    let mut tus_by_partition: HashMap<(String, Strand), Vec<ExistingTu>> = HashMap::new();
    for tu in existing_tus {
        tus_by_partition
            .entry((tu.contig.clone(), tu.strand))
            .or_default()
            .push(tu.clone());
    }

    let mut genes_by_partition: HashMap<(String, Strand), Vec<GeneRecord>> = HashMap::new();
    for gene in genes {
        genes_by_partition
            .entry((gene.contig.clone(), gene.strand))
            .or_default()
            .push(gene.clone());
    }

    let families = build_three_prime_families(
        reads,
        params.three_prime_window_bp,
        params.max_three_prime_family_diameter_bp,
    );
    let mut candidates: Vec<CandidateMode> = Vec::new();

    for family in families {
        let family_support = family.read_indices.len();
        if family_support < params.min_family_support {
            continue;
        }

        let partition_key = (family.contig.clone(), family.strand);
        let partition_tus = tus_by_partition
            .get(&partition_key)
            .map(Vec::as_slice)
            .unwrap_or(&[]);
        let partition_genes = genes_by_partition
            .get(&partition_key)
            .map(Vec::as_slice)
            .unwrap_or(&[]);

        let family_reads: Vec<&ReadRecord> =
            family.read_indices.iter().map(|&idx| &reads[idx]).collect();
        let mut peaks =
            detect_five_prime_modes(&family_reads, family.strand, params.five_prime_window_bp);
        peaks.sort_by(|a, b| {
            b.support.cmp(&a.support).then_with(|| {
                boundary_priority(family.strand, BoundaryKind::FivePrime, a.five_prime_coord)
                    .cmp(&boundary_priority(
                        family.strand,
                        BoundaryKind::FivePrime,
                        b.five_prime_coord,
                    ))
                    .reverse()
            })
        });

        for peak in peaks.into_iter().take(params.max_candidates_per_family) {
            let support_read_indices: Vec<usize> = family
                .read_indices
                .iter()
                .copied()
                .filter(|&read_idx| {
                    let read = &reads[read_idx];
                    three_prime_coord(read).abs_diff(family.three_prime_coord)
                        <= params.three_prime_window_bp
                        && five_prime_coord(read).abs_diff(peak.five_prime_coord)
                            <= params.five_prime_window_bp
                })
                .collect();
            let mode_support = support_read_indices.len();
            if mode_support < params.min_mode_support {
                continue;
            }

            let mode_fraction = mode_support as f64 / family_support as f64;
            if mode_fraction < params.min_mode_fraction {
                continue;
            }

            let Some(interval) = interval_from_boundaries(
                family.strand,
                peak.five_prime_coord,
                family.three_prime_coord,
            ) else {
                continue;
            };

            let best_match = best_existing_tu_match(
                &family.contig,
                family.strand,
                interval,
                partition_tus,
                params.three_prime_window_bp,
                params.five_prime_window_bp,
            );
            if best_match.represented {
                continue;
            }

            candidates.push(CandidateMode {
                report_id: String::new(),
                contig: family.contig.clone(),
                strand: family.strand,
                interval,
                three_prime_coord: family.three_prime_coord,
                five_prime_coord: peak.five_prime_coord,
                family_support,
                mode_support,
                mode_fraction,
                reason: infer_reason(&best_match),
                nearest_tu: best_match.nearest,
                gene_hint: gene_hint(interval, partition_genes),
                support_read_indices,
            });
        }
    }

    candidates
}

fn build_three_prime_families(
    reads: &[ReadRecord],
    window_bp: u32,
    max_diameter_bp: u32,
) -> Vec<ThreePrimeFamily> {
    let mut indices: Vec<usize> = (0..reads.len()).collect();
    indices.sort_unstable_by(|&a, &b| {
        reads[a]
            .contig
            .cmp(&reads[b].contig)
            .then_with(|| reads[a].strand.cmp(&reads[b].strand))
            .then_with(|| three_prime_coord(&reads[a]).cmp(&three_prime_coord(&reads[b])))
            .then_with(|| five_prime_coord(&reads[a]).cmp(&five_prime_coord(&reads[b])))
            .then_with(|| reads[a].id.cmp(&reads[b].id))
    });

    let mut families: Vec<ThreePrimeFamily> = Vec::new();
    let mut family_indices: Vec<usize> = Vec::new();
    let mut family_contig: Option<String> = None;
    let mut family_strand = Strand::Unknown;
    let mut first_three_prime: Option<u32> = None;
    let mut prev_three_prime: Option<u32> = None;

    for read_idx in indices {
        let read = &reads[read_idx];
        let coord = three_prime_coord(read);
        let same_partition = family_contig
            .as_deref()
            .map(|contig| contig == read.contig.as_str() && family_strand == read.strand)
            .unwrap_or(false);
        let same_family = same_partition
            && prev_three_prime
                .map(|prev| coord.saturating_sub(prev) <= window_bp)
                .unwrap_or(false)
            && first_three_prime
                .map(|first| coord.saturating_sub(first) <= max_diameter_bp)
                .unwrap_or(false);

        if !family_indices.is_empty() && !same_family {
            families.push(finalize_three_prime_family(
                reads,
                std::mem::take(&mut family_indices),
            ));
        }

        if family_indices.is_empty() {
            family_contig = Some(read.contig.clone());
            family_strand = read.strand;
            first_three_prime = Some(coord);
        }

        family_indices.push(read_idx);
        prev_three_prime = Some(coord);
    }

    if !family_indices.is_empty() {
        families.push(finalize_three_prime_family(reads, family_indices));
    }

    families
}

fn finalize_three_prime_family(reads: &[ReadRecord], read_indices: Vec<usize>) -> ThreePrimeFamily {
    let first = &reads[read_indices[0]];
    let three_prime_coord = choose_mode_coordinate(
        read_indices
            .iter()
            .map(|&idx| three_prime_coord(&reads[idx])),
        first.strand,
        BoundaryKind::ThreePrime,
    );

    ThreePrimeFamily {
        contig: first.contig.clone(),
        strand: first.strand,
        three_prime_coord,
        read_indices,
    }
}

fn detect_five_prime_modes(
    reads: &[&ReadRecord],
    strand: Strand,
    window_bp: u32,
) -> Vec<BoundaryPeak> {
    let direction = match strand {
        Strand::Plus | Strand::Unknown => TieDirection::PreferLower,
        Strand::Minus => TieDirection::PreferHigher,
    };
    detect_endpoint_modes(
        reads.iter().map(|read| five_prime_coord(read)),
        window_bp,
        direction,
    )
    .into_iter()
    .map(|mode| BoundaryPeak {
        five_prime_coord: mode.coordinate(),
        support: mode.support(),
    })
    .collect()
}

#[derive(Default)]
struct MatchSummary {
    represented: bool,
    any_same_three_prime: bool,
    any_same_five_prime: bool,
    nearest: Option<BestExistingTu>,
}

fn best_existing_tu_match(
    contig: &str,
    strand: Strand,
    candidate: Interval,
    existing_tus: &[ExistingTu],
    three_prime_window_bp: u32,
    five_prime_window_bp: u32,
) -> MatchSummary {
    let candidate_three_prime = three_prime_coord_interval(candidate, strand);
    let candidate_five_prime = five_prime_coord_interval(candidate, strand);
    let mut summary = MatchSummary::default();

    for tu in existing_tus {
        if tu.contig != contig || tu.strand != strand {
            continue;
        }

        let three_prime_delta_bp =
            candidate_three_prime.abs_diff(three_prime_coord_interval(tu.interval, strand));
        let five_prime_delta_bp =
            candidate_five_prime.abs_diff(five_prime_coord_interval(tu.interval, strand));
        let score1 = score1_interval(candidate, tu.interval);
        let score2 = score2_interval(candidate, tu.interval);

        if three_prime_delta_bp <= three_prime_window_bp {
            summary.any_same_three_prime = true;
        }
        if five_prime_delta_bp <= five_prime_window_bp {
            summary.any_same_five_prime = true;
        }
        if three_prime_delta_bp <= three_prime_window_bp
            && five_prime_delta_bp <= five_prime_window_bp
        {
            summary.represented = true;
        }

        let candidate_match = BestExistingTu {
            id: tu.id.clone(),
            interval: tu.interval,
            three_prime_delta_bp,
            five_prime_delta_bp,
            score1,
            score2,
        };

        let replace = match summary.nearest.as_ref() {
            None => true,
            Some(current) => {
                (
                    candidate_match.three_prime_delta_bp,
                    candidate_match.five_prime_delta_bp,
                ) < (current.three_prime_delta_bp, current.five_prime_delta_bp)
                    || ((
                        candidate_match.three_prime_delta_bp,
                        candidate_match.five_prime_delta_bp,
                    ) == (current.three_prime_delta_bp, current.five_prime_delta_bp)
                        && candidate_match.score2.total_cmp(&current.score2) == Ordering::Greater)
            }
        };
        if replace {
            summary.nearest = Some(candidate_match);
        }
    }

    summary
}

fn infer_reason(match_summary: &MatchSummary) -> String {
    if match_summary.any_same_three_prime && !match_summary.any_same_five_prime {
        "same_3p_missing_5p_mode".to_owned()
    } else if match_summary.any_same_five_prime && !match_summary.any_same_three_prime {
        "same_5p_missing_3p_mode".to_owned()
    } else if match_summary.any_same_three_prime && match_summary.any_same_five_prime {
        "boundary_mode_shifted".to_owned()
    } else {
        "no_matching_tu".to_owned()
    }
}

fn gene_hint(candidate: Interval, genes: &[GeneRecord]) -> String {
    let overlaps: Vec<&str> = genes
        .iter()
        .filter(|gene| gene.interval.overlaps(candidate))
        .map(|gene| gene.id.as_str())
        .collect();
    if !overlaps.is_empty() {
        return overlaps.join(",");
    }

    let mut nearest: Option<(&str, u32)> = None;
    for gene in genes {
        let distance = interval_distance(candidate, gene.interval);
        let replace = match nearest {
            None => true,
            Some((_, best_distance)) => distance < best_distance,
        };
        if replace {
            nearest = Some((gene.id.as_str(), distance));
        }
    }

    nearest
        .map(|(id, _)| id.to_owned())
        .unwrap_or_else(|| ".".to_owned())
}

fn interval_distance(a: Interval, b: Interval) -> u32 {
    if a.overlaps(b) {
        0
    } else if a.end() <= b.start() {
        b.start().get().saturating_sub(a.end().get())
    } else {
        a.start().get().saturating_sub(b.end().get())
    }
}

fn interval_from_boundaries(strand: Strand, five_prime: u32, three_prime: u32) -> Option<Interval> {
    let start = match strand {
        Strand::Plus | Strand::Unknown => five_prime,
        Strand::Minus => three_prime,
    };
    let end = match strand {
        Strand::Plus | Strand::Unknown => three_prime,
        Strand::Minus => five_prime,
    };
    Interval::new(Coord::new(start), Coord::new(end)).ok()
}

pub(super) fn three_prime_coord(read: &ReadRecord) -> u32 {
    three_prime_coord_interval(read.interval, read.strand)
}

pub(super) fn five_prime_coord(read: &ReadRecord) -> u32 {
    five_prime_coord_interval(read.interval, read.strand)
}

pub(super) fn three_prime_coord_interval(interval: Interval, strand: Strand) -> u32 {
    match strand {
        Strand::Plus | Strand::Unknown => interval.end().get(),
        Strand::Minus => interval.start().get(),
    }
}

pub(super) fn five_prime_coord_interval(interval: Interval, strand: Strand) -> u32 {
    match strand {
        Strand::Plus | Strand::Unknown => interval.start().get(),
        Strand::Minus => interval.end().get(),
    }
}

fn choose_mode_coordinate<I>(coords: I, strand: Strand, boundary: BoundaryKind) -> u32
where
    I: IntoIterator<Item = u32>,
{
    let mut counts: BTreeMap<u32, usize> = BTreeMap::new();
    for coord in coords {
        *counts.entry(coord).or_default() += 1;
    }

    counts
        .into_iter()
        .max_by(|(coord_a, count_a), (coord_b, count_b)| {
            count_a.cmp(count_b).then_with(|| {
                boundary_priority(strand, boundary, *coord_a)
                    .cmp(&boundary_priority(strand, boundary, *coord_b))
            })
        })
        .map(|(coord, _)| coord)
        .unwrap_or(0)
}

fn boundary_priority(strand: Strand, boundary: BoundaryKind, coord: u32) -> i64 {
    let coord = coord as i64;
    match (strand, boundary) {
        (Strand::Plus, BoundaryKind::ThreePrime) => coord,
        (Strand::Minus, BoundaryKind::ThreePrime) => -coord,
        (Strand::Unknown, BoundaryKind::ThreePrime) => coord,
        (Strand::Plus, BoundaryKind::FivePrime) => -coord,
        (Strand::Minus, BoundaryKind::FivePrime) => coord,
        (Strand::Unknown, BoundaryKind::FivePrime) => -coord,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn plus_read(id: &str, three_prime: u32) -> ReadRecord {
        ReadRecord {
            id: id.to_owned(),
            contig: "chr1".to_owned(),
            strand: Strand::Plus,
            interval: Interval::new(Coord::new(0), Coord::new(three_prime)).unwrap(),
        }
    }

    fn read_with_five_prime(id: String, five_prime: u32, strand: Strand) -> ReadRecord {
        let interval = match strand {
            Strand::Plus => Interval::new(Coord::new(five_prime), Coord::new(1_000)).unwrap(),
            Strand::Minus => Interval::new(Coord::new(0), Coord::new(five_prime)).unwrap(),
            Strand::Unknown => unreachable!(),
        };
        ReadRecord {
            id,
            contig: "chr1".to_owned(),
            strand,
            interval,
        }
    }

    fn detect_five_prime_modes_reference(
        reads: &[&ReadRecord],
        strand: Strand,
        window_bp: u32,
    ) -> Vec<BoundaryPeak> {
        let mut counts: BTreeMap<u32, usize> = BTreeMap::new();
        for read in reads {
            *counts.entry(five_prime_coord(read)).or_default() += 1;
        }
        let mut remaining: Vec<(u32, usize)> = counts.into_iter().collect();
        let mut peaks = Vec::new();
        while !remaining.is_empty() {
            let mut best_idx = 0usize;
            let mut best_support = 0usize;
            let mut best_count = 0usize;
            for i in 0..remaining.len() {
                let center = remaining[i].0;
                let support = remaining
                    .iter()
                    .filter(|(coord, _)| center.abs_diff(*coord) <= window_bp)
                    .map(|(_, count)| *count)
                    .sum::<usize>();
                let count = remaining[i].1;
                if support > best_support
                    || (support == best_support && count > best_count)
                    || (support == best_support
                        && count == best_count
                        && boundary_priority(strand, BoundaryKind::FivePrime, center)
                            > boundary_priority(
                                strand,
                                BoundaryKind::FivePrime,
                                remaining[best_idx].0,
                            ))
                {
                    best_idx = i;
                    best_support = support;
                    best_count = count;
                }
            }
            let center = remaining[best_idx].0;
            let mut support = 0usize;
            remaining.retain(|(coord, count)| {
                if center.abs_diff(*coord) <= window_bp {
                    support += *count;
                    false
                } else {
                    true
                }
            });
            peaks.push(BoundaryPeak {
                five_prime_coord: center,
                support,
            });
        }
        peaks
    }

    #[test]
    fn invalid_detection_configuration_is_typed() {
        assert!(matches!(
            DetectionParams::try_new(12, None, 10, 20, 20, f64::NAN, 3),
            Err(DetectionParamsError::InvalidModeFraction { .. })
        ));
        assert_eq!(
            DetectionParams::try_new(12, None, 10, 20, 20, 0.2, 0).unwrap_err(),
            DetectionParamsError::NoCandidatesPerFamily
        );
    }

    #[test]
    fn three_prime_families_do_not_single_link_across_the_diameter() {
        let reads = vec![
            plus_read("r0", 100),
            plus_read("r1", 112),
            plus_read("r2", 124),
            plus_read("r3", 136),
        ];

        let families = build_three_prime_families(&reads, 12, 12);
        assert_eq!(families.len(), 2);
        assert_eq!(families[0].read_indices.len(), 2);
        assert_eq!(families[1].read_indices.len(), 2);
    }

    #[test]
    fn three_prime_family_diameter_can_be_explicitly_relaxed() {
        let reads = vec![
            plus_read("r0", 100),
            plus_read("r1", 112),
            plus_read("r2", 124),
            plus_read("r3", 136),
        ];

        let families = build_three_prime_families(&reads, 12, 36);
        assert_eq!(families.len(), 1);
        assert_eq!(families[0].read_indices.len(), 4);
    }

    #[test]
    fn indexed_five_prime_modes_match_reference_algorithm() {
        let mut state = 0x9e37_79b9_u32;
        for strand in [Strand::Plus, Strand::Minus] {
            for window_bp in [0, 1, 5, 12, 50] {
                for case_idx in 0..40 {
                    let mut reads = Vec::new();
                    for read_idx in 0..80 {
                        state = state.wrapping_mul(1_664_525).wrapping_add(1_013_904_223);
                        let coordinate = 100 + (state % 200);
                        reads.push(read_with_five_prime(
                            format!("{case_idx}_{read_idx}"),
                            coordinate,
                            strand,
                        ));
                    }
                    let references: Vec<&ReadRecord> = reads.iter().collect();
                    assert_eq!(
                        detect_five_prime_modes(&references, strand, window_bp),
                        detect_five_prime_modes_reference(&references, strand, window_bp)
                    );
                }
            }
        }
    }
}
