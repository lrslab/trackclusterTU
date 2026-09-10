//! Deterministic, strand-aware transcription-unit clustering and read assignment.

pub(crate) mod multi;

use std::cmp::{Ordering, Reverse};
use std::collections::{BTreeMap, BTreeSet, BinaryHeap, HashMap, HashSet};

use rayon::prelude::*;
use thiserror::Error;

use crate::model::{Coord, Interval, Strand};
use crate::score::{score1_interval, score2_interval};

#[derive(Clone, Debug, PartialEq, Eq)]
/// One boundary-bearing molecule supplied to TU clustering.
pub struct ReadRecord {
    /// Reference sequence name.
    pub contig: String,
    /// Molecule alignment strand.
    pub strand: Strand,
    /// Half-open aligned span.
    pub interval: Interval,
    /// Stable molecule identifier.
    pub id: String,
}

#[derive(Clone, Debug, PartialEq, Eq)]
/// One emitted transcription-unit consensus.
pub struct Tu {
    /// TU identifier.
    pub id: String,
    /// Reference sequence name.
    pub contig: String,
    /// TU strand.
    pub strand: Strand,
    /// Half-open consensus span.
    pub interval: Interval,
    /// Index of the example molecule in the original read slice.
    pub rep_read_index: usize,
}

#[derive(Clone, Debug, PartialEq, Eq)]
/// TU calls and index-aligned assignments/statistics from a clustering run.
pub struct TuClusteringResult {
    /// Emitted TUs in deterministic order.
    tus: Vec<Tu>,
    /// TU index for each element of the input read slice.
    read_to_tu: Vec<usize>,
    /// Endpoint consensus and spread statistics aligned one-to-one with [`Self::tus`].
    endpoint_stats: Vec<TuEndpointStats>,
}

impl TuClusteringResult {
    fn new(tus: Vec<Tu>, read_to_tu: Vec<usize>, endpoint_stats: Vec<TuEndpointStats>) -> Self {
        assert_eq!(
            tus.len(),
            endpoint_stats.len(),
            "TU and endpoint-stat vectors must remain aligned"
        );
        assert!(
            read_to_tu.iter().all(|&tu_index| tu_index < tus.len()),
            "every read assignment must reference an emitted TU"
        );
        Self {
            tus,
            read_to_tu,
            endpoint_stats,
        }
    }

    /// Emitted TUs in deterministic order.
    pub fn tus(&self) -> &[Tu] {
        &self.tus
    }

    /// TU index for each input read.
    pub fn read_to_tu(&self) -> &[usize] {
        &self.read_to_tu
    }

    /// Endpoint statistics aligned one-to-one with [`Self::tus`].
    pub fn endpoint_stats(&self) -> &[TuEndpointStats] {
        &self.endpoint_stats
    }

    /// Return the TU assigned to one input-read index, if the index exists.
    pub fn tu_for_read(&self, read_index: usize) -> Option<&Tu> {
        self.read_to_tu
            .get(read_index)
            .and_then(|&tu_index| self.tus.get(tu_index))
    }

    /// Consume the result into its three index-aligned components.
    pub fn into_parts(self) -> (Vec<Tu>, Vec<usize>, Vec<TuEndpointStats>) {
        (self.tus, self.read_to_tu, self.endpoint_stats)
    }
}

/// Auditable endpoint evidence for one emitted transcript unit.
///
/// The consensus coordinates define the emitted TU interval. Min/max coordinates are
/// genomic coordinates even on the minus strand, so each spread is simply `max - min`.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct TuEndpointStats {
    /// Number of reads supporting the TU.
    pub support: usize,
    /// Strand-aware 5-prime consensus coordinate.
    pub five_prime_consensus: u32,
    /// Strand-aware 3-prime consensus coordinate.
    pub three_prime_consensus: u32,
    /// Minimum observed genomic 5-prime coordinate.
    pub five_prime_min: u32,
    /// Maximum observed genomic 5-prime coordinate.
    pub five_prime_max: u32,
    /// Minimum observed genomic 3-prime coordinate.
    pub three_prime_min: u32,
    /// Maximum observed genomic 3-prime coordinate.
    pub three_prime_max: u32,
}

/// How confidently a read can be assigned to an emitted transcript unit.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum AssignmentStatus {
    /// Exactly one candidate is supported outside the configured ambiguity margin.
    Unique,
    /// Two or more candidates have assignment scores within the configured margin.
    Ambiguous,
    /// Exactly one candidate qualifies only through the contained-read attachment rule.
    Partial,
    /// No emitted TU satisfies the direct or contained-read assignment rules.
    Unassigned,
}

impl std::fmt::Display for AssignmentStatus {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.write_str(match self {
            Self::Unique => "unique",
            Self::Ambiguous => "ambiguous",
            Self::Partial => "partial",
            Self::Unassigned => "unassigned",
        })
    }
}

/// Auditable biological evidence for one read-to-TU candidate.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct AssignmentCandidate {
    /// Index into the TU slice passed to [`assign_reads_to_tus`].
    pub tu_index: usize,
    /// Span-Jaccard similarity.
    pub score1: f64,
    /// Overlap-over-longer similarity.
    pub score2: f64,
    /// Absolute strand-aware 5-prime coordinate delta.
    pub five_prime_delta_bp: u32,
    /// Absolute strand-aware 3-prime coordinate delta.
    pub three_prime_delta_bp: u32,
    /// The scalar used for ambiguity margins. This is currently overlap-over-longer
    /// (overlap divided by the longer span).
    pub assignment_score: f64,
    /// True when the candidate qualifies only through the contained-read rule.
    pub partial: bool,
}

/// One exact fractional contribution made by a read.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct FractionalAssignment {
    /// Index into the candidate TU slice.
    pub tu_index: usize,
    /// Exact fractional molecule contribution.
    pub weight: f64,
}

/// Best/runner-up evidence and any count contributions for one read.
#[derive(Clone, Debug, PartialEq)]
pub struct ReadAssignment {
    /// Assignment confidence classification.
    pub status: AssignmentStatus,
    /// Highest-ranked qualifying candidate, if any.
    pub best: Option<AssignmentCandidate>,
    /// Second-highest qualifying candidate, if any.
    pub second_best: Option<AssignmentCandidate>,
    /// Non-negative best-minus-runner-up assignment score.
    pub score_margin: Option<f64>,
    /// Count contributions when fractional ambiguity handling is enabled.
    pub fractional_assignments: Vec<FractionalAssignment>,
}

impl TuEndpointStats {
    /// Return the 5-prime genomic coordinate spread in base pairs.
    pub fn five_prime_spread_bp(self) -> u32 {
        self.five_prime_max - self.five_prime_min
    }

    /// Return the 3-prime genomic coordinate spread in base pairs.
    pub fn three_prime_spread_bp(self) -> u32 {
        self.three_prime_max - self.three_prime_min
    }
}

#[derive(Clone, Debug, Default, PartialEq, Eq)]
/// Peak partition/region sizes observed during clustering.
pub struct TuClusteringStats {
    /// Number of reference/strand partitions.
    pub partition_count: usize,
    /// Largest number of reads in one partition.
    pub max_partition_reads: usize,
    /// Number of overlap-connected regions processed.
    pub region_count: usize,
    /// Largest number of reads in one region.
    pub max_region_reads: usize,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
/// Controls for TU attachment, fragment retention, and read assignment.
pub struct TuClusteringOptions {
    /// Enable second-pass pooling, contained-fragment retention, and partial assignment.
    ///
    /// When false, span-Jaccard seed components still undergo direct-consensus refinement.
    pub attach_contained_reads: bool,
    /// Maximum strand-aware 3-prime overhang beyond the parent for a second-pass edge.
    ///
    /// Internal 3-prime ends are allowed. This is also the minimum jitter window for
    /// fragment retention and partial assignment, which scale with read length.
    pub three_prime_tolerance_bp: u32,
    /// Optional maximum absolute 5-prime mismatch override.
    pub max_five_prime_delta_bp: Option<u32>,
}

impl Default for TuClusteringOptions {
    fn default() -> Self {
        Self {
            attach_contained_reads: true,
            three_prime_tolerance_bp: 12,
            max_five_prime_delta_bp: None,
        }
    }
}

#[derive(Error, Debug)]
/// Validation or molecule-domain error returned by TU algorithms.
pub enum TuError {
    /// Span-Jaccard threshold is non-finite or outside `[0, 1]`.
    #[error("span_jaccard_threshold must be between 0 and 1, got {value}")]
    InvalidScore1Threshold {
        /// Invalid threshold.
        value: f64,
    },

    /// Overlap-over-longer threshold is non-finite or outside `[0, 1]`.
    #[error("overlap_over_longer_threshold must be between 0 and 1, got {value}")]
    InvalidScore2Threshold {
        /// Invalid threshold.
        value: f64,
    },

    /// Ambiguity margin is non-finite or outside `[0, 1]`.
    #[error("ambiguity_margin must be between 0 and 1, got {value}")]
    InvalidAmbiguityMargin {
        /// Invalid margin.
        value: f64,
    },

    /// A read lacks the strand required for endpoint-aware comparison.
    #[error(
        "read {read_id:?} has unknown strand; boundary-aware TU clustering requires '+' or '-'"
    )]
    UnknownStrand {
        /// Identifier of the offending read.
        read_id: String,
    },
}

#[derive(Clone, Debug)]
struct Dsu {
    parent: Vec<usize>,
    rank: Vec<u8>,
}

impl Dsu {
    fn new(n: usize) -> Self {
        Self {
            parent: (0..n).collect(),
            rank: vec![0; n],
        }
    }

    fn find(&mut self, mut x: usize) -> usize {
        let mut root = x;
        while self.parent[root] != root {
            root = self.parent[root];
        }

        while self.parent[x] != x {
            let next = self.parent[x];
            self.parent[x] = root;
            x = next;
        }

        root
    }

    fn union(&mut self, a: usize, b: usize) {
        let ra = self.find(a);
        let rb = self.find(b);
        if ra == rb {
            return;
        }

        let rank_a = self.rank[ra];
        let rank_b = self.rank[rb];
        match rank_a.cmp(&rank_b) {
            Ordering::Less => self.parent[ra] = rb,
            Ordering::Greater => self.parent[rb] = ra,
            Ordering::Equal => {
                self.parent[rb] = ra;
                self.rank[ra] = rank_a.saturating_add(1);
            }
        }
    }
}

#[derive(Clone, Debug)]
struct SeedCluster {
    rep_read_index: usize,
    members: Vec<usize>,
}

fn cmp_read_idx(reads: &[ReadRecord], a: usize, b: usize) -> Ordering {
    cmp_read(&reads[a], &reads[b]).then_with(|| a.cmp(&b))
}

fn cmp_read(a: &ReadRecord, b: &ReadRecord) -> Ordering {
    a.contig
        .cmp(&b.contig)
        .then_with(|| a.strand.cmp(&b.strand))
        .then_with(|| a.interval.start().cmp(&b.interval.start()))
        .then_with(|| a.interval.end().cmp(&b.interval.end()))
        .then_with(|| a.id.cmp(&b.id))
}

fn read_len(read: &ReadRecord) -> u32 {
    read.interval.len()
}

fn better_representative(candidate: &ReadRecord, current: &ReadRecord) -> bool {
    let candidate_len = read_len(candidate);
    let current_len = read_len(current);

    if candidate_len != current_len {
        return candidate_len > current_len;
    }

    candidate.id < current.id
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

fn cluster_partition_score1(
    indices_sorted: &[usize],
    reads: &[ReadRecord],
    dsu: &mut Dsu,
    score1_threshold: f64,
) {
    let mut active: Vec<usize> = Vec::new();
    let mut next_active: Vec<usize> = Vec::new();
    let mut prev_unique: Option<(u32, u32, usize)> = None;

    for local_pos in 0..indices_sorted.len() {
        let idx = indices_sorted[local_pos];
        let read = &reads[idx];
        if read.interval.is_empty() {
            continue;
        }

        let current_start = read.interval.start().get();
        let current_end = read.interval.end().get();
        if let Some((prev_start, prev_end, prev_pos)) = prev_unique {
            if current_start == prev_start && current_end == prev_end {
                dsu.union(local_pos, prev_pos);
                continue;
            }
        }
        prev_unique = Some((current_start, current_end, local_pos));

        let current_len_u32 = read_len(read);
        let current_len_u64 = current_len_u32 as u64;
        let current_len_f64 = current_len_u32 as f64;
        next_active.clear();
        next_active.reserve(active.len());

        for &other_pos in &active {
            let other_idx = indices_sorted[other_pos];
            let other = &reads[other_idx];
            if other.interval.is_empty() {
                continue;
            }
            if other.interval.end().get() <= current_start {
                continue;
            }

            let other_len_u32 = read_len(other);
            if other_len_u32 == 0 {
                continue;
            }

            let start_delta = current_start.saturating_sub(other.interval.start().get());
            let max_keep_start_delta = ((other_len_u32 as f64) * (1.0 - score1_threshold))
                .ceil()
                .max(0.0) as u32;
            if start_delta > max_keep_start_delta {
                continue;
            }

            next_active.push(other_pos);

            let len_a = current_len_u64;
            let len_b = other_len_u32 as u64;
            if len_a == 0 || len_b == 0 {
                continue;
            }
            let max_possible = (len_a.min(len_b) as f64) / (len_a.max(len_b) as f64);
            if max_possible < score1_threshold {
                continue;
            }

            let max_pair_start_delta = (((other_len_u32 as f64)
                - (score1_threshold * current_len_f64))
                / (1.0 + score1_threshold))
                .ceil()
                .max(0.0) as u32;
            if start_delta > max_pair_start_delta {
                continue;
            }

            if score1_interval(read.interval, other.interval) >= score1_threshold {
                dsu.union(local_pos, other_pos);
            }
        }

        std::mem::swap(&mut active, &mut next_active);
        active.push(local_pos);
    }
}

fn format_tu_id(index: usize, width: usize) -> String {
    format!("TU{:0width$}", index + 1, width = width)
}

fn resolve_root(cluster_idx: usize, parents: &mut [Option<usize>]) -> usize {
    let mut current = cluster_idx;
    let mut path: Vec<usize> = Vec::new();

    while let Some(parent) = parents[current] {
        path.push(current);
        current = parent;
    }

    for node in path {
        parents[node] = Some(current);
    }

    current
}

fn five_prime_coord(read: &ReadRecord) -> u32 {
    match read.strand {
        Strand::Plus | Strand::Unknown => read.interval.start().get(),
        Strand::Minus => read.interval.end().get(),
    }
}

fn three_prime_coord(read: &ReadRecord) -> u32 {
    match read.strand {
        Strand::Plus | Strand::Unknown => read.interval.end().get(),
        Strand::Minus => read.interval.start().get(),
    }
}

fn endpoint_coord(interval: Interval, strand: Strand, five_prime: bool) -> u32 {
    match (strand, five_prime) {
        (Strand::Plus | Strand::Unknown, true) | (Strand::Minus, false) => interval.start().get(),
        (Strand::Plus | Strand::Unknown, false) | (Strand::Minus, true) => interval.end().get(),
    }
}

fn cmp_in_transcription_direction(strand: Strand, a: u32, b: u32) -> Ordering {
    match strand {
        Strand::Plus | Strand::Unknown => a.cmp(&b),
        Strand::Minus => b.cmp(&a),
    }
}

fn cmp_local_read_in_transcription_direction(
    indices_sorted: &[usize],
    reads: &[ReadRecord],
    a: usize,
    b: usize,
) -> Ordering {
    let a_idx = indices_sorted[a];
    let b_idx = indices_sorted[b];
    let a_read = &reads[a_idx];
    let b_read = &reads[b_idx];

    cmp_in_transcription_direction(
        a_read.strand,
        five_prime_coord(a_read),
        five_prime_coord(b_read),
    )
    .then_with(|| {
        cmp_in_transcription_direction(
            a_read.strand,
            three_prime_coord(a_read),
            three_prime_coord(b_read),
        )
    })
    .then_with(|| a_read.id.cmp(&b_read.id))
    .then_with(|| a_idx.cmp(&b_idx))
}

fn weighted_median_endpoint(mut coordinates: Vec<u32>, strand: Strand) -> u32 {
    coordinates.sort_unstable_by(|a, b| cmp_in_transcription_direction(strand, *a, *b));
    coordinates[(coordinates.len() - 1) / 2]
}

fn modal_consensus_endpoint(coordinates: &[u32], strand: Strand, five_prime: bool) -> u32 {
    let median = weighted_median_endpoint(coordinates.to_vec(), strand);
    let mut counts: BTreeMap<u32, usize> = BTreeMap::new();
    for &coordinate in coordinates {
        *counts.entry(coordinate).or_default() += 1;
    }

    counts
        .into_iter()
        .max_by(|(coordinate_a, count_a), (coordinate_b, count_b)| {
            count_a
                .cmp(count_b)
                // Among equally supported coordinates, stay closest to the median
                // so ordinary jitter cannot pull the consensus to an extreme.
                .then_with(|| {
                    median
                        .abs_diff(*coordinate_b)
                        .cmp(&median.abs_diff(*coordinate_a))
                })
                // A remaining exact tie uses the transcript's outer boundary.
                .then_with(|| match (strand, five_prime) {
                    (Strand::Plus | Strand::Unknown, true) | (Strand::Minus, false) => {
                        coordinate_b.cmp(coordinate_a)
                    }
                    (Strand::Plus | Strand::Unknown, false) | (Strand::Minus, true) => {
                        coordinate_a.cmp(coordinate_b)
                    }
                })
        })
        .map(|(coordinate, _)| coordinate)
        .expect("members are non-empty")
}

fn consensus_for_members(
    members: &[usize],
    indices_sorted: &[usize],
    reads: &[ReadRecord],
) -> (Interval, TuEndpointStats) {
    debug_assert!(!members.is_empty());
    let strand = reads[indices_sorted[members[0]]].strand;
    let mut five_prime: Vec<u32> = Vec::with_capacity(members.len());
    let mut three_prime: Vec<u32> = Vec::with_capacity(members.len());

    for &local_pos in members {
        let read = &reads[indices_sorted[local_pos]];
        debug_assert_eq!(read.strand, strand);
        five_prime.push(five_prime_coord(read));
        three_prime.push(three_prime_coord(read));
    }

    let five_prime_min = *five_prime.iter().min().expect("members are non-empty");
    let five_prime_max = *five_prime.iter().max().expect("members are non-empty");
    let three_prime_min = *three_prime.iter().min().expect("members are non-empty");
    let three_prime_max = *three_prime.iter().max().expect("members are non-empty");
    let five_prime_consensus = modal_consensus_endpoint(&five_prime, strand, true);
    let three_prime_consensus = modal_consensus_endpoint(&three_prime, strand, false);

    let (start, end) = match strand {
        Strand::Plus | Strand::Unknown => (five_prime_consensus, three_prime_consensus),
        Strand::Minus => (three_prime_consensus, five_prime_consensus),
    };
    let interval = Interval::new(Coord::new(start), Coord::new(end))
        .expect("coordinate-wise endpoint modes of valid intervals form a valid interval");

    (
        interval,
        TuEndpointStats {
            support: members.len(),
            five_prime_consensus,
            three_prime_consensus,
            five_prime_min,
            five_prime_max,
            three_prime_min,
            three_prime_max,
        },
    )
}

/// Boundary consensus and endpoint stats for one final family.
///
/// The emitted interval is the modal endpoint pair of the family's core
/// members (anchor-quality reads, judged by span-Jaccard against the family
/// anchor at `members[0]`), matching how [`ConsensusFamily`] grows the
/// consensus incrementally. Support and endpoint min/max still describe every
/// member, so absorbed fragments remain visible in the stats.
fn family_consensus(
    members: &[usize],
    indices_sorted: &[usize],
    reads: &[ReadRecord],
    score1_threshold: f64,
) -> (Interval, TuEndpointStats) {
    debug_assert!(!members.is_empty());
    let anchor = reads[indices_sorted[members[0]]].interval;
    let core: Vec<usize> = members
        .iter()
        .copied()
        .filter(|&local_pos| {
            score1_interval(reads[indices_sorted[local_pos]].interval, anchor) >= score1_threshold
        })
        .collect();
    debug_assert!(!core.is_empty(), "the anchor itself is always core");
    let (interval, mut stats) = consensus_for_members(&core, indices_sorted, reads);

    stats.support = members.len();
    for &local_pos in members {
        let read = &reads[indices_sorted[local_pos]];
        let five_prime = five_prime_coord(read);
        let three_prime = three_prime_coord(read);
        stats.five_prime_min = stats.five_prime_min.min(five_prime);
        stats.five_prime_max = stats.five_prime_max.max(five_prime);
        stats.three_prime_min = stats.three_prime_min.min(three_prime);
        stats.three_prime_max = stats.three_prime_max.max(three_prime);
    }

    (interval, stats)
}

/// Strand-aware 3' agreement window for one read against a consensus.
///
/// Span-Jaccard grants two boundary-matched reads a combined endpoint budget of
/// `(1 - score1_threshold)` of their span, so full-length molecules absorb 3'
/// jitter proportional to their length. Truncated fragments qualify through the
/// contained-read paths instead and would otherwise be held to the fixed
/// attachment tolerance; scaling the window by read length keeps one TU's short
/// and long molecules under the same relative noise budget.
fn three_prime_window_bp(
    read_len: u32,
    score1_threshold: f64,
    options: TuClusteringOptions,
) -> u32 {
    let scaled = ((1.0 - score1_threshold).max(0.0) * read_len as f64).floor() as u32;
    options.three_prime_tolerance_bp.max(scaled)
}

fn interval_satisfies_consensus(
    read: &ReadRecord,
    consensus: Interval,
    score1_threshold: f64,
    score2_threshold: f64,
    options: TuClusteringOptions,
) -> bool {
    if score1_interval(read.interval, consensus) >= score1_threshold {
        return true;
    }
    if !options.attach_contained_reads {
        return false;
    }

    let window = three_prime_window_bp(read.interval.len(), score1_threshold, options);
    let read_five_prime = five_prime_coord(read);
    let read_three_prime = three_prime_coord(read);
    let consensus_five_prime = endpoint_coord(consensus, read.strand, true);
    let consensus_three_prime = endpoint_coord(consensus, read.strand, false);

    // A fragment contained in the consensus span (up to boundary jitter) is a
    // truncated molecule of this TU, not evidence for a new one. Keeping it
    // restores the 0.1.x absorption semantics that read-level span-Jaccard
    // chains and the overlap-over-longer ladders pooled, without their
    // boundary drift: an absorbed member never votes on the consensus, cannot
    // extend it, and qualification is always against the final consensus,
    // never against an intermediate read.
    let five_prime_overhang =
        match cmp_in_transcription_direction(read.strand, read_five_prime, consensus_five_prime) {
            Ordering::Less => read_five_prime.abs_diff(consensus_five_prime),
            Ordering::Greater | Ordering::Equal => 0,
        };
    let three_prime_overhang = match cmp_in_transcription_direction(
        read.strand,
        read_three_prime,
        consensus_three_prime,
    ) {
        Ordering::Greater => read_three_prime.abs_diff(consensus_three_prime),
        Ordering::Less | Ordering::Equal => 0,
    };
    if five_prime_overhang <= window && three_prime_overhang <= window {
        return true;
    }

    // 3'-anchored paths for reads reaching upstream of the consensus 5'.
    if read_three_prime.abs_diff(consensus_three_prime) > window {
        return false;
    }

    if score2_interval(read.interval, consensus) >= score2_threshold {
        return true;
    }

    options.max_five_prime_delta_bp.is_some_and(|max_delta| {
        five_prime_coord(read).abs_diff(endpoint_coord(consensus, read.strand, true)) <= max_delta
    })
}

/// An interval tree over the TUs in one reference/strand partition.
///
/// Assignment used to compare every read with every emitted TU. At direct-RNA
/// scale that made the post-clustering pass quadratic even though almost all
/// pairs are on different loci. Nodes keep the intervals crossing a median
/// start coordinate in two orders, allowing a query to visit only spatially
/// overlapping TUs (plus the logarithmic search path).
#[derive(Debug)]
struct TuIntervalIndexNode {
    center: u32,
    crossing_by_start: Vec<usize>,
    crossing_by_end: Vec<usize>,
    left: Option<Box<Self>>,
    right: Option<Box<Self>>,
}

impl TuIntervalIndexNode {
    /// Build from indices already ordered by interval start/end and original index.
    fn build_sorted(tu_indices: Vec<usize>, tus: &[Tu]) -> Option<Box<Self>> {
        if tu_indices.is_empty() {
            return None;
        }

        // A median start keeps both recursive sides balanced. Empty intervals
        // are omitted by TuIntervalIndex::build, so the interval supplying the
        // median always crosses its own start and recursion must make progress.
        let center = tus[tu_indices[tu_indices.len() / 2]].interval.start().get();

        let mut left = Vec::new();
        let mut right = Vec::new();
        let mut crossing_by_start = Vec::new();
        for tu_index in tu_indices {
            let interval = tus[tu_index].interval;
            if interval.end().get() <= center {
                left.push(tu_index);
            } else if interval.start().get() > center {
                right.push(tu_index);
            } else {
                crossing_by_start.push(tu_index);
            }
        }

        // Partitioning preserves the incoming start order.
        let mut crossing_by_end = crossing_by_start.clone();
        crossing_by_end.sort_unstable_by(|&a, &b| {
            tus[b]
                .interval
                .end()
                .cmp(&tus[a].interval.end())
                .then_with(|| tus[a].interval.start().cmp(&tus[b].interval.start()))
                .then_with(|| a.cmp(&b))
        });

        Some(Box::new(Self {
            center,
            crossing_by_start,
            crossing_by_end,
            left: Self::build_sorted(left, tus),
            right: Self::build_sorted(right, tus),
        }))
    }

    fn for_each_overlap<F>(&self, query: Interval, tus: &[Tu], visit: &mut F)
    where
        F: FnMut(usize),
    {
        let query_start = query.start().get();
        let query_end = query.end().get();

        if query_end <= self.center {
            // Every centered interval ends beyond `center`; start is the only
            // remaining overlap condition for a query entirely to the left.
            for &tu_index in &self.crossing_by_start {
                if tus[tu_index].interval.start().get() >= query_end {
                    break;
                }
                visit(tu_index);
            }
            if let Some(left) = self.left.as_deref() {
                left.for_each_overlap(query, tus, visit);
            }
        } else if query_start > self.center {
            // Symmetrically, centered starts are left of the query and only
            // their end coordinates can exclude them.
            for &tu_index in &self.crossing_by_end {
                if tus[tu_index].interval.end().get() <= query_start {
                    break;
                }
                visit(tu_index);
            }
            if let Some(right) = self.right.as_deref() {
                right.for_each_overlap(query, tus, visit);
            }
        } else {
            // The query contains the center, so every centered interval
            // overlaps it. Either recursive side may contain more overlaps.
            for &tu_index in &self.crossing_by_start {
                visit(tu_index);
            }
            if let Some(left) = self.left.as_deref() {
                left.for_each_overlap(query, tus, visit);
            }
            if let Some(right) = self.right.as_deref() {
                right.for_each_overlap(query, tus, visit);
            }
        }
    }
}

#[derive(Debug, Default)]
struct TuIntervalIndex {
    root: Option<Box<TuIntervalIndexNode>>,
}

impl TuIntervalIndex {
    fn build(tu_indices: Vec<usize>, tus: &[Tu]) -> Self {
        // Empty TUs can never have positive overlap with a read. Removing them
        // also guarantees progress when constructing the recursive tree.
        let mut non_empty: Vec<usize> = tu_indices
            .into_iter()
            .filter(|&tu_index| !tus[tu_index].interval.is_empty())
            .collect();
        non_empty.sort_unstable_by(|&a, &b| {
            tus[a]
                .interval
                .start()
                .cmp(&tus[b].interval.start())
                .then_with(|| tus[a].interval.end().cmp(&tus[b].interval.end()))
                .then_with(|| a.cmp(&b))
        });
        Self {
            root: TuIntervalIndexNode::build_sorted(non_empty, tus),
        }
    }

    fn for_each_overlap<F>(&self, query: Interval, tus: &[Tu], mut visit: F)
    where
        F: FnMut(usize),
    {
        if query.is_empty() {
            return;
        }
        if let Some(root) = self.root.as_deref() {
            root.for_each_overlap(query, tus, &mut visit);
        }
    }
}

fn build_tu_interval_indices(tus: &[Tu]) -> HashMap<String, HashMap<Strand, TuIntervalIndex>> {
    let mut partitions: HashMap<String, HashMap<Strand, Vec<usize>>> = HashMap::new();
    for (tu_index, tu) in tus.iter().enumerate() {
        if let Some(strands) = partitions.get_mut(tu.contig.as_str()) {
            strands.entry(tu.strand).or_default().push(tu_index);
        } else {
            let mut strands = HashMap::new();
            strands.insert(tu.strand, vec![tu_index]);
            partitions.insert(tu.contig.clone(), strands);
        }
    }

    partitions
        .into_iter()
        .map(|(contig, strands)| {
            let indices = strands
                .into_iter()
                .map(|(strand, tu_indices)| (strand, TuIntervalIndex::build(tu_indices, tus)))
                .collect();
            (contig, indices)
        })
        .collect()
}

fn assignment_candidate(
    read: &ReadRecord,
    tu: &Tu,
    tu_index: usize,
    score1_threshold: f64,
    score2_threshold: f64,
    options: TuClusteringOptions,
) -> Option<AssignmentCandidate> {
    if read.contig != tu.contig
        || read.strand != tu.strand
        || read.interval.overlap_len(tu.interval) == 0
    {
        return None;
    }

    let score1 = score1_interval(read.interval, tu.interval);
    let score2 = score2_interval(read.interval, tu.interval);
    let five_prime_delta_bp =
        five_prime_coord(read).abs_diff(endpoint_coord(tu.interval, read.strand, true));
    let three_prime_delta_bp =
        three_prime_coord(read).abs_diff(endpoint_coord(tu.interval, read.strand, false));

    let direct = score1 >= score1_threshold;
    // Mirror of `interval_satisfies_consensus`: contained fragments of a TU and
    // 3'-anchored near-matches count toward it as partial evidence.
    let window = three_prime_window_bp(read.interval.len(), score1_threshold, options);
    let five_prime_overhang = match cmp_in_transcription_direction(
        read.strand,
        five_prime_coord(read),
        endpoint_coord(tu.interval, read.strand, true),
    ) {
        Ordering::Less => five_prime_delta_bp,
        Ordering::Greater | Ordering::Equal => 0,
    };
    let three_prime_overhang = match cmp_in_transcription_direction(
        read.strand,
        three_prime_coord(read),
        endpoint_coord(tu.interval, read.strand, false),
    ) {
        Ordering::Greater => three_prime_delta_bp,
        Ordering::Less | Ordering::Equal => 0,
    };
    let contained = five_prime_overhang <= window && three_prime_overhang <= window;
    let anchored = three_prime_delta_bp <= window
        && (score2 >= score2_threshold
            || options
                .max_five_prime_delta_bp
                .is_some_and(|maximum| five_prime_delta_bp <= maximum));
    let partial = !direct && options.attach_contained_reads && (contained || anchored);
    if !direct && !partial {
        return None;
    }

    Some(AssignmentCandidate {
        tu_index,
        score1,
        score2,
        five_prime_delta_bp,
        three_prime_delta_bp,
        assignment_score: score2,
        partial,
    })
}

fn compare_assignment_candidates(
    left: &AssignmentCandidate,
    right: &AssignmentCandidate,
    tus: &[Tu],
) -> Ordering {
    // Biological evidence is compared before IDs. IDs make the output ordering total, but the
    // ambiguity decision below uses only the score margin and can therefore never be resolved by
    // this lexical tie-breaker.
    right
        .assignment_score
        .total_cmp(&left.assignment_score)
        .then_with(|| right.score1.total_cmp(&left.score1))
        .then_with(|| left.three_prime_delta_bp.cmp(&right.three_prime_delta_bp))
        .then_with(|| left.five_prime_delta_bp.cmp(&right.five_prime_delta_bp))
        .then_with(|| tus[left.tu_index].id.cmp(&tus[right.tu_index].id))
        .then_with(|| left.tu_index.cmp(&right.tu_index))
}

fn exact_equal_weights(candidates: &[AssignmentCandidate]) -> Vec<FractionalAssignment> {
    debug_assert!(!candidates.is_empty());
    let equal_weight = 1.0 / candidates.len() as f64;
    let mut accumulated = 0.0;
    let mut assignments = Vec::with_capacity(candidates.len());
    for (position, candidate) in candidates.iter().enumerate() {
        let weight = if position + 1 == candidates.len() {
            // Computing the final residual guarantees that summing the emitted weights in this
            // deterministic order produces exactly 1.0, including for non-power-of-two counts.
            1.0 - accumulated
        } else {
            equal_weight
        };
        accumulated += weight;
        assignments.push(FractionalAssignment {
            tu_index: candidate.tu_index,
            weight,
        });
    }
    debug_assert_eq!(
        assignments
            .iter()
            .map(|assignment| assignment.weight)
            .sum::<f64>(),
        1.0
    );
    assignments
}

fn classify_assignment_candidates(
    mut candidates: Vec<AssignmentCandidate>,
    tus: &[Tu],
    ambiguity_margin: f64,
    fractional_assignment: bool,
) -> ReadAssignment {
    candidates.sort_by(|left, right| compare_assignment_candidates(left, right, tus));

    let Some(best) = candidates.first().copied() else {
        return ReadAssignment {
            status: AssignmentStatus::Unassigned,
            best: None,
            second_best: None,
            score_margin: None,
            fractional_assignments: Vec::new(),
        };
    };
    let second_best = candidates.get(1).copied();
    let score_margin =
        second_best.map(|runner_up| (best.assignment_score - runner_up.assignment_score).max(0.0));
    let ambiguous = score_margin.is_some_and(|margin| margin <= ambiguity_margin);
    let status = if ambiguous {
        AssignmentStatus::Ambiguous
    } else if best.partial {
        AssignmentStatus::Partial
    } else {
        AssignmentStatus::Unique
    };

    let fractional_assignments = if ambiguous {
        if fractional_assignment {
            let candidates_in_margin: Vec<AssignmentCandidate> = candidates
                .iter()
                .copied()
                .take_while(|candidate| {
                    best.assignment_score - candidate.assignment_score <= ambiguity_margin
                })
                .collect();
            exact_equal_weights(&candidates_in_margin)
        } else {
            Vec::new()
        }
    } else {
        exact_equal_weights(&[best])
    };

    ReadAssignment {
        status,
        best: Some(best),
        second_best,
        score_margin,
        fractional_assignments,
    }
}

/// Classify every read against all emitted TU consensuses.
///
/// A candidate qualifies directly when score1 reaches `score1_threshold`, or as a partial
/// candidate when it satisfies the same contained-read rule used by clustering. Candidates are
/// ordered by score2, score1, endpoint deltas, and finally TU ID/index. The lexical fields only
/// stabilize reporting: a runner-up whose score2 is within `ambiguity_margin` always makes the
/// read ambiguous. When fractional assignment is enabled, every candidate within that margin
/// receives an equal share and the shares sum exactly to `1.0`.
pub fn assign_reads_to_tus(
    reads: &[ReadRecord],
    tus: &[Tu],
    score1_threshold: f64,
    score2_threshold: f64,
    options: TuClusteringOptions,
    ambiguity_margin: f64,
    fractional_assignment: bool,
) -> Result<Vec<ReadAssignment>, TuError> {
    if !ambiguity_margin.is_finite() || !(0.0..=1.0).contains(&ambiguity_margin) {
        return Err(TuError::InvalidAmbiguityMargin {
            value: ambiguity_margin,
        });
    }
    if !score1_threshold.is_finite() || !(0.0..=1.0).contains(&score1_threshold) {
        return Err(TuError::InvalidScore1Threshold {
            value: score1_threshold,
        });
    }
    if options.attach_contained_reads
        && (!score2_threshold.is_finite() || !(0.0..=1.0).contains(&score2_threshold))
    {
        return Err(TuError::InvalidScore2Threshold {
            value: score2_threshold,
        });
    }

    let tu_indices = build_tu_interval_indices(tus);

    reads
        .par_iter()
        .map(|read| {
            if read.strand == Strand::Unknown {
                return Err(TuError::UnknownStrand {
                    read_id: read.id.clone(),
                });
            }

            let mut candidates: Vec<AssignmentCandidate> = Vec::new();
            if let Some(index) = tu_indices
                .get(read.contig.as_str())
                .and_then(|strands| strands.get(&read.strand))
            {
                index.for_each_overlap(read.interval, tus, |tu_index| {
                    if let Some(candidate) = assignment_candidate(
                        read,
                        &tus[tu_index],
                        tu_index,
                        score1_threshold,
                        score2_threshold,
                        options,
                    ) {
                        candidates.push(candidate);
                    }
                });
            }

            Ok(classify_assignment_candidates(
                candidates,
                tus,
                ambiguity_margin,
                fractional_assignment,
            ))
        })
        .collect()
}

#[derive(Clone, Debug, Default)]
struct CoordinateMedian {
    lower: BTreeMap<u32, usize>,
    upper: BTreeMap<u32, usize>,
    lower_len: usize,
    upper_len: usize,
    reverse: bool,
}

impl CoordinateMedian {
    fn new(reverse: bool) -> Self {
        Self {
            reverse,
            ..Self::default()
        }
    }

    fn increment(map: &mut BTreeMap<u32, usize>, coordinate: u32) {
        *map.entry(coordinate).or_default() += 1;
    }

    fn decrement(map: &mut BTreeMap<u32, usize>, coordinate: u32) {
        let remove = {
            let count = map
                .get_mut(&coordinate)
                .expect("coordinate must exist in median multiset");
            *count -= 1;
            *count == 0
        };
        if remove {
            map.remove(&coordinate);
        }
    }

    fn desired_lower_len(&self) -> usize {
        let total = self.lower_len + self.upper_len;
        if total == 0 {
            0
        } else if self.reverse {
            // The lower median in descending transcription order is the
            // upper genomic median when the observation count is even.
            total / 2 + 1
        } else {
            total.div_ceil(2)
        }
    }

    fn rebalance(&mut self) {
        let desired = self.desired_lower_len();
        while self.lower_len > desired {
            let coordinate = *self
                .lower
                .keys()
                .next_back()
                .expect("an oversized lower partition cannot be empty");
            Self::decrement(&mut self.lower, coordinate);
            self.lower_len -= 1;
            Self::increment(&mut self.upper, coordinate);
            self.upper_len += 1;
        }
        while self.lower_len < desired {
            let coordinate = *self
                .upper
                .keys()
                .next()
                .expect("an undersized lower partition requires an upper value");
            Self::decrement(&mut self.upper, coordinate);
            self.upper_len -= 1;
            Self::increment(&mut self.lower, coordinate);
            self.lower_len += 1;
        }
    }

    fn insert(&mut self, coordinate: u32) {
        let belongs_in_lower = self
            .lower
            .keys()
            .next_back()
            .is_none_or(|&lower_max| coordinate <= lower_max);
        if belongs_in_lower {
            Self::increment(&mut self.lower, coordinate);
            self.lower_len += 1;
        } else {
            Self::increment(&mut self.upper, coordinate);
            self.upper_len += 1;
        }
        self.rebalance();
    }

    fn remove(&mut self, coordinate: u32) {
        if self.lower.contains_key(&coordinate) {
            Self::decrement(&mut self.lower, coordinate);
            self.lower_len -= 1;
        } else {
            Self::decrement(&mut self.upper, coordinate);
            self.upper_len -= 1;
        }
        self.rebalance();
    }

    fn median(&self) -> u32 {
        *self
            .lower
            .keys()
            .next_back()
            .expect("endpoint consensus requires at least one observation")
    }
}

#[derive(Clone, Debug)]
struct EndpointConsensusState {
    median: CoordinateMedian,
    counts: BTreeMap<u32, usize>,
    coordinates_by_count: BTreeMap<usize, BTreeSet<u32>>,
    prefer_lower: bool,
}

impl EndpointConsensusState {
    fn new(strand: Strand, prefer_lower: bool) -> Self {
        Self {
            median: CoordinateMedian::new(strand == Strand::Minus),
            counts: BTreeMap::new(),
            coordinates_by_count: BTreeMap::new(),
            prefer_lower,
        }
    }

    fn remove_count_group(&mut self, count: usize, coordinate: u32) {
        let remove_group = {
            let coordinates = self
                .coordinates_by_count
                .get_mut(&count)
                .expect("endpoint frequency group must exist");
            assert!(coordinates.remove(&coordinate));
            coordinates.is_empty()
        };
        if remove_group {
            self.coordinates_by_count.remove(&count);
        }
    }

    fn insert(&mut self, coordinate: u32) {
        self.median.insert(coordinate);
        let old_count = self.counts.get(&coordinate).copied().unwrap_or(0);
        if old_count > 0 {
            self.remove_count_group(old_count, coordinate);
        }
        let new_count = old_count + 1;
        self.counts.insert(coordinate, new_count);
        self.coordinates_by_count
            .entry(new_count)
            .or_default()
            .insert(coordinate);
    }

    fn remove(&mut self, coordinate: u32) {
        self.median.remove(coordinate);
        let old_count = self.counts[&coordinate];
        self.remove_count_group(old_count, coordinate);
        if old_count == 1 {
            self.counts.remove(&coordinate);
        } else {
            let new_count = old_count - 1;
            self.counts.insert(coordinate, new_count);
            self.coordinates_by_count
                .entry(new_count)
                .or_default()
                .insert(coordinate);
        }
    }

    fn mode(&self) -> u32 {
        let median = self.median.median();
        let modal_coordinates = self
            .coordinates_by_count
            .last_key_value()
            .map(|(_, coordinates)| coordinates)
            .expect("endpoint consensus requires at least one frequency group");
        let below = modal_coordinates.range(..=median).next_back().copied();
        let above = modal_coordinates.range(median..).next().copied();

        match (below, above) {
            (Some(left), Some(right)) => {
                let left_distance = median.abs_diff(left);
                let right_distance = median.abs_diff(right);
                match left_distance.cmp(&right_distance) {
                    Ordering::Less => left,
                    Ordering::Greater => right,
                    Ordering::Equal if self.prefer_lower => left.min(right),
                    Ordering::Equal => left.max(right),
                }
            }
            (Some(coordinate), None) | (None, Some(coordinate)) => coordinate,
            (None, None) => unreachable!("a modal frequency group cannot be empty"),
        }
    }
}

/// One growing family during the final consensus split.
///
/// Boundary consensus is voted on exclusively by "core" members: reads whose
/// span-Jaccard against the fixed family anchor clears `score1_threshold`.
/// Absorbed fragments (contained reads and 3'-anchored near-matches) join the
/// family and count toward support, but never shift the boundary. This keeps
/// the emitted TU anchored on its full-length molecules: a dense hotspot of
/// truncation fragments cannot outvote the boundary mode and evict the
/// full-length cloud, which mirrors how 0.1.x derived the TU span from the
/// representative read while absorbing fragments silently.
#[derive(Clone, Debug)]
struct ConsensusFamily {
    members: Vec<usize>,
    five_prime: EndpointConsensusState,
    three_prime: EndpointConsensusState,
    consensus: Interval,
}

#[derive(Clone, Copy, Debug)]
struct ConsensusSplitConfig {
    score1_threshold: f64,
    score2_threshold: f64,
    options: TuClusteringOptions,
}

impl ConsensusFamily {
    fn new(local_pos: usize, indices_sorted: &[usize], reads: &[ReadRecord]) -> Self {
        let read = &reads[indices_sorted[local_pos]];
        let five_prime_prefers_lower = read.strand != Strand::Minus;
        let three_prime_prefers_lower = read.strand == Strand::Minus;
        let mut five_prime = EndpointConsensusState::new(read.strand, five_prime_prefers_lower);
        let mut three_prime = EndpointConsensusState::new(read.strand, three_prime_prefers_lower);
        five_prime.insert(five_prime_coord(read));
        three_prime.insert(three_prime_coord(read));
        Self {
            members: vec![local_pos],
            five_prime,
            three_prime,
            consensus: read.interval,
        }
    }

    fn consensus_interval(&self, strand: Strand) -> Interval {
        let five_prime = self.five_prime.mode();
        let three_prime = self.three_prime.mode();
        let (start, end) = match strand {
            Strand::Plus | Strand::Unknown => (five_prime, three_prime),
            Strand::Minus => (three_prime, five_prime),
        };
        Interval::new(Coord::new(start), Coord::new(end))
            .expect("coordinate-wise endpoint modes of valid intervals form a valid interval")
    }

    fn try_push(
        &mut self,
        local_pos: usize,
        indices_sorted: &[usize],
        reads: &[ReadRecord],
        config: ConsensusSplitConfig,
    ) -> bool {
        let candidate = &reads[indices_sorted[local_pos]];
        let anchor = &reads[indices_sorted[self.members[0]]];
        if !interval_satisfies_consensus(
            candidate,
            anchor.interval,
            config.score1_threshold,
            config.score2_threshold,
            config.options,
        ) {
            return false;
        }

        // Only anchor-quality reads vote on the boundary; absorbed fragments
        // must qualify against the standing consensus but cannot move it.
        let is_core =
            score1_interval(candidate.interval, anchor.interval) >= config.score1_threshold;
        if !is_core {
            if interval_satisfies_consensus(
                candidate,
                self.consensus,
                config.score1_threshold,
                config.score2_threshold,
                config.options,
            ) {
                self.members.push(local_pos);
                return true;
            }
            return false;
        }

        let five_prime = five_prime_coord(candidate);
        let three_prime = three_prime_coord(candidate);
        self.five_prime.insert(five_prime);
        self.three_prime.insert(three_prime);
        let prospective_consensus = self.consensus_interval(candidate.strand);

        // Existing members were already checked against the current consensus.
        // A stable mode therefore needs only one new check; the former code
        // cloned, re-sorted, and rescanned the whole growing family here.
        let all_direct = if prospective_consensus == self.consensus {
            interval_satisfies_consensus(
                candidate,
                prospective_consensus,
                config.score1_threshold,
                config.score2_threshold,
                config.options,
            )
        } else {
            self.members
                .iter()
                .copied()
                .chain(std::iter::once(local_pos))
                .all(|member_pos| {
                    interval_satisfies_consensus(
                        &reads[indices_sorted[member_pos]],
                        prospective_consensus,
                        config.score1_threshold,
                        config.score2_threshold,
                        config.options,
                    )
                })
        };

        if all_direct {
            self.members.push(local_pos);
            self.consensus = prospective_consensus;
            true
        } else {
            self.five_prime.remove(five_prime);
            self.three_prime.remove(three_prime);
            false
        }
    }
}

/// Order component members so the best-supported exact boundary pair is processed first.
///
/// Members are grouped by identical read intervals (equal strand-aware endpoint pairs inside
/// one partition) and groups are emitted by descending support, breaking ties in transcription
/// direction. Family anchors therefore sit on modal boundary evidence instead of on whichever
/// read happens to start a component in genomic order.
fn members_in_support_order(
    members: &[usize],
    indices_sorted: &[usize],
    reads: &[ReadRecord],
) -> Vec<usize> {
    let mut ordered = members.to_vec();
    ordered
        .sort_by(|&a, &b| cmp_local_read_in_transcription_direction(indices_sorted, reads, a, b));

    let mut groups: Vec<Vec<usize>> = Vec::new();
    for local_pos in ordered {
        let interval = reads[indices_sorted[local_pos]].interval;
        match groups.last_mut() {
            Some(group) if reads[indices_sorted[group[0]]].interval == interval => {
                group.push(local_pos);
            }
            _ => groups.push(vec![local_pos]),
        }
    }

    let mut group_order: Vec<usize> = (0..groups.len()).collect();
    group_order.sort_by(|&a, &b| groups[b].len().cmp(&groups[a].len()).then(a.cmp(&b)));

    group_order
        .into_iter()
        .flat_map(|group_idx| std::mem::take(&mut groups[group_idx]))
        .collect()
}

/// Split a connected component into deterministic anchor-bounded consensus families.
///
/// Members are processed in support order, so the strongest exact boundary pair founds the
/// first family and anchors it. A family can grow only when the candidate qualifies directly
/// against its fixed anchor and every prospective member qualifies directly against the
/// prospective support-mode consensus. This is the final-cluster invariant: no membership can
/// be justified solely by a path through intermediate reads.
fn split_members_by_consensus(
    members: &[usize],
    indices_sorted: &[usize],
    reads: &[ReadRecord],
    score1_threshold: f64,
    score2_threshold: f64,
    options: TuClusteringOptions,
) -> Vec<Vec<usize>> {
    let config = ConsensusSplitConfig {
        score1_threshold,
        score2_threshold,
        options,
    };
    let mut families: Vec<ConsensusFamily> = Vec::new();
    for local_pos in members_in_support_order(members, indices_sorted, reads) {
        let mut selected_family = None;

        for (family_idx, family) in families.iter_mut().enumerate() {
            if family.try_push(local_pos, indices_sorted, reads, config) {
                selected_family = Some(family_idx);
                break;
            }
        }

        if selected_family.is_none() {
            families.push(ConsensusFamily::new(local_pos, indices_sorted, reads));
        }
    }

    families.into_iter().map(|family| family.members).collect()
}

fn representative_read_index(
    members: &[usize],
    indices_sorted: &[usize],
    reads: &[ReadRecord],
) -> usize {
    let mut ordered = members.to_vec();
    ordered
        .sort_by(|&a, &b| cmp_local_read_in_transcription_direction(indices_sorted, reads, a, b));
    let mut representative = indices_sorted[ordered[0]];
    for &local_pos in &ordered[1..] {
        let candidate = indices_sorted[local_pos];
        if better_representative(&reads[candidate], &reads[representative]) {
            representative = candidate;
        }
    }
    representative
}

/// One-sided 3' gate for cluster attachment, as in 0.1.x.
///
/// The child may terminate anywhere inside the candidate's span (a truncated
/// molecule), but may not overhang the candidate's 3' end by more than the
/// tolerance. Attachment chains through `resolve_root`, so this is what lets
/// degradation ladders pool with their full-length parent; the final consensus
/// split then decides membership against the family consensus, so pooling
/// breadth cannot leak into TU boundaries.
fn candidate_within_three_prime_tolerance(
    child: &ReadRecord,
    candidate: &ReadRecord,
    tolerance_bp: u32,
) -> bool {
    match cmp_in_transcription_direction(
        child.strand,
        three_prime_coord(child),
        three_prime_coord(candidate),
    ) {
        Ordering::Greater => {
            three_prime_coord(child).abs_diff(three_prime_coord(candidate)) <= tolerance_bp
        }
        Ordering::Less | Ordering::Equal => true,
    }
}

fn candidate_five_prime_delta(child: &ReadRecord, candidate: &ReadRecord) -> u32 {
    five_prime_coord(child).abs_diff(five_prime_coord(candidate))
}

#[derive(Clone, Copy, Debug)]
struct AttachmentCandidate {
    parent_cluster_idx: usize,
    score2: f64,
    score1: f64,
    three_prime_delta: u32,
    five_prime_delta: u32,
}

fn orient_attachment_pair(
    a_cluster_idx: usize,
    b_cluster_idx: usize,
    clusters: &[SeedCluster],
    reads: &[ReadRecord],
) -> (usize, usize) {
    let a_rep_idx = clusters[a_cluster_idx].rep_read_index;
    let b_rep_idx = clusters[b_cluster_idx].rep_read_index;
    let a = &reads[a_rep_idx];
    let b = &reads[b_rep_idx];

    match read_len(a).cmp(&read_len(b)) {
        Ordering::Less => (a_cluster_idx, b_cluster_idx),
        Ordering::Greater => (b_cluster_idx, a_cluster_idx),
        Ordering::Equal => {
            // Equal-length pairs attach toward the lexically smaller representative ID.
            // The read index makes the orientation total even for duplicate IDs.
            if (a.id.as_str(), a_rep_idx) <= (b.id.as_str(), b_rep_idx) {
                (b_cluster_idx, a_cluster_idx)
            } else {
                (a_cluster_idx, b_cluster_idx)
            }
        }
    }
}

fn compare_attachment_candidates(
    a: &AttachmentCandidate,
    b: &AttachmentCandidate,
    clusters: &[SeedCluster],
    reads: &[ReadRecord],
) -> Ordering {
    let a_rep_idx = clusters[a.parent_cluster_idx].rep_read_index;
    let b_rep_idx = clusters[b.parent_cluster_idx].rep_read_index;
    let a_id = reads[a_rep_idx].id.as_str();
    let b_id = reads[b_rep_idx].id.as_str();

    // Preferred candidates sort first. Scores are descending; endpoint distances,
    // representative IDs, and indices are ascending deterministic tie-breakers.
    b.score2
        .total_cmp(&a.score2)
        .then_with(|| b.score1.total_cmp(&a.score1))
        .then_with(|| a.three_prime_delta.cmp(&b.three_prime_delta))
        .then_with(|| a.five_prime_delta.cmp(&b.five_prime_delta))
        .then_with(|| a_id.cmp(b_id))
        .then_with(|| a_rep_idx.cmp(&b_rep_idx))
        .then_with(|| a.parent_cluster_idx.cmp(&b.parent_cluster_idx))
}

#[derive(Clone, Debug)]
struct TuNoId {
    contig: String,
    strand: Strand,
    interval: Interval,
    rep_read_index: usize,
    endpoint_stats: TuEndpointStats,
}

#[derive(Clone, Debug)]
struct PartitionClusteringResult {
    tus: Vec<TuNoId>,
    local_read_to_tu: Vec<usize>,
}

#[derive(Clone, Debug)]
struct FinalConsensusFamily {
    rep_read_index: usize,
    members: Vec<usize>,
    interval: Interval,
    endpoint_stats: TuEndpointStats,
}

fn final_consensus_family(
    members: Vec<usize>,
    indices_sorted: &[usize],
    reads: &[ReadRecord],
    score1_threshold: f64,
    score2_threshold: f64,
    options: TuClusteringOptions,
) -> FinalConsensusFamily {
    let rep_read_index = representative_read_index(&members, indices_sorted, reads);
    let (interval, endpoint_stats) =
        family_consensus(&members, indices_sorted, reads, score1_threshold);
    debug_assert!(members.iter().all(|&local_pos| {
        interval_satisfies_consensus(
            &reads[indices_sorted[local_pos]],
            interval,
            score1_threshold,
            score2_threshold,
            options,
        )
    }));
    FinalConsensusFamily {
        rep_read_index,
        members,
        interval,
        endpoint_stats,
    }
}

/// Coalesce families that emit the same consensus interval.
///
/// Separate anchor-bounded families can converge on identical endpoint modes. Since every
/// member already qualifies directly against that shared consensus, retaining both families
/// would emit duplicate biological TUs and duplicate coordinate-derived IDs.
fn coalesce_identical_consensus_families(
    families: Vec<FinalConsensusFamily>,
    indices_sorted: &[usize],
    reads: &[ReadRecord],
    score1_threshold: f64,
    score2_threshold: f64,
    options: TuClusteringOptions,
) -> Vec<FinalConsensusFamily> {
    let mut by_interval: BTreeMap<Interval, Vec<FinalConsensusFamily>> = BTreeMap::new();
    for family in families {
        by_interval.entry(family.interval).or_default().push(family);
    }

    let mut coalesced = Vec::new();
    for (shared_interval, mut matching) in by_interval {
        if matching.len() == 1 {
            coalesced.push(matching.pop().expect("one matching family exists"));
            continue;
        }

        let merged_members: Vec<usize> = matching
            .iter()
            .flat_map(|family| family.members.iter().copied())
            .collect();
        let (merged_interval, _) =
            family_consensus(&merged_members, indices_sorted, reads, score1_threshold);
        let merged_is_valid = merged_interval == shared_interval
            && merged_members.iter().all(|&local_pos| {
                interval_satisfies_consensus(
                    &reads[indices_sorted[local_pos]],
                    merged_interval,
                    score1_threshold,
                    score2_threshold,
                    options,
                )
            });

        if merged_is_valid {
            coalesced.push(final_consensus_family(
                merged_members,
                indices_sorted,
                reads,
                score1_threshold,
                score2_threshold,
                options,
            ));
        } else {
            // Preserve the validated families if their union changes the modal consensus.
            coalesced.extend(matching);
        }
    }

    coalesced.sort_by(|a, b| cmp_read_idx(reads, a.rep_read_index, b.rep_read_index));
    coalesced
}

fn partition_into_regions(indices_sorted: &[usize], reads: &[ReadRecord]) -> Vec<(usize, usize)> {
    if indices_sorted.is_empty() {
        return Vec::new();
    }

    let mut regions: Vec<(usize, usize)> = Vec::new();
    let mut start = 0usize;
    let mut current_end = reads[indices_sorted[0]].interval.end().get();

    for pos in 1..indices_sorted.len() {
        let read = &reads[indices_sorted[pos]];
        let read_start = read.interval.start().get();
        if read_start >= current_end {
            regions.push((start, pos));
            start = pos;
            current_end = read.interval.end().get();
        } else {
            current_end = current_end.max(read.interval.end().get());
        }
    }
    regions.push((start, indices_sorted.len()));

    regions
}

fn cluster_tus_region(
    indices_sorted: &[usize],
    reads: &[ReadRecord],
    score1_threshold: f64,
    score2_threshold: f64,
    options: TuClusteringOptions,
) -> PartitionClusteringResult {
    let mut dsu = Dsu::new(indices_sorted.len());
    cluster_partition_score1(indices_sorted, reads, &mut dsu, score1_threshold);

    let mut clusters_by_root: HashMap<usize, Vec<usize>> = HashMap::new();
    for local_pos in 0..indices_sorted.len() {
        let root = dsu.find(local_pos);
        clusters_by_root.entry(root).or_default().push(local_pos);
    }

    // Seed clusters are the raw span-Jaccard components, exactly as in 0.1.x, so
    // the attachment pass sees the same merge topology: dense truncation ladders
    // stay connected at read level and containment chains can reach the parent.
    // The consensus invariant is enforced once, after pooling, by the final
    // split below.
    let mut clusters: Vec<SeedCluster> = Vec::with_capacity(clusters_by_root.len());
    for members in clusters_by_root.into_values() {
        clusters.push(SeedCluster {
            rep_read_index: representative_read_index(&members, indices_sorted, reads),
            members,
        });
    }

    clusters.sort_by(|a, b| cmp_read_idx(reads, a.rep_read_index, b.rep_read_index));

    let mut parents: Vec<Option<usize>> = vec![None; clusters.len()];
    if options.attach_contained_reads {
        let mut active: HashSet<usize> = HashSet::new();
        let mut ends: BinaryHeap<Reverse<(u32, usize)>> = BinaryHeap::new();
        let mut best_parents: Vec<Option<AttachmentCandidate>> = vec![None; clusters.len()];

        for cluster_idx in 0..clusters.len() {
            let rep_idx = clusters[cluster_idx].rep_read_index;
            let current = &reads[rep_idx];
            if current.interval.is_empty() {
                continue;
            }

            expire_active(current.interval.start().get(), &mut active, &mut ends);

            for &other_cluster_idx in &active {
                let (child_cluster_idx, parent_cluster_idx) =
                    orient_attachment_pair(cluster_idx, other_cluster_idx, &clusters, reads);
                let child = &reads[clusters[child_cluster_idx].rep_read_index];
                let candidate = &reads[clusters[parent_cluster_idx].rep_read_index];

                if !candidate_within_three_prime_tolerance(
                    child,
                    candidate,
                    options.three_prime_tolerance_bp,
                ) {
                    continue;
                }

                let five_prime_delta = candidate_five_prime_delta(child, candidate);
                let five_prime_override = options
                    .max_five_prime_delta_bp
                    .map(|max_delta| five_prime_delta <= max_delta)
                    .unwrap_or(false);

                let s2 = score2_interval(child.interval, candidate.interval);
                if s2 < score2_threshold && !five_prime_override {
                    continue;
                }
                let s1 = score1_interval(child.interval, candidate.interval);

                let attachment = AttachmentCandidate {
                    parent_cluster_idx,
                    score2: s2,
                    score1: s1,
                    three_prime_delta: three_prime_coord(child)
                        .abs_diff(three_prime_coord(candidate)),
                    five_prime_delta,
                };
                let replace = best_parents[child_cluster_idx]
                    .as_ref()
                    .map(|best| {
                        compare_attachment_candidates(&attachment, best, &clusters, reads)
                            == Ordering::Less
                    })
                    .unwrap_or(true);
                if replace {
                    best_parents[child_cluster_idx] = Some(attachment);
                }
            }

            add_active(
                current.interval.end().get(),
                cluster_idx,
                &mut active,
                &mut ends,
            );
        }

        for (child_cluster_idx, best) in best_parents.into_iter().enumerate() {
            if let Some(best) = best {
                parents[child_cluster_idx] = Some(best.parent_cluster_idx);
            }
        }
    }

    let mut cluster_to_root: Vec<usize> = Vec::with_capacity(clusters.len());
    for idx in 0..clusters.len() {
        cluster_to_root.push(resolve_root(idx, &mut parents));
    }

    let mut members_by_root: HashMap<usize, Vec<usize>> = HashMap::new();
    for (cluster_idx, cluster) in clusters.iter().enumerate() {
        let root = cluster_to_root[cluster_idx];
        members_by_root
            .entry(root)
            .or_default()
            .extend(cluster.members.iter().copied());
    }

    let mut final_families: Vec<FinalConsensusFamily> = Vec::new();
    for members in members_by_root.into_values() {
        for family in split_members_by_consensus(
            &members,
            indices_sorted,
            reads,
            score1_threshold,
            score2_threshold,
            options,
        ) {
            final_families.push(final_consensus_family(
                family,
                indices_sorted,
                reads,
                score1_threshold,
                score2_threshold,
                options,
            ));
        }
    }
    final_families = coalesce_identical_consensus_families(
        final_families,
        indices_sorted,
        reads,
        score1_threshold,
        score2_threshold,
        options,
    );

    let mut tus: Vec<TuNoId> = Vec::with_capacity(final_families.len());
    let mut local_read_to_tu: Vec<usize> = vec![usize::MAX; indices_sorted.len()];
    for (tu_index, family) in final_families.into_iter().enumerate() {
        let rep = &reads[family.rep_read_index];
        tus.push(TuNoId {
            contig: rep.contig.clone(),
            strand: rep.strand,
            interval: family.interval,
            rep_read_index: family.rep_read_index,
            endpoint_stats: family.endpoint_stats,
        });
        for local_pos in family.members {
            local_read_to_tu[local_pos] = tu_index;
        }
    }

    PartitionClusteringResult {
        tus,
        local_read_to_tu,
    }
}

#[derive(Clone, Copy, Debug, Default)]
struct PartitionStats {
    region_count: usize,
    max_region_reads: usize,
}

fn cluster_tus_partition(
    indices_sorted: &[usize],
    reads: &[ReadRecord],
    score1_threshold: f64,
    score2_threshold: f64,
    options: TuClusteringOptions,
) -> (PartitionClusteringResult, PartitionStats) {
    let regions = partition_into_regions(indices_sorted, reads);
    let stats = PartitionStats {
        region_count: regions.len(),
        max_region_reads: regions
            .iter()
            .map(|&(start, end)| end - start)
            .max()
            .unwrap_or(0),
    };

    if regions.is_empty() {
        return (
            PartitionClusteringResult {
                tus: Vec::new(),
                local_read_to_tu: Vec::new(),
            },
            stats,
        );
    }

    let region_results: Vec<PartitionClusteringResult> = regions
        .par_iter()
        .map(|&(start, end)| {
            cluster_tus_region(
                &indices_sorted[start..end],
                reads,
                score1_threshold,
                score2_threshold,
                options,
            )
        })
        .collect();

    let total_tus: usize = region_results.iter().map(|r| r.tus.len()).sum();
    let mut tus: Vec<TuNoId> = Vec::with_capacity(total_tus);
    let mut local_read_to_tu: Vec<usize> = Vec::with_capacity(indices_sorted.len());

    let mut tu_offset: usize = 0;
    for region_result in region_results {
        let PartitionClusteringResult {
            tus: region_tus,
            local_read_to_tu: region_read_to_tu,
        } = region_result;

        let region_tu_count = region_tus.len();
        tus.extend(region_tus);

        local_read_to_tu.extend(
            region_read_to_tu
                .into_iter()
                .map(|local_tu_idx| tu_offset + local_tu_idx),
        );

        tu_offset += region_tu_count;
    }

    (
        PartitionClusteringResult {
            tus,
            local_read_to_tu,
        },
        stats,
    )
}

/// Cluster reads with explicit second-pass options and return scale statistics.
pub fn cluster_tus_with_stats_options(
    reads: &[ReadRecord],
    score1_threshold: f64,
    score2_threshold: f64,
    options: TuClusteringOptions,
) -> Result<(TuClusteringResult, TuClusteringStats), TuError> {
    if !score1_threshold.is_finite() || !(0.0..=1.0).contains(&score1_threshold) {
        return Err(TuError::InvalidScore1Threshold {
            value: score1_threshold,
        });
    }
    if options.attach_contained_reads
        && (!score2_threshold.is_finite() || !(0.0..=1.0).contains(&score2_threshold))
    {
        return Err(TuError::InvalidScore2Threshold {
            value: score2_threshold,
        });
    }
    if let Some(read) = reads.iter().find(|read| read.strand == Strand::Unknown) {
        return Err(TuError::UnknownStrand {
            read_id: read.id.clone(),
        });
    }

    if reads.is_empty() {
        return Ok((
            TuClusteringResult::new(Vec::new(), Vec::new(), Vec::new()),
            TuClusteringStats::default(),
        ));
    }

    let mut indices: Vec<usize> = (0..reads.len()).collect();
    indices.par_sort_unstable_by(|&i, &j| cmp_read_idx(reads, i, j));

    let mut partitions: Vec<(usize, usize)> = Vec::new();
    let mut start = 0;
    while start < indices.len() {
        let first_idx = indices[start];
        let contig = reads[first_idx].contig.as_str();
        let strand = reads[first_idx].strand;

        let mut end = start + 1;
        while end < indices.len() {
            let idx = indices[end];
            if reads[idx].strand != strand || reads[idx].contig.as_str() != contig {
                break;
            }
            end += 1;
        }

        partitions.push((start, end));
        start = end;
    }

    let partition_count = partitions.len();
    let max_partition_reads = partitions
        .iter()
        .map(|&(start, end)| end - start)
        .max()
        .unwrap_or(0);

    let partition_results: Vec<(PartitionClusteringResult, PartitionStats)> = partitions
        .par_iter()
        .map(|&(start, end)| {
            cluster_tus_partition(
                &indices[start..end],
                reads,
                score1_threshold,
                score2_threshold,
                options,
            )
        })
        .collect();

    let total_tus: usize = partition_results.iter().map(|(r, _)| r.tus.len()).sum();
    let width = 6usize.max(total_tus.to_string().len());

    let region_count: usize = partition_results.iter().map(|(_, s)| s.region_count).sum();
    let max_region_reads: usize = partition_results
        .iter()
        .map(|(_, s)| s.max_region_reads)
        .max()
        .unwrap_or(0);

    let mut tus: Vec<Tu> = Vec::with_capacity(total_tus);
    let mut endpoint_stats: Vec<TuEndpointStats> = Vec::with_capacity(total_tus);
    let mut read_to_tu: Vec<usize> = vec![usize::MAX; reads.len()];
    let mut tu_offset: usize = 0;

    for (partition_idx, (result, _)) in partition_results.into_iter().enumerate() {
        let PartitionClusteringResult {
            tus: tus_no_id,
            local_read_to_tu,
        } = result;

        let partition_tu_count = tus_no_id.len();
        for (local_tu_idx, tu) in tus_no_id.into_iter().enumerate() {
            endpoint_stats.push(tu.endpoint_stats);
            tus.push(Tu {
                id: format_tu_id(tu_offset + local_tu_idx, width),
                contig: tu.contig,
                strand: tu.strand,
                interval: tu.interval,
                rep_read_index: tu.rep_read_index,
            });
        }

        let (start, end) = partitions[partition_idx];
        let indices_sorted = &indices[start..end];
        for (local_pos, &read_idx) in indices_sorted.iter().enumerate() {
            read_to_tu[read_idx] = tu_offset + local_read_to_tu[local_pos];
        }

        tu_offset += partition_tu_count;
    }

    Ok((
        TuClusteringResult::new(tus, read_to_tu, endpoint_stats),
        TuClusteringStats {
            partition_count,
            max_partition_reads,
            region_count,
            max_region_reads,
        },
    ))
}

/// Cluster reads with default options and return scale statistics.
pub fn cluster_tus_with_stats(
    reads: &[ReadRecord],
    score1_threshold: f64,
    score2_threshold: f64,
) -> Result<(TuClusteringResult, TuClusteringStats), TuError> {
    cluster_tus_with_stats_options(
        reads,
        score1_threshold,
        score2_threshold,
        TuClusteringOptions::default(),
    )
}

/// Cluster reads with explicit second-pass options.
pub fn cluster_tus_with_options(
    reads: &[ReadRecord],
    score1_threshold: f64,
    score2_threshold: f64,
    options: TuClusteringOptions,
) -> Result<TuClusteringResult, TuError> {
    Ok(cluster_tus_with_stats_options(reads, score1_threshold, score2_threshold, options)?.0)
}

/// Cluster reads with default second-pass options.
///
/// # Examples
///
/// ```
/// use trackclustertu::model::{Coord, Interval, Strand};
/// use trackclustertu::tu::{cluster_tus, ReadRecord};
///
/// let reads = vec![ReadRecord {
///     contig: "chromosome".into(),
///     strand: Strand::Plus,
///     interval: Interval::new(Coord::new(100), Coord::new(500))?,
///     id: "read-1".into(),
/// }];
/// let result = cluster_tus(&reads, 0.95, 0.80)?;
/// assert_eq!(result.tus().len(), 1);
/// assert_eq!(result.read_to_tu(), &[0]);
/// assert_eq!(result.tu_for_read(0), result.tus().first());
/// # Ok::<(), Box<dyn std::error::Error>>(())
/// ```
pub fn cluster_tus(
    reads: &[ReadRecord],
    score1_threshold: f64,
    score2_threshold: f64,
) -> Result<TuClusteringResult, TuError> {
    cluster_tus_with_options(
        reads,
        score1_threshold,
        score2_threshold,
        TuClusteringOptions::default(),
    )
}

#[cfg(test)]
mod tests {
    use std::collections::{BTreeMap, HashMap};

    use crate::model::Coord;
    use proptest::prelude::*;

    use super::*;

    fn read(contig: &str, strand: Strand, start: u32, end: u32, id: &str) -> ReadRecord {
        ReadRecord {
            contig: contig.to_owned(),
            strand,
            interval: Interval::new(Coord::new(start), Coord::new(end)).unwrap(),
            id: id.to_owned(),
        }
    }

    fn tu(contig: &str, strand: Strand, start: u32, end: u32, id: &str) -> Tu {
        Tu {
            id: id.to_owned(),
            contig: contig.to_owned(),
            strand,
            interval: Interval::new(Coord::new(start), Coord::new(end)).unwrap(),
            rep_read_index: 0,
        }
    }

    fn assign_reads_to_tus_naive_for_test(
        reads: &[ReadRecord],
        tus: &[Tu],
        score1_threshold: f64,
        score2_threshold: f64,
        options: TuClusteringOptions,
        ambiguity_margin: f64,
        fractional_assignment: bool,
    ) -> Vec<ReadAssignment> {
        reads
            .iter()
            .map(|read| {
                let candidates = tus
                    .iter()
                    .enumerate()
                    .filter_map(|(tu_index, candidate_tu)| {
                        assignment_candidate(
                            read,
                            candidate_tu,
                            tu_index,
                            score1_threshold,
                            score2_threshold,
                            options,
                        )
                    })
                    .collect();
                classify_assignment_candidates(
                    candidates,
                    tus,
                    ambiguity_margin,
                    fractional_assignment,
                )
            })
            .collect()
    }

    fn split_members_by_consensus_naive_for_test(
        members: &[usize],
        indices_sorted: &[usize],
        reads: &[ReadRecord],
        score1_threshold: f64,
        score2_threshold: f64,
        options: TuClusteringOptions,
    ) -> Vec<Vec<usize>> {
        let mut families: Vec<Vec<usize>> = Vec::new();
        for local_pos in members_in_support_order(members, indices_sorted, reads) {
            let candidate = &reads[indices_sorted[local_pos]];
            let mut selected_family = None;
            for (family_idx, family) in families.iter().enumerate() {
                let anchor = &reads[indices_sorted[family[0]]];
                if !interval_satisfies_consensus(
                    candidate,
                    anchor.interval,
                    score1_threshold,
                    score2_threshold,
                    options,
                ) {
                    continue;
                }

                let mut prospective = family.clone();
                prospective.push(local_pos);
                let (consensus, _) =
                    family_consensus(&prospective, indices_sorted, reads, score1_threshold);
                if prospective.iter().all(|&member_pos| {
                    interval_satisfies_consensus(
                        &reads[indices_sorted[member_pos]],
                        consensus,
                        score1_threshold,
                        score2_threshold,
                        options,
                    )
                }) {
                    selected_family = Some(family_idx);
                    break;
                }
            }

            if let Some(family_idx) = selected_family {
                families[family_idx].push(local_pos);
            } else {
                families.push(vec![local_pos]);
            }
        }
        families
    }

    #[test]
    fn clustering_result_accessors_preserve_alignment() {
        let reads = vec![
            read("chr", Strand::Plus, 100, 200, "r1"),
            read("chr", Strand::Plus, 101, 201, "r2"),
        ];
        let result = cluster_tus(&reads, 0.95, 0.80).unwrap();

        assert_eq!(result.tus().len(), result.endpoint_stats().len());
        assert_eq!(result.read_to_tu().len(), reads.len());
        assert_eq!(result.tu_for_read(0), result.tus().first());
        assert_eq!(result.tu_for_read(1), result.tus().first());
        assert!(result.tu_for_read(reads.len()).is_none());

        let (tus, read_to_tu, endpoint_stats) = result.into_parts();
        assert_eq!(tus.len(), endpoint_stats.len());
        assert!(read_to_tu.iter().all(|&tu_index| tu_index < tus.len()));
    }

    #[test]
    fn dense_duplicate_locus_does_not_rebuild_the_growing_family() {
        const READ_COUNT: usize = 50_000;
        let reads: Vec<ReadRecord> = (0..READ_COUNT)
            .map(|index| read("chr", Strand::Plus, 100, 1_100, &format!("r{index}")))
            .collect();

        let result = cluster_tus(&reads, 0.95, 0.80).unwrap();

        assert_eq!(result.tus().len(), 1);
        assert_eq!(result.endpoint_stats()[0].support, READ_COUNT);
        assert!(result.read_to_tu().iter().all(|&tu_index| tu_index == 0));
    }

    #[test]
    fn coordinate_identical_final_families_are_coalesced() {
        let reads = vec![
            read("chr", Strand::Plus, 100, 200, "family_a"),
            read("chr", Strand::Plus, 100, 200, "family_b"),
        ];
        let indices_sorted = vec![0, 1];
        let options = TuClusteringOptions::default();
        let families = vec![
            final_consensus_family(vec![0], &indices_sorted, &reads, 0.95, 0.80, options),
            final_consensus_family(vec![1], &indices_sorted, &reads, 0.95, 0.80, options),
        ];

        let coalesced = coalesce_identical_consensus_families(
            families,
            &indices_sorted,
            &reads,
            0.95,
            0.80,
            options,
        );

        assert_eq!(coalesced.len(), 1);
        assert_eq!(coalesced[0].members, vec![0, 1]);
        assert_eq!(coalesced[0].endpoint_stats.support, 2);
        assert_eq!(
            coalesced[0].interval,
            Interval::new(Coord::new(100), Coord::new(200)).unwrap()
        );
    }

    #[test]
    fn public_clustering_coalesces_coordinate_identical_final_families() {
        // These reads split into two anchor-bounded final families whose endpoint
        // modes both converge on [0, 60). The public result must merge them before
        // IDs, assignments, or output records are derived.
        let reads = vec![
            read("chr", Strand::Plus, 10, 60, "r0"),
            read("chr", Strand::Plus, 0, 70, "r1"),
            read("chr", Strand::Plus, 0, 60, "r2"),
            read("chr", Strand::Plus, 0, 60, "r3"),
            read("chr", Strand::Plus, 0, 30, "r4"),
        ];
        let result = cluster_tus_with_options(
            &reads,
            0.5,
            0.5,
            TuClusteringOptions {
                attach_contained_reads: true,
                three_prime_tolerance_bp: 20,
                max_five_prime_delta_bp: Some(20),
            },
        )
        .unwrap();

        assert_eq!(result.tus().len(), 1);
        assert_eq!(
            result.tus()[0].interval,
            Interval::new(Coord::new(0), Coord::new(60)).unwrap()
        );
        assert_eq!(result.endpoint_stats()[0].support, reads.len());
        assert!(result.read_to_tu().iter().all(|&tu_index| tu_index == 0));
    }

    #[test]
    fn exact_assignment_tie_is_ambiguous_not_lexically_resolved() {
        let reads = vec![read("chr", Strand::Plus, 5, 105, "query")];
        let tus = vec![
            tu("chr", Strand::Plus, 10, 110, "TU_B"),
            tu("chr", Strand::Plus, 0, 100, "TU_A"),
        ];

        let assignments = assign_reads_to_tus(
            &reads,
            &tus,
            0.80,
            0.80,
            TuClusteringOptions::default(),
            0.0,
            false,
        )
        .unwrap();
        let assignment = &assignments[0];
        assert_eq!(assignment.status, AssignmentStatus::Ambiguous);
        assert_eq!(
            &tus[assignment.best.unwrap().tu_index].id,
            "TU_A",
            "lexical ordering may stabilize the best field but not resolve the status"
        );
        assert_eq!(&tus[assignment.second_best.unwrap().tu_index].id, "TU_B");
        assert_eq!(assignment.score_margin, Some(0.0));
        assert!(assignment.fractional_assignments.is_empty());
    }

    #[test]
    fn near_tie_respects_configured_ambiguity_margin() {
        let reads = vec![read("chr", Strand::Plus, 2, 102, "query")];
        let tus = vec![
            tu("chr", Strand::Plus, 0, 100, "TU_A"),
            tu("chr", Strand::Plus, 6, 106, "TU_B"),
        ];

        let strict = assign_reads_to_tus(
            &reads,
            &tus,
            0.80,
            0.80,
            TuClusteringOptions::default(),
            0.01,
            false,
        )
        .unwrap();
        let relaxed = assign_reads_to_tus(
            &reads,
            &tus,
            0.80,
            0.80,
            TuClusteringOptions::default(),
            0.05,
            false,
        )
        .unwrap();

        assert_eq!(strict[0].status, AssignmentStatus::Unique);
        assert_eq!(relaxed[0].status, AssignmentStatus::Ambiguous);
        assert!(strict[0].score_margin.unwrap() > 0.01);
        assert!(strict[0].score_margin.unwrap() <= 0.05);
    }

    #[test]
    fn fractional_tie_weights_sum_exactly_to_one() {
        let reads = vec![read("chr", Strand::Plus, 5, 105, "query")];
        let tus = vec![
            tu("chr", Strand::Plus, 0, 110, "TU_A"),
            tu("chr", Strand::Plus, 0, 110, "TU_B"),
            tu("chr", Strand::Plus, 0, 110, "TU_C"),
        ];

        let assignments = assign_reads_to_tus(
            &reads,
            &tus,
            0.80,
            0.80,
            TuClusteringOptions::default(),
            0.0,
            true,
        )
        .unwrap();
        assert_eq!(assignments[0].status, AssignmentStatus::Ambiguous);
        assert_eq!(assignments[0].fractional_assignments.len(), 3);
        assert_eq!(
            assignments[0]
                .fractional_assignments
                .iter()
                .map(|assignment| assignment.weight)
                .sum::<f64>(),
            1.0
        );
    }

    #[test]
    fn contained_only_match_is_partial_and_filtered_candidates_are_unassigned() {
        let reads = vec![read("chr", Strand::Plus, 10, 100, "query")];
        let tus = vec![tu("chr", Strand::Plus, 0, 100, "TU_A")];
        let partial = assign_reads_to_tus(
            &reads,
            &tus,
            0.95,
            0.80,
            TuClusteringOptions::default(),
            0.0,
            false,
        )
        .unwrap();
        let filtered = assign_reads_to_tus(
            &reads,
            &[],
            0.95,
            0.80,
            TuClusteringOptions::default(),
            0.0,
            false,
        )
        .unwrap();

        assert_eq!(partial[0].status, AssignmentStatus::Partial);
        assert!(partial[0].best.unwrap().partial);
        assert_eq!(filtered[0].status, AssignmentStatus::Unassigned);
        assert!(filtered[0].best.is_none());
    }

    #[test]
    fn assignment_evidence_is_deterministic_under_tu_permutations() {
        let reads = vec![read("chr", Strand::Plus, 5, 105, "query")];
        let canonical = vec![
            tu("chr", Strand::Plus, 0, 100, "TU_A"),
            tu("chr", Strand::Plus, 10, 110, "TU_B"),
            tu("chr", Strand::Plus, 20, 120, "TU_C"),
        ];
        let permuted = vec![
            canonical[2].clone(),
            canonical[1].clone(),
            canonical[0].clone(),
        ];

        let summarize = |tus: &[Tu]| {
            let assignment = assign_reads_to_tus(
                &reads,
                tus,
                0.80,
                0.80,
                TuClusteringOptions::default(),
                0.0,
                true,
            )
            .unwrap()
            .remove(0);
            (
                assignment.status,
                assignment.best.map(|candidate| {
                    (
                        tus[candidate.tu_index].id.clone(),
                        candidate.score1,
                        candidate.score2,
                        candidate.five_prime_delta_bp,
                        candidate.three_prime_delta_bp,
                    )
                }),
                assignment.second_best.map(|candidate| {
                    (
                        tus[candidate.tu_index].id.clone(),
                        candidate.score1,
                        candidate.score2,
                        candidate.five_prime_delta_bp,
                        candidate.three_prime_delta_bp,
                    )
                }),
            )
        };

        assert_eq!(summarize(&canonical), summarize(&permuted));
    }

    #[test]
    fn assignment_interval_index_skips_distant_tus() {
        let tus: Vec<Tu> = (0..20_000u32)
            .map(|index| {
                let start = index * 10;
                tu("chr", Strand::Plus, start, start + 5, &format!("TU{index}"))
            })
            .collect();
        let index = TuIntervalIndex::build((0..tus.len()).collect(), &tus);
        let query = Interval::new(Coord::new(100_002), Coord::new(100_004)).unwrap();
        let mut visited = Vec::new();
        index.for_each_overlap(query, &tus, |tu_index| visited.push(tu_index));

        assert_eq!(visited, vec![10_000]);
        assert!(
            visited.len() * 1_000 < tus.len(),
            "a local assignment query must not degrade to scanning the TU collection"
        );
    }

    #[test]
    fn indexed_assignment_matches_naive_with_unsorted_tus() {
        let reads = vec![
            read("chr2", Strand::Minus, 190, 260, "r0"),
            read("chr1", Strand::Plus, 95, 205, "r1"),
            read("chr1", Strand::Minus, 100, 200, "r2"),
            read("chr1", Strand::Plus, 500, 550, "r3"),
        ];
        let tus = vec![
            tu("chr1", Strand::Plus, 490, 560, "TU_D"),
            tu("chr2", Strand::Minus, 200, 270, "TU_C"),
            tu("chr1", Strand::Plus, 90, 200, "TU_B"),
            tu("chr1", Strand::Minus, 90, 210, "TU_E"),
            tu("chr1", Strand::Plus, 100, 210, "TU_A"),
            tu("chr3", Strand::Plus, 0, 1_000, "TU_F"),
        ];
        let options = TuClusteringOptions {
            max_five_prime_delta_bp: Some(20),
            ..TuClusteringOptions::default()
        };

        let indexed = assign_reads_to_tus(&reads, &tus, 0.80, 0.70, options, 0.05, true).unwrap();
        let naive =
            assign_reads_to_tus_naive_for_test(&reads, &tus, 0.80, 0.70, options, 0.05, true);
        assert_eq!(indexed, naive);
    }

    fn assignments_by_read_id(
        reads: &[ReadRecord],
        result: &TuClusteringResult,
    ) -> BTreeMap<String, (String, String, Interval)> {
        reads
            .iter()
            .enumerate()
            .map(|(read_idx, read)| {
                let tu = &result.tus()[result.read_to_tu()[read_idx]];
                (
                    read.id.clone(),
                    (
                        tu.id.clone(),
                        reads[tu.rep_read_index].id.clone(),
                        tu.interval,
                    ),
                )
            })
            .collect()
    }

    fn mirror_reads(reads: &[ReadRecord], mirror_axis: u32) -> Vec<ReadRecord> {
        reads
            .iter()
            .map(|read| {
                assert!(read.interval.end().get() <= mirror_axis);
                ReadRecord {
                    contig: read.contig.clone(),
                    strand: match read.strand {
                        Strand::Plus => Strand::Minus,
                        Strand::Minus => Strand::Plus,
                        Strand::Unknown => Strand::Unknown,
                    },
                    interval: Interval::new(
                        Coord::new(mirror_axis - read.interval.end().get()),
                        Coord::new(mirror_axis - read.interval.start().get()),
                    )
                    .unwrap(),
                    id: read.id.clone(),
                }
            })
            .collect()
    }

    fn representative_ids_by_read_id(
        reads: &[ReadRecord],
        result: &TuClusteringResult,
    ) -> BTreeMap<String, String> {
        assignments_by_read_id(reads, result)
            .into_iter()
            .map(|(read_id, (_, representative_id, _))| (read_id, representative_id))
            .collect()
    }

    #[test]
    fn default_score2_avoids_overmerging_short_contained_reads() {
        let reads = vec![
            read("chr1", Strand::Plus, 120, 180, "r3"),
            read("chr1", Strand::Plus, 100, 200, "r1"),
            read("chr1", Strand::Plus, 101, 201, "r2"),
            read("chr1", Strand::Plus, 300, 400, "r4"),
            read("chr1", Strand::Plus, 301, 401, "r5"),
            read("chr1", Strand::Plus, 320, 360, "r6"),
            read("chr1", Strand::Minus, 100, 200, "r7"),
        ];

        let result = cluster_tus(&reads, 0.95, 0.99).unwrap();
        assert_eq!(result.tus().len(), 5);

        assert_eq!(
            result.tus()[0],
            Tu {
                id: "TU000001".to_owned(),
                contig: "chr1".to_owned(),
                strand: Strand::Plus,
                interval: Interval::new(Coord::new(100), Coord::new(200)).unwrap(),
                rep_read_index: 1,
            }
        );
        assert_eq!(
            result.tus()[1],
            Tu {
                id: "TU000002".to_owned(),
                contig: "chr1".to_owned(),
                strand: Strand::Plus,
                interval: Interval::new(Coord::new(120), Coord::new(180)).unwrap(),
                rep_read_index: 0,
            }
        );
        assert_eq!(
            result.tus()[2],
            Tu {
                id: "TU000003".to_owned(),
                contig: "chr1".to_owned(),
                strand: Strand::Plus,
                interval: Interval::new(Coord::new(300), Coord::new(400)).unwrap(),
                rep_read_index: 3,
            }
        );
        assert_eq!(
            result.tus()[3],
            Tu {
                id: "TU000004".to_owned(),
                contig: "chr1".to_owned(),
                strand: Strand::Plus,
                interval: Interval::new(Coord::new(320), Coord::new(360)).unwrap(),
                rep_read_index: 5,
            }
        );
        assert_eq!(
            result.tus()[4],
            Tu {
                id: "TU000005".to_owned(),
                contig: "chr1".to_owned(),
                strand: Strand::Minus,
                interval: Interval::new(Coord::new(100), Coord::new(200)).unwrap(),
                rep_read_index: 6,
            }
        );

        let mut by_id: HashMap<String, String> = HashMap::new();
        for (read_idx, read) in reads.iter().enumerate() {
            let tu_id = &result.tus()[result.read_to_tu()[read_idx]].id;
            by_id.insert(read.id.clone(), tu_id.clone());
        }

        assert_eq!(by_id["r1"], "TU000001");
        assert_eq!(by_id["r2"], "TU000001");
        assert_eq!(by_id["r3"], "TU000002");
        assert_eq!(by_id["r4"], "TU000003");
        assert_eq!(by_id["r5"], "TU000003");
        assert_eq!(by_id["r6"], "TU000004");
        assert_eq!(by_id["r7"], "TU000005");
    }

    #[test]
    fn score2_can_attach_near_full_length_reads_when_threshold_is_relaxed() {
        let reads = vec![
            read("chr1", Strand::Plus, 108, 200, "r3"),
            read("chr1", Strand::Plus, 100, 200, "r1"),
            read("chr1", Strand::Plus, 101, 201, "r2"),
        ];

        let result = cluster_tus(&reads, 0.95, 0.90).unwrap();
        assert_eq!(result.tus().len(), 1);

        let mut by_id: HashMap<String, String> = HashMap::new();
        for (read_idx, read) in reads.iter().enumerate() {
            let tu_id = &result.tus()[result.read_to_tu()[read_idx]].id;
            by_id.insert(read.id.clone(), tu_id.clone());
        }

        assert_eq!(by_id["r1"], "TU000001");
        assert_eq!(by_id["r2"], "TU000001");
        assert_eq!(by_id["r3"], "TU000001");
    }

    #[test]
    fn default_three_prime_tolerance_merges_plus_strand_near_matches() {
        let reads = vec![
            read("chr1", Strand::Plus, 100, 205, "parent"),
            read("chr1", Strand::Plus, 112, 210, "child"),
        ];

        let result = cluster_tus(&reads, 0.95, 0.60).unwrap();
        assert_eq!(result.tus().len(), 1);
    }

    #[test]
    fn default_three_prime_tolerance_merges_minus_strand_near_matches() {
        let reads = vec![
            read("chr1", Strand::Minus, 100, 205, "parent"),
            read("chr1", Strand::Minus, 108, 210, "child"),
        ];

        let result = cluster_tus(&reads, 0.95, 0.60).unwrap();
        assert_eq!(result.tus().len(), 1);
    }

    #[test]
    fn three_prime_gap_beyond_tolerance_stays_split() {
        let reads = vec![
            read("chr1", Strand::Plus, 100, 205, "parent"),
            read("chr1", Strand::Plus, 120, 218, "child"),
        ];

        let result = cluster_tus(&reads, 0.95, 0.60).unwrap();
        assert_eq!(result.tus().len(), 2);
    }

    #[test]
    fn second_pass_attachment_is_independent_of_genomic_and_input_order() {
        let forward = vec![
            read("chr1", Strand::Plus, 0, 95, "shorter"),
            read("chr1", Strand::Plus, 5, 105, "longer"),
        ];
        let reversed = vec![forward[1].clone(), forward[0].clone()];
        let options = TuClusteringOptions {
            three_prime_tolerance_bp: 12,
            ..TuClusteringOptions::default()
        };

        let forward_result = cluster_tus_with_options(&forward, 0.95, 0.90, options).unwrap();
        let reversed_result = cluster_tus_with_options(&reversed, 0.95, 0.90, options).unwrap();

        assert_eq!(forward_result.tus().len(), 1);
        assert_eq!(reversed_result.tus().len(), 1);
        assert_eq!(
            assignments_by_read_id(&forward, &forward_result),
            assignments_by_read_id(&reversed, &reversed_result)
        );
        assert_eq!(forward[forward_result.tus()[0].rep_read_index].id, "longer");
    }

    #[test]
    fn mirrored_minus_strand_attachment_matches_plus_strand() {
        let plus = vec![
            read("chr1", Strand::Plus, 0, 95, "shorter"),
            read("chr1", Strand::Plus, 5, 105, "longer"),
        ];
        let minus = mirror_reads(&plus, 105);
        let options = TuClusteringOptions {
            three_prime_tolerance_bp: 12,
            ..TuClusteringOptions::default()
        };

        let plus_result = cluster_tus_with_options(&plus, 0.95, 0.90, options).unwrap();
        let minus_result = cluster_tus_with_options(&minus, 0.95, 0.90, options).unwrap();

        assert_eq!(plus_result.tus().len(), 1);
        assert_eq!(minus_result.tus().len(), 1);
        assert_eq!(
            representative_ids_by_read_id(&plus, &plus_result),
            representative_ids_by_read_id(&minus, &minus_result)
        );
    }

    #[test]
    fn five_prime_override_absorbs_truncated_but_not_overhanging_reads() {
        let options = TuClusteringOptions {
            three_prime_tolerance_bp: 12,
            max_five_prime_delta_bp: Some(1),
            ..TuClusteringOptions::default()
        };

        // A 5'-anchored fragment terminating inside the longer read is a
        // truncated molecule of it: the override pools the pair and
        // containment absorbs the fragment without moving the boundary.
        let truncated = vec![
            read("chr1", Strand::Plus, 0, 1000, "longer"),
            read("chr1", Strand::Plus, 1, 101, "shorter"),
        ];
        let result = cluster_tus_with_options(&truncated, 0.95, 0.99, options).unwrap();
        assert_eq!(result.tus().len(), 1);
        assert_eq!(
            result.tus()[0].interval,
            Interval::new(Coord::new(0), Coord::new(1000)).unwrap()
        );
        assert_eq!(result.endpoint_stats()[0].support, 2);

        // A read overhanging the 3' end beyond the tolerance is not a
        // fragment: the one-sided attachment gate rejects it even though the
        // 5' override matches.
        let overhanging = vec![
            read("chr1", Strand::Plus, 0, 200, "longer"),
            read("chr1", Strand::Plus, 1, 220, "overhang"),
        ];
        let result = cluster_tus_with_options(&overhanging, 0.95, 0.99, options).unwrap();
        assert_eq!(result.tus().len(), 2);
    }

    #[test]
    fn five_prime_override_uses_absolute_strand_aware_delta() {
        let plus = vec![
            read("chr1", Strand::Plus, 0, 95, "shorter"),
            read("chr1", Strand::Plus, 5, 105, "longer"),
        ];
        let minus = mirror_reads(&plus, 105);

        for reads in [&plus, &minus] {
            let result = cluster_tus_with_options(
                reads,
                0.95,
                0.99,
                TuClusteringOptions {
                    three_prime_tolerance_bp: 12,
                    max_five_prime_delta_bp: Some(0),
                    ..TuClusteringOptions::default()
                },
            )
            .unwrap();

            assert_eq!(result.tus().len(), 2);
        }
    }

    #[test]
    fn same_genomic_start_containment_attaches_on_both_strands() {
        for strand in [Strand::Plus, Strand::Minus] {
            let reads = vec![
                read("chr1", strand, 0, 80, "shorter"),
                read("chr1", strand, 0, 100, "longer"),
            ];
            let result = cluster_tus_with_options(
                &reads,
                0.95,
                0.80,
                TuClusteringOptions {
                    three_prime_tolerance_bp: 20,
                    ..TuClusteringOptions::default()
                },
            )
            .unwrap();

            assert_eq!(result.tus().len(), 1, "strand={strand:?}");
            assert_eq!(reads[result.tus()[0].rep_read_index].id, "longer");
        }
    }

    #[test]
    fn equal_length_attachment_chooses_representative_deterministically() {
        let forward = vec![
            read("chr1", Strand::Plus, 0, 100, "z_read"),
            read("chr1", Strand::Plus, 5, 105, "a_read"),
        ];
        let reversed = vec![forward[1].clone(), forward[0].clone()];

        for reads in [&forward, &reversed] {
            let result = cluster_tus_with_options(
                reads,
                0.99,
                0.95,
                TuClusteringOptions {
                    three_prime_tolerance_bp: 5,
                    ..TuClusteringOptions::default()
                },
            )
            .unwrap();
            assert_eq!(result.tus().len(), 1);
            assert_eq!(reads[result.tus()[0].rep_read_index].id, "a_read");
        }
    }

    #[test]
    fn best_parent_ties_use_representative_id_independent_of_input_order() {
        let canonical = vec![
            read("chr1", Strand::Plus, 90, 190, "a_parent"),
            read("chr1", Strand::Plus, 100, 200, "z_child"),
            read("chr1", Strand::Plus, 110, 210, "b_parent"),
        ];
        let permuted = vec![
            canonical[2].clone(),
            canonical[1].clone(),
            canonical[0].clone(),
        ];
        let options = TuClusteringOptions {
            three_prime_tolerance_bp: 10,
            ..TuClusteringOptions::default()
        };

        for reads in [&canonical, &permuted] {
            let result = cluster_tus_with_options(reads, 0.95, 0.90, options).unwrap();
            let assignments = representative_ids_by_read_id(reads, &result);

            assert_eq!(result.tus().len(), 2);
            assert_eq!(assignments["z_child"], "a_parent");
            assert_eq!(assignments["a_parent"], "a_parent");
            assert_eq!(assignments["b_parent"], "b_parent");
        }
    }

    #[test]
    fn clustering_is_identical_across_rayon_thread_counts() {
        let reads = vec![
            read("chr2", Strand::Minus, 500, 600, "m_long"),
            read("chr1", Strand::Plus, 0, 95, "p_short"),
            read("chr1", Strand::Plus, 5, 105, "p_long"),
            read("chr2", Strand::Minus, 510, 605, "m_short"),
            read("chr3", Strand::Plus, 20, 80, "isolated"),
        ];
        let run = || cluster_tus(&reads, 0.95, 0.90).unwrap();
        let single = rayon::ThreadPoolBuilder::new()
            .num_threads(1)
            .build()
            .unwrap()
            .install(run);
        let multiple = rayon::ThreadPoolBuilder::new()
            .num_threads(4)
            .build()
            .unwrap()
            .install(run);

        assert_eq!(
            assignments_by_read_id(&reads, &single),
            assignments_by_read_id(&reads, &multiple)
        );
        assert_eq!(single.tus(), multiple.tus());
        assert_eq!(single.read_to_tu(), multiple.read_to_tu());
    }

    #[test]
    fn score1_daisy_chain_is_split_by_the_final_consensus_invariant() {
        let reads = vec![
            read("chr1", Strand::Plus, 0, 1000, "left"),
            read("chr1", Strand::Plus, 25, 1025, "bridge"),
            read("chr1", Strand::Plus, 50, 1050, "right"),
        ];
        let options = TuClusteringOptions {
            attach_contained_reads: false,
            ..TuClusteringOptions::default()
        };

        let result = cluster_tus_with_options(&reads, 0.95, 0.99, options).unwrap();

        assert_eq!(result.tus().len(), 2);
        assert_eq!(result.endpoint_stats().len(), result.tus().len());
        for (read_idx, read) in reads.iter().enumerate() {
            let consensus = result.tus()[result.read_to_tu()[read_idx]].interval;
            assert!(
                score1_interval(read.interval, consensus) >= 0.95,
                "{} must qualify directly against {consensus:?}",
                read.id
            );
        }
        assert_eq!(
            result.tus()[0].interval,
            Interval::new(Coord::new(0), Coord::new(1000)).unwrap()
        );
        assert_eq!(result.endpoint_stats()[0].support, 2);
        assert_eq!(result.endpoint_stats()[0].five_prime_spread_bp(), 25);
        assert_eq!(result.endpoint_stats()[0].three_prime_spread_bp(), 25);
    }

    #[test]
    fn adjacent_endpoint_windows_do_not_chain_without_a_matching_anchor() {
        let reads = vec![
            read("chr1", Strand::Plus, 0, 100, "r0"),
            read("chr1", Strand::Plus, 0, 112, "r12"),
            read("chr1", Strand::Plus, 0, 124, "r24"),
            read("chr1", Strand::Plus, 0, 136, "r36"),
        ];
        let options = TuClusteringOptions {
            three_prime_tolerance_bp: 12,
            ..TuClusteringOptions::default()
        };

        let result = cluster_tus_with_options(&reads, 0.95, 0.80, options).unwrap();

        assert_eq!(result.tus().len(), 2);
        assert_eq!(
            result
                .endpoint_stats()
                .iter()
                .map(|stats| stats.support)
                .collect::<Vec<_>>(),
            vec![2, 2]
        );
        assert!(result
            .endpoint_stats()
            .iter()
            .all(|stats| stats.three_prime_spread_bp() <= 12));
    }

    #[test]
    fn low_support_read_cannot_bridge_two_high_support_modes() {
        let mut reads = Vec::new();
        for id in ["left1", "left2", "left3"] {
            reads.push(read("chr1", Strand::Plus, 0, 1000, id));
        }
        reads.push(read("chr1", Strand::Plus, 25, 1025, "bridge"));
        for id in ["right1", "right2", "right3"] {
            reads.push(read("chr1", Strand::Plus, 50, 1050, id));
        }

        let result = cluster_tus_with_options(
            &reads,
            0.95,
            0.99,
            TuClusteringOptions {
                attach_contained_reads: false,
                ..TuClusteringOptions::default()
            },
        )
        .unwrap();

        assert_eq!(result.tus().len(), 2);
        let mut support: Vec<usize> = result
            .endpoint_stats()
            .iter()
            .map(|stats| stats.support)
            .collect();
        support.sort_unstable();
        assert_eq!(support, vec![3, 4]);
        assert_ne!(result.read_to_tu()[0], result.read_to_tu()[4]);
    }

    #[test]
    fn emitted_boundary_is_consensus_and_longest_read_is_only_the_example() {
        let reads = vec![
            read("chr1", Strand::Plus, 0, 100, "longest"),
            read("chr1", Strand::Plus, 10, 100, "mode_a"),
            read("chr1", Strand::Plus, 10, 100, "mode_b"),
        ];

        let result = cluster_tus_with_options(
            &reads,
            0.80,
            0.99,
            TuClusteringOptions {
                attach_contained_reads: false,
                ..TuClusteringOptions::default()
            },
        )
        .unwrap();

        assert_eq!(result.tus().len(), 1);
        assert_eq!(reads[result.tus()[0].rep_read_index].id, "longest");
        assert_eq!(
            result.tus()[0].interval,
            Interval::new(Coord::new(10), Coord::new(100)).unwrap()
        );
        assert_eq!(result.endpoint_stats()[0].five_prime_consensus, 10);
        assert_eq!(result.endpoint_stats()[0].support, 3);
    }

    #[test]
    fn endpoint_mode_preserves_repeated_full_boundaries_amid_distinct_truncations() {
        let plus_five_prime = [100, 100, 120, 180, 250];
        let plus_three_prime = [500, 500, 470, 450, 430];
        assert_eq!(
            modal_consensus_endpoint(&plus_five_prime, Strand::Plus, true),
            100
        );
        assert_eq!(
            modal_consensus_endpoint(&plus_three_prime, Strand::Plus, false),
            500
        );

        let minus_five_prime = [500, 500, 480, 420, 350];
        let minus_three_prime = [100, 100, 130, 150, 170];
        assert_eq!(
            modal_consensus_endpoint(&minus_five_prime, Strand::Minus, true),
            500
        );
        assert_eq!(
            modal_consensus_endpoint(&minus_three_prime, Strand::Minus, false),
            100
        );
    }

    #[test]
    fn bounded_consensus_families_are_permutation_and_mirror_deterministic() {
        let canonical = vec![
            read("chr1", Strand::Plus, 0, 1000, "left"),
            read("chr1", Strand::Plus, 25, 1025, "bridge"),
            read("chr1", Strand::Plus, 50, 1050, "right"),
        ];
        let permuted = vec![
            canonical[2].clone(),
            canonical[0].clone(),
            canonical[1].clone(),
        ];
        let mirrored = mirror_reads(&canonical, 1050);
        let options = TuClusteringOptions {
            attach_contained_reads: false,
            ..TuClusteringOptions::default()
        };

        let canonical_result = cluster_tus_with_options(&canonical, 0.95, 0.99, options).unwrap();
        let permuted_result = cluster_tus_with_options(&permuted, 0.95, 0.99, options).unwrap();
        let mirrored_result = cluster_tus_with_options(&mirrored, 0.95, 0.99, options).unwrap();

        assert_eq!(
            representative_ids_by_read_id(&canonical, &canonical_result),
            representative_ids_by_read_id(&permuted, &permuted_result)
        );
        assert_eq!(
            representative_ids_by_read_id(&canonical, &canonical_result),
            representative_ids_by_read_id(&mirrored, &mirrored_result)
        );
        for (read_idx, read) in canonical.iter().enumerate() {
            let mirrored_idx = mirrored
                .iter()
                .position(|candidate| candidate.id == read.id)
                .unwrap();
            let plus_interval =
                canonical_result.tus()[canonical_result.read_to_tu()[read_idx]].interval;
            let minus_interval =
                mirrored_result.tus()[mirrored_result.read_to_tu()[mirrored_idx]].interval;
            assert_eq!(
                minus_interval.start().get(),
                1050 - plus_interval.end().get()
            );
            assert_eq!(
                minus_interval.end().get(),
                1050 - plus_interval.start().get()
            );
        }
    }

    proptest! {
        #[test]
        fn incremental_endpoint_consensus_matches_batch_modes(
            coordinates in prop::collection::vec(0u32..1_000, 1..80),
            reverse in any::<bool>(),
            five_prime in any::<bool>(),
        ) {
            let strand = if reverse { Strand::Minus } else { Strand::Plus };
            let prefer_lower = matches!(
                (strand, five_prime),
                (Strand::Plus, true) | (Strand::Minus, false)
            );
            let mut state = EndpointConsensusState::new(strand, prefer_lower);
            let mut observed = Vec::new();

            for &coordinate in &coordinates {
                state.insert(coordinate);
                observed.push(coordinate);
                prop_assert_eq!(
                    state.mode(),
                    modal_consensus_endpoint(&observed, strand, five_prime)
                );
            }

            for &coordinate in coordinates.iter().rev().take(coordinates.len() - 1) {
                state.remove(coordinate);
                observed.pop();
                prop_assert_eq!(
                    state.mode(),
                    modal_consensus_endpoint(&observed, strand, five_prime)
                );
            }
        }

        #[test]
        fn incremental_family_splitting_matches_batch_reference(
            interval_specs in prop::collection::vec((0u32..500, 1u32..200), 1..24),
            reverse in any::<bool>(),
            attach_contained_reads in any::<bool>(),
        ) {
            let strand = if reverse { Strand::Minus } else { Strand::Plus };
            let reads: Vec<ReadRecord> = interval_specs
                .into_iter()
                .enumerate()
                .map(|(index, (start, len))| {
                    read("chr", strand, start, start + len, &format!("r{index:03}"))
                })
                .collect();
            let mut indices_sorted: Vec<usize> = (0..reads.len()).collect();
            indices_sorted.sort_unstable_by(|&a, &b| cmp_read_idx(&reads, a, b));
            let members: Vec<usize> = (0..reads.len()).collect();
            let options = TuClusteringOptions {
                attach_contained_reads,
                three_prime_tolerance_bp: 20,
                max_five_prime_delta_bp: Some(15),
            };

            let incremental = split_members_by_consensus(
                &members,
                &indices_sorted,
                &reads,
                0.75,
                0.60,
                options,
            );
            let reference = split_members_by_consensus_naive_for_test(
                &members,
                &indices_sorted,
                &reads,
                0.75,
                0.60,
                options,
            );
            prop_assert_eq!(incremental, reference);
        }

        #[test]
        fn indexed_assignment_matches_naive_randomized_unordered_inputs(
            read_specs in prop::collection::vec(
                (0u8..2, 0u8..2, 0u32..500, 0u32..100),
                0..20,
            ),
            tu_specs in prop::collection::vec(
                (0u8..2, 0u8..2, 0u32..500, 0u32..100),
                0..20,
            ),
            attach_contained_reads in any::<bool>(),
            fractional_assignment in any::<bool>(),
        ) {
            let choose_contig = |key| if key == 0 { "chr1" } else { "chr2" };
            let choose_strand = |key| if key == 0 { Strand::Plus } else { Strand::Minus };
            let reads: Vec<ReadRecord> = read_specs
                .into_iter()
                .enumerate()
                .map(|(index, (contig, strand, start, len))| {
                    read(
                        choose_contig(contig),
                        choose_strand(strand),
                        start,
                        start + len,
                        &format!("r{index}"),
                    )
                })
                .collect();
            let tus: Vec<Tu> = tu_specs
                .into_iter()
                .enumerate()
                .map(|(index, (contig, strand, start, len))| {
                    tu(
                        choose_contig(contig),
                        choose_strand(strand),
                        start,
                        start + len,
                        &format!("TU{index:06}"),
                    )
                })
                .collect();
            let options = TuClusteringOptions {
                attach_contained_reads,
                three_prime_tolerance_bp: 20,
                max_five_prime_delta_bp: Some(15),
            };

            let indexed = assign_reads_to_tus(
                &reads,
                &tus,
                0.75,
                0.60,
                options,
                0.05,
                fractional_assignment,
            ).unwrap();
            let naive = assign_reads_to_tus_naive_for_test(
                &reads,
                &tus,
                0.75,
                0.60,
                options,
                0.05,
                fractional_assignment,
            );

            prop_assert_eq!(indexed, naive);
        }

        #[test]
        fn random_input_permutations_preserve_assignments(keys in prop::collection::vec(any::<u64>(), 8)) {
            let canonical = vec![
                read("chr1", Strand::Plus, 0, 95, "p_short"),
                read("chr1", Strand::Plus, 5, 105, "p_long"),
                read("chr1", Strand::Plus, 300, 380, "p_same_start_short"),
                read("chr1", Strand::Plus, 300, 400, "p_same_start_long"),
                read("chr2", Strand::Minus, 100, 195, "m_short"),
                read("chr2", Strand::Minus, 90, 190, "m_long"),
                read("chr3", Strand::Plus, 0, 1000, "far_three_prime_long"),
                read("chr3", Strand::Plus, 1, 101, "far_three_prime_short"),
            ];
            let options = TuClusteringOptions {
                three_prime_tolerance_bp: 20,
                max_five_prime_delta_bp: Some(1),
                ..TuClusteringOptions::default()
            };
            let canonical_result = cluster_tus_with_options(&canonical, 0.95, 0.80, options).unwrap();
            let expected = assignments_by_read_id(&canonical, &canonical_result);

            let mut order: Vec<usize> = (0..canonical.len()).collect();
            order.sort_by_key(|&idx| (keys[idx], idx));
            let permuted: Vec<ReadRecord> = order
                .into_iter()
                .map(|idx| canonical[idx].clone())
                .collect();
            let permuted_result = cluster_tus_with_options(&permuted, 0.95, 0.80, options).unwrap();

            prop_assert_eq!(assignments_by_read_id(&permuted, &permuted_result), expected);
        }

        #[test]
        fn plus_minus_mirroring_preserves_assignments(
            intervals in prop::collection::vec((0u32..800, 1u32..150), 1..8)
        ) {
            let plus: Vec<ReadRecord> = intervals
                .into_iter()
                .enumerate()
                .map(|(idx, (start, len))| {
                    read("chr1", Strand::Plus, start, start + len, &format!("r{idx}"))
                })
                .collect();
            let minus = mirror_reads(&plus, 1000);
            let options = TuClusteringOptions {
                three_prime_tolerance_bp: 20,
                max_five_prime_delta_bp: Some(10),
                ..TuClusteringOptions::default()
            };
            let plus_result = cluster_tus_with_options(&plus, 0.95, 0.70, options).unwrap();
            let minus_result = cluster_tus_with_options(&minus, 0.95, 0.70, options).unwrap();

            prop_assert_eq!(
                representative_ids_by_read_id(&plus, &plus_result),
                representative_ids_by_read_id(&minus, &minus_result)
            );
        }
    }

    #[test]
    fn five_prime_delta_override_can_merge_plus_strand_near_matches() {
        let reads = vec![
            read("chr1", Strand::Plus, 100, 205, "parent"),
            read("chr1", Strand::Plus, 150, 210, "child"),
        ];

        let default_result = cluster_tus(&reads, 0.95, 0.60).unwrap();
        assert_eq!(default_result.tus().len(), 2);

        let relaxed_result = cluster_tus_with_options(
            &reads,
            0.95,
            0.60,
            TuClusteringOptions {
                max_five_prime_delta_bp: Some(50),
                ..TuClusteringOptions::default()
            },
        )
        .unwrap();
        assert_eq!(relaxed_result.tus().len(), 1);
    }

    #[test]
    fn five_prime_delta_override_can_merge_minus_strand_near_matches() {
        let reads = vec![
            read("chr1", Strand::Minus, 100, 205, "parent"),
            read("chr1", Strand::Minus, 105, 155, "child"),
        ];

        let default_result = cluster_tus(&reads, 0.95, 0.60).unwrap();
        assert_eq!(default_result.tus().len(), 2);

        let relaxed_result = cluster_tus_with_options(
            &reads,
            0.95,
            0.60,
            TuClusteringOptions {
                max_five_prime_delta_bp: Some(50),
                ..TuClusteringOptions::default()
            },
        )
        .unwrap();
        assert_eq!(relaxed_result.tus().len(), 1);
    }

    #[test]
    fn five_prime_delta_override_respects_the_requested_cap() {
        let reads = vec![
            read("chr1", Strand::Plus, 100, 205, "parent"),
            read("chr1", Strand::Plus, 150, 210, "child"),
        ];

        let result = cluster_tus_with_options(
            &reads,
            0.95,
            0.60,
            TuClusteringOptions {
                max_five_prime_delta_bp: Some(40),
                ..TuClusteringOptions::default()
            },
        )
        .unwrap();
        assert_eq!(result.tus().len(), 2);
    }

    #[test]
    fn five_prime_delta_override_does_not_block_normal_score2_merges() {
        let reads = vec![
            read("chr1", Strand::Plus, 100, 205, "parent"),
            read("chr1", Strand::Plus, 120, 210, "child"),
        ];

        let result = cluster_tus_with_options(
            &reads,
            0.95,
            0.60,
            TuClusteringOptions {
                max_five_prime_delta_bp: Some(10),
                ..TuClusteringOptions::default()
            },
        )
        .unwrap();
        assert_eq!(result.tus().len(), 1);
    }

    #[test]
    fn score1_only_keeps_contained_reads_as_separate_tus() {
        let reads = vec![
            read("chr1", Strand::Plus, 120, 180, "r3"),
            read("chr1", Strand::Plus, 100, 200, "r1"),
            read("chr1", Strand::Plus, 101, 201, "r2"),
            read("chr1", Strand::Plus, 300, 400, "r4"),
            read("chr1", Strand::Plus, 301, 401, "r5"),
            read("chr1", Strand::Plus, 320, 360, "r6"),
            read("chr1", Strand::Minus, 100, 200, "r7"),
        ];

        let result = cluster_tus_with_options(
            &reads,
            0.95,
            0.99,
            TuClusteringOptions {
                attach_contained_reads: false,
                ..TuClusteringOptions::default()
            },
        )
        .unwrap();
        assert_eq!(result.tus().len(), 5);

        assert_eq!(
            result.tus()[0],
            Tu {
                id: "TU000001".to_owned(),
                contig: "chr1".to_owned(),
                strand: Strand::Plus,
                interval: Interval::new(Coord::new(100), Coord::new(200)).unwrap(),
                rep_read_index: 1,
            }
        );
        assert_eq!(
            result.tus()[1],
            Tu {
                id: "TU000002".to_owned(),
                contig: "chr1".to_owned(),
                strand: Strand::Plus,
                interval: Interval::new(Coord::new(120), Coord::new(180)).unwrap(),
                rep_read_index: 0,
            }
        );
        assert_eq!(
            result.tus()[2],
            Tu {
                id: "TU000003".to_owned(),
                contig: "chr1".to_owned(),
                strand: Strand::Plus,
                interval: Interval::new(Coord::new(300), Coord::new(400)).unwrap(),
                rep_read_index: 3,
            }
        );
        assert_eq!(
            result.tus()[3],
            Tu {
                id: "TU000004".to_owned(),
                contig: "chr1".to_owned(),
                strand: Strand::Plus,
                interval: Interval::new(Coord::new(320), Coord::new(360)).unwrap(),
                rep_read_index: 5,
            }
        );
        assert_eq!(
            result.tus()[4],
            Tu {
                id: "TU000005".to_owned(),
                contig: "chr1".to_owned(),
                strand: Strand::Minus,
                interval: Interval::new(Coord::new(100), Coord::new(200)).unwrap(),
                rep_read_index: 6,
            }
        );

        let mut by_id: HashMap<String, String> = HashMap::new();
        for (read_idx, read) in reads.iter().enumerate() {
            let tu_id = &result.tus()[result.read_to_tu()[read_idx]].id;
            by_id.insert(read.id.clone(), tu_id.clone());
        }

        assert_eq!(by_id["r1"], "TU000001");
        assert_eq!(by_id["r2"], "TU000001");
        assert_eq!(by_id["r3"], "TU000002");
        assert_eq!(by_id["r4"], "TU000003");
        assert_eq!(by_id["r5"], "TU000003");
        assert_eq!(by_id["r6"], "TU000004");
        assert_eq!(by_id["r7"], "TU000005");
    }

    #[test]
    fn five_prime_degradation_ladder_collapses_into_one_tu() {
        // Direct-RNA style locus: a full-length boundary mode plus a dense 5'
        // truncation ladder sharing the 3' end. The overlap-over-longer pass
        // chains the ladder to the full-length family; the final consensus
        // must retain every 3'-anchored contained fragment instead of packing
        // the ladder into separate low-support TUs.
        let reads = vec![
            read("chr1", Strand::Plus, 0, 3000, "full1"),
            read("chr1", Strand::Plus, 0, 3000, "full2"),
            read("chr1", Strand::Plus, 0, 3000, "full3"),
            read("chr1", Strand::Plus, 450, 3004, "trunc450"),
            read("chr1", Strand::Plus, 832, 2998, "trunc832"),
            read("chr1", Strand::Plus, 1157, 3006, "trunc1157"),
            read("chr1", Strand::Plus, 1433, 3002, "trunc1433"),
            read("chr1", Strand::Plus, 1668, 2996, "trunc1668"),
        ];

        let result = cluster_tus(&reads, 0.95, 0.80).unwrap();

        assert_eq!(result.tus().len(), 1);
        assert_eq!(
            result.tus()[0].interval,
            Interval::new(Coord::new(0), Coord::new(3000)).unwrap()
        );
        assert_eq!(result.endpoint_stats()[0].support, reads.len());
        assert!(result.read_to_tu().iter().all(|&tu_index| tu_index == 0));

        let mirrored = mirror_reads(&reads, 3006);
        let mirrored_result = cluster_tus(&mirrored, 0.95, 0.80).unwrap();
        assert_eq!(mirrored_result.tus().len(), 1);
        assert_eq!(mirrored_result.endpoint_stats()[0].support, reads.len());
        assert_eq!(
            mirrored_result.tus()[0].interval,
            Interval::new(Coord::new(6), Coord::new(3006)).unwrap()
        );
    }

    #[test]
    fn contained_reads_with_distant_three_prime_stay_separate() {
        // Containment absorption only applies within a pooled component, and
        // pooling still demands 3' agreement or a score1 chain: a nested
        // population terminating 500 bp early never pools with the covering
        // family, so it stays a separate TU.
        let reads = vec![
            read("chr1", Strand::Plus, 0, 3000, "full1"),
            read("chr1", Strand::Plus, 0, 3000, "full2"),
            read("chr1", Strand::Plus, 500, 2500, "nested1"),
            read("chr1", Strand::Plus, 500, 2500, "nested2"),
        ];

        let result = cluster_tus(&reads, 0.95, 0.80).unwrap();

        assert_eq!(result.tus().len(), 2);
        let mut supports: Vec<usize> = result
            .endpoint_stats()
            .iter()
            .map(|stats| stats.support)
            .collect();
        supports.sort_unstable();
        assert_eq!(supports, vec![2, 2]);
        assert_ne!(result.read_to_tu()[0], result.read_to_tu()[2]);
    }

    #[test]
    fn modal_boundary_pair_seeds_the_family_anchor() {
        // A jitter cloud chained by score1: the family anchor must sit on the
        // modal boundary pair (three reads at [4, 100)) so both edge reads
        // qualify directly, instead of anchoring on the 5'-most edge read and
        // splitting the far edge into a second TU.
        let reads = vec![
            read("chr1", Strand::Plus, 0, 100, "edge_low"),
            read("chr1", Strand::Plus, 4, 100, "mode1"),
            read("chr1", Strand::Plus, 4, 100, "mode2"),
            read("chr1", Strand::Plus, 4, 100, "mode3"),
            read("chr1", Strand::Plus, 8, 100, "edge_high"),
        ];

        let result = cluster_tus_with_options(
            &reads,
            0.95,
            0.80,
            TuClusteringOptions {
                attach_contained_reads: false,
                ..TuClusteringOptions::default()
            },
        )
        .unwrap();

        assert_eq!(result.tus().len(), 1);
        assert_eq!(
            result.tus()[0].interval,
            Interval::new(Coord::new(4), Coord::new(100)).unwrap()
        );
        assert_eq!(result.endpoint_stats()[0].support, reads.len());
    }

    #[test]
    fn contained_fragments_assign_partial_to_the_covering_tu() {
        // Degradation fragments of a TU count toward it: both the 3'-anchored
        // fragment and the fully internal fragment are partial evidence for
        // TU_A, while a read overhanging the TU 3' end beyond the jitter
        // window is not explained by it and stays unassigned.
        let tus = vec![tu("chr", Strand::Plus, 0, 3000, "TU_A")];
        let reads = vec![
            read("chr", Strand::Plus, 1800, 3005, "fragment"),
            read("chr", Strand::Plus, 1800, 2500, "internal"),
            read("chr", Strand::Plus, 1800, 3400, "overhang"),
        ];

        let assignments = assign_reads_to_tus(
            &reads,
            &tus,
            0.95,
            0.80,
            TuClusteringOptions::default(),
            0.0,
            false,
        )
        .unwrap();

        assert_eq!(assignments[0].status, AssignmentStatus::Partial);
        assert_eq!(assignments[0].best.unwrap().tu_index, 0);
        assert_eq!(assignments[1].status, AssignmentStatus::Partial);
        assert_eq!(assignments[1].best.unwrap().tu_index, 0);
        assert_eq!(assignments[2].status, AssignmentStatus::Unassigned);
    }
}
