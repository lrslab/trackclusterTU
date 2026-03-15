pub mod multi;

use std::cmp::{Ordering, Reverse};
use std::collections::{BinaryHeap, HashMap, HashSet};

use rayon::prelude::*;
use thiserror::Error;

use crate::model::{Interval, Strand};
use crate::score::{score1_interval, score2_interval};

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct ReadRecord {
    pub contig: String,
    pub strand: Strand,
    pub interval: Interval,
    pub id: String,
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Tu {
    pub id: String,
    pub contig: String,
    pub strand: Strand,
    pub interval: Interval,
    pub rep_read_index: usize,
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct TuClusteringResult {
    pub tus: Vec<Tu>,
    pub read_to_tu: Vec<usize>,
}

#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct TuClusteringStats {
    pub partition_count: usize,
    pub max_partition_reads: usize,
    pub region_count: usize,
    pub max_region_reads: usize,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct TuClusteringOptions {
    pub attach_contained_reads: bool,
}

impl Default for TuClusteringOptions {
    fn default() -> Self {
        Self {
            attach_contained_reads: true,
        }
    }
}

#[derive(Error, Debug)]
pub enum TuError {
    #[error("score1_threshold must be between 0 and 1, got {value}")]
    InvalidScore1Threshold { value: f64 },

    #[error("score2_threshold must be between 0 and 1, got {value}")]
    InvalidScore2Threshold { value: f64 },
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
        .then_with(|| a.interval.start.cmp(&b.interval.start))
        .then_with(|| a.interval.end.cmp(&b.interval.end))
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

        let current_start = read.interval.start.get();
        let current_end = read.interval.end.get();
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
            if other.interval.end.get() <= current_start {
                continue;
            }

            let other_len_u32 = read_len(other);
            if other_len_u32 == 0 {
                continue;
            }

            let start_delta = current_start.saturating_sub(other.interval.start.get());
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

#[derive(Clone, Debug)]
struct TuNoId {
    contig: String,
    strand: Strand,
    interval: Interval,
    rep_read_index: usize,
}

#[derive(Clone, Debug)]
struct PartitionClusteringResult {
    tus: Vec<TuNoId>,
    local_read_to_tu: Vec<usize>,
}

fn partition_into_regions(indices_sorted: &[usize], reads: &[ReadRecord]) -> Vec<(usize, usize)> {
    if indices_sorted.is_empty() {
        return Vec::new();
    }

    let mut regions: Vec<(usize, usize)> = Vec::new();
    let mut start = 0usize;
    let mut current_end = reads[indices_sorted[0]].interval.end.get();

    for pos in 1..indices_sorted.len() {
        let read = &reads[indices_sorted[pos]];
        let read_start = read.interval.start.get();
        if read_start >= current_end {
            regions.push((start, pos));
            start = pos;
            current_end = read.interval.end.get();
        } else {
            current_end = current_end.max(read.interval.end.get());
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

    let mut clusters: Vec<SeedCluster> = Vec::with_capacity(clusters_by_root.len());
    for (_, mut members) in clusters_by_root {
        members.sort_by(|&a, &b| cmp_read_idx(reads, indices_sorted[a], indices_sorted[b]));

        let mut rep_pos = members[0];
        for &candidate_pos in &members[1..] {
            let candidate = &reads[indices_sorted[candidate_pos]];
            let current = &reads[indices_sorted[rep_pos]];
            if better_representative(candidate, current) {
                rep_pos = candidate_pos;
            }
        }

        clusters.push(SeedCluster {
            rep_read_index: indices_sorted[rep_pos],
            members,
        });
    }

    clusters.sort_by(|a, b| cmp_read_idx(reads, a.rep_read_index, b.rep_read_index));

    let mut parents: Vec<Option<usize>> = vec![None; clusters.len()];
    if options.attach_contained_reads {
        let mut active: HashSet<usize> = HashSet::new();
        let mut ends: BinaryHeap<Reverse<(u32, usize)>> = BinaryHeap::new();

        for cluster_idx in 0..clusters.len() {
            let rep_idx = clusters[cluster_idx].rep_read_index;
            let child = &reads[rep_idx];
            if child.interval.is_empty() {
                continue;
            }

            expire_active(child.interval.start.get(), &mut active, &mut ends);

            let mut best: Option<(usize, f64, f64, &str)> = None;
            for &cand_cluster_idx in &active {
                let cand_rep_idx = clusters[cand_cluster_idx].rep_read_index;
                let candidate = &reads[cand_rep_idx];

                let cand_len = read_len(candidate) as u64;
                let child_len = read_len(child) as u64;
                if cand_len <= child_len || child_len == 0 {
                    continue;
                }

                if candidate.interval.end < child.interval.end {
                    continue;
                }

                let s2 = score2_interval(child.interval, candidate.interval);
                if s2 < score2_threshold {
                    continue;
                }
                let s1 = score1_interval(child.interval, candidate.interval);

                match best {
                    None => best = Some((cand_cluster_idx, s2, s1, candidate.id.as_str())),
                    Some((best_idx, best_s2, best_s1, best_id)) => {
                        let better = match s2.total_cmp(&best_s2) {
                            Ordering::Greater => true,
                            Ordering::Less => false,
                            Ordering::Equal => match s1.total_cmp(&best_s1) {
                                Ordering::Greater => true,
                                Ordering::Less => false,
                                Ordering::Equal => {
                                    candidate.id.as_str() < best_id
                                        || (candidate.id.as_str() == best_id
                                            && cand_cluster_idx < best_idx)
                                }
                            },
                        };
                        if better {
                            best = Some((cand_cluster_idx, s2, s1, candidate.id.as_str()));
                        }
                    }
                }
            }

            if let Some((parent_idx, _, _, _)) = best {
                parents[cluster_idx] = Some(parent_idx);
            }

            add_active(
                child.interval.end.get(),
                cluster_idx,
                &mut active,
                &mut ends,
            );
        }
    }

    let mut cluster_to_root: Vec<usize> = Vec::with_capacity(clusters.len());
    for idx in 0..clusters.len() {
        cluster_to_root.push(resolve_root(idx, &mut parents));
    }

    let mut root_clusters: Vec<usize> = Vec::new();
    for (idx, &root) in cluster_to_root.iter().enumerate() {
        if idx == root {
            root_clusters.push(idx);
        }
    }

    let mut tus: Vec<TuNoId> = Vec::with_capacity(root_clusters.len());
    let mut root_to_tu: HashMap<usize, usize> = HashMap::with_capacity(root_clusters.len());

    for (tu_index, &cluster_idx) in root_clusters.iter().enumerate() {
        let rep_read_index = clusters[cluster_idx].rep_read_index;
        let rep = &reads[rep_read_index];
        tus.push(TuNoId {
            contig: rep.contig.clone(),
            strand: rep.strand,
            interval: rep.interval,
            rep_read_index,
        });
        root_to_tu.insert(cluster_idx, tu_index);
    }

    let mut local_read_to_tu: Vec<usize> = vec![usize::MAX; indices_sorted.len()];
    for (cluster_idx, cluster) in clusters.iter().enumerate() {
        let root = cluster_to_root[cluster_idx];
        let tu_index = *root_to_tu
            .get(&root)
            .expect("root cluster must map to TU index");
        for &local_read_pos in &cluster.members {
            local_read_to_tu[local_read_pos] = tu_index;
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

    if reads.is_empty() {
        return Ok((
            TuClusteringResult {
                tus: Vec::new(),
                read_to_tu: Vec::new(),
            },
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
    let mut read_to_tu: Vec<usize> = vec![usize::MAX; reads.len()];
    let mut tu_offset: usize = 0;

    for (partition_idx, (result, _)) in partition_results.into_iter().enumerate() {
        let PartitionClusteringResult {
            tus: tus_no_id,
            local_read_to_tu,
        } = result;

        let partition_tu_count = tus_no_id.len();
        for (local_tu_idx, tu) in tus_no_id.into_iter().enumerate() {
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
        TuClusteringResult { tus, read_to_tu },
        TuClusteringStats {
            partition_count,
            max_partition_reads,
            region_count,
            max_region_reads,
        },
    ))
}

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

pub fn cluster_tus_with_options(
    reads: &[ReadRecord],
    score1_threshold: f64,
    score2_threshold: f64,
    options: TuClusteringOptions,
) -> Result<TuClusteringResult, TuError> {
    Ok(cluster_tus_with_stats_options(reads, score1_threshold, score2_threshold, options)?.0)
}

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
    use std::collections::HashMap;

    use crate::model::Coord;

    use super::*;

    fn read(contig: &str, strand: Strand, start: u32, end: u32, id: &str) -> ReadRecord {
        ReadRecord {
            contig: contig.to_owned(),
            strand,
            interval: Interval::new(Coord::new(start), Coord::new(end)).unwrap(),
            id: id.to_owned(),
        }
    }

    #[test]
    fn clusters_and_attaches_truncations() {
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
        assert_eq!(result.tus.len(), 3);

        assert_eq!(
            result.tus[0],
            Tu {
                id: "TU000001".to_owned(),
                contig: "chr1".to_owned(),
                strand: Strand::Plus,
                interval: Interval::new(Coord::new(100), Coord::new(200)).unwrap(),
                rep_read_index: 1,
            }
        );
        assert_eq!(
            result.tus[1],
            Tu {
                id: "TU000002".to_owned(),
                contig: "chr1".to_owned(),
                strand: Strand::Plus,
                interval: Interval::new(Coord::new(300), Coord::new(400)).unwrap(),
                rep_read_index: 3,
            }
        );
        assert_eq!(
            result.tus[2],
            Tu {
                id: "TU000003".to_owned(),
                contig: "chr1".to_owned(),
                strand: Strand::Minus,
                interval: Interval::new(Coord::new(100), Coord::new(200)).unwrap(),
                rep_read_index: 6,
            }
        );

        let mut by_id: HashMap<String, String> = HashMap::new();
        for (read_idx, read) in reads.iter().enumerate() {
            let tu_id = &result.tus[result.read_to_tu[read_idx]].id;
            by_id.insert(read.id.clone(), tu_id.clone());
        }

        assert_eq!(by_id["r1"], "TU000001");
        assert_eq!(by_id["r2"], "TU000001");
        assert_eq!(by_id["r3"], "TU000001");
        assert_eq!(by_id["r4"], "TU000002");
        assert_eq!(by_id["r5"], "TU000002");
        assert_eq!(by_id["r6"], "TU000002");
        assert_eq!(by_id["r7"], "TU000003");
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
            },
        )
        .unwrap();
        assert_eq!(result.tus.len(), 5);

        assert_eq!(
            result.tus[0],
            Tu {
                id: "TU000001".to_owned(),
                contig: "chr1".to_owned(),
                strand: Strand::Plus,
                interval: Interval::new(Coord::new(100), Coord::new(200)).unwrap(),
                rep_read_index: 1,
            }
        );
        assert_eq!(
            result.tus[1],
            Tu {
                id: "TU000002".to_owned(),
                contig: "chr1".to_owned(),
                strand: Strand::Plus,
                interval: Interval::new(Coord::new(120), Coord::new(180)).unwrap(),
                rep_read_index: 0,
            }
        );
        assert_eq!(
            result.tus[2],
            Tu {
                id: "TU000003".to_owned(),
                contig: "chr1".to_owned(),
                strand: Strand::Plus,
                interval: Interval::new(Coord::new(300), Coord::new(400)).unwrap(),
                rep_read_index: 3,
            }
        );
        assert_eq!(
            result.tus[3],
            Tu {
                id: "TU000004".to_owned(),
                contig: "chr1".to_owned(),
                strand: Strand::Plus,
                interval: Interval::new(Coord::new(320), Coord::new(360)).unwrap(),
                rep_read_index: 5,
            }
        );
        assert_eq!(
            result.tus[4],
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
            let tu_id = &result.tus[result.read_to_tu[read_idx]].id;
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
}
