use std::cmp::Reverse;
use std::collections::{BinaryHeap, HashSet};

use crate::model::{Interval, Transcript};

use super::{partition, StrandMode};

#[derive(Clone, Copy, Debug, Default)]
pub struct IntersectOpts {
    pub strand_mode: StrandMode,
    pub min_overlap_bp: Option<u32>,
}

fn span(transcript: &Transcript) -> Interval {
    Interval {
        start: transcript.tx_start,
        end: transcript.tx_end,
    }
}

fn span_overlap_len(a: &Transcript, b: &Transcript) -> u32 {
    span(a).overlap_len(span(b))
}

fn add_active(
    transcript_end: u32,
    idx: usize,
    active: &mut HashSet<usize>,
    ends: &mut BinaryHeap<Reverse<(u32, usize)>>,
) {
    active.insert(idx);
    ends.push(Reverse((transcript_end, idx)));
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

fn sweep_partition_pairs(
    a: &[Transcript],
    a_indices: &[usize],
    b: &[Transcript],
    b_indices: &[usize],
    opts: &IntersectOpts,
    out: &mut Vec<(usize, usize)>,
) {
    let mut ai: usize = 0;
    let mut bi: usize = 0;

    let mut active_a: HashSet<usize> = HashSet::new();
    let mut active_b: HashSet<usize> = HashSet::new();

    let mut a_ends: BinaryHeap<Reverse<(u32, usize)>> = BinaryHeap::new();
    let mut b_ends: BinaryHeap<Reverse<(u32, usize)>> = BinaryHeap::new();

    loop {
        let next_a_start = a_indices
            .get(ai)
            .map(|&idx| a[idx].tx_start.get())
            .unwrap_or(u32::MAX);
        let next_b_start = b_indices
            .get(bi)
            .map(|&idx| b[idx].tx_start.get())
            .unwrap_or(u32::MAX);

        if next_a_start == u32::MAX && next_b_start == u32::MAX {
            break;
        }

        let current_start = next_a_start.min(next_b_start);
        expire_active(current_start, &mut active_a, &mut a_ends);
        expire_active(current_start, &mut active_b, &mut b_ends);

        if next_a_start <= next_b_start {
            let a_idx = a_indices[ai];
            ai += 1;
            let a_tx = &a[a_idx];
            if a_tx.tx_start == a_tx.tx_end {
                continue;
            }

            add_active(a_tx.tx_end.get(), a_idx, &mut active_a, &mut a_ends);
            for &b_idx in &active_b {
                let b_tx = &b[b_idx];

                if opts.strand_mode == StrandMode::Match && a_tx.strand != b_tx.strand {
                    continue;
                }

                let overlap = span_overlap_len(a_tx, b_tx);
                if let Some(min_overlap) = opts.min_overlap_bp {
                    if overlap < min_overlap {
                        continue;
                    }
                }
                out.push((a_idx, b_idx));
            }
        } else {
            let b_idx = b_indices[bi];
            bi += 1;
            let b_tx = &b[b_idx];
            if b_tx.tx_start == b_tx.tx_end {
                continue;
            }

            add_active(b_tx.tx_end.get(), b_idx, &mut active_b, &mut b_ends);
            for &a_idx in &active_a {
                let a_tx = &a[a_idx];

                if opts.strand_mode == StrandMode::Match && a_tx.strand != b_tx.strand {
                    continue;
                }

                let overlap = span_overlap_len(a_tx, b_tx);
                if let Some(min_overlap) = opts.min_overlap_bp {
                    if overlap < min_overlap {
                        continue;
                    }
                }
                out.push((a_idx, b_idx));
            }
        }
    }
}

/// Returns all `(a_index, b_index)` pairs whose transcript spans overlap (half-open).
///
/// Requirements:
/// - Inputs must be sorted by `(chrom, start, end, strand)` for linear-time per partition.
pub fn sweep_intersect_pairs(
    a: &[Transcript],
    b: &[Transcript],
    opts: &IntersectOpts,
) -> Vec<(usize, usize)> {
    let a_parts = partition(a, opts.strand_mode);
    let b_parts = partition(b, opts.strand_mode);

    let mut pairs: Vec<(usize, usize)> = Vec::new();
    for (key, a_indices) in a_parts {
        let Some(b_indices) = b_parts.get(&key) else {
            continue;
        };
        sweep_partition_pairs(a, &a_indices, b, b_indices, opts, &mut pairs);
    }

    pairs.sort_unstable();
    pairs
}

#[cfg(test)]
mod tests {
    use crate::interval::sort::sort_by_coord;
    use crate::model::{Bed12Attrs, Coord, Interval, Strand, Transcript};
    use proptest::prelude::*;

    use super::*;

    fn make_tx(chrom: &str, strand: Strand, start: u32, end: u32, name: &str) -> Transcript {
        Transcript::new(
            chrom.to_owned(),
            strand,
            Coord::new(start),
            Coord::new(end),
            name.to_owned(),
            vec![Interval::new(Coord::new(start), Coord::new(end)).unwrap()],
            Bed12Attrs {
                score: 0,
                thick_start: Coord::new(start),
                thick_end: Coord::new(end),
                item_rgb: "0".to_owned(),
                extra_fields: Vec::new(),
            },
        )
        .unwrap()
    }

    fn naive_pairs(
        a: &[Transcript],
        b: &[Transcript],
        opts: &IntersectOpts,
    ) -> Vec<(usize, usize)> {
        let mut out = Vec::new();
        for (ai, atx) in a.iter().enumerate() {
            for (bi, btx) in b.iter().enumerate() {
                if atx.chrom != btx.chrom {
                    continue;
                }
                if opts.strand_mode == StrandMode::Match && atx.strand != btx.strand {
                    continue;
                }
                let overlap = span_overlap_len(atx, btx);
                if overlap == 0 {
                    continue;
                }
                if let Some(min_overlap) = opts.min_overlap_bp {
                    if overlap < min_overlap {
                        continue;
                    }
                }
                out.push((ai, bi));
            }
        }
        out.sort_unstable();
        out
    }

    #[test]
    fn sweep_matches_naive_simple() {
        let mut a = vec![
            make_tx("chr1", Strand::Plus, 10, 20, "a1"),
            make_tx("chr1", Strand::Plus, 30, 40, "a2"),
        ];
        let mut b = vec![
            make_tx("chr1", Strand::Plus, 15, 16, "b1"),
            make_tx("chr1", Strand::Minus, 18, 19, "b2"),
            make_tx("chr1", Strand::Plus, 35, 45, "b3"),
        ];
        sort_by_coord(&mut a);
        sort_by_coord(&mut b);

        let opts = IntersectOpts {
            strand_mode: StrandMode::Ignore,
            min_overlap_bp: None,
        };
        assert_eq!(
            sweep_intersect_pairs(&a, &b, &opts),
            naive_pairs(&a, &b, &opts)
        );

        let opts = IntersectOpts {
            strand_mode: StrandMode::Match,
            min_overlap_bp: None,
        };
        assert_eq!(
            sweep_intersect_pairs(&a, &b, &opts),
            naive_pairs(&a, &b, &opts)
        );
    }

    #[test]
    fn sweep_respects_min_overlap() {
        let mut a = vec![make_tx("chr1", Strand::Plus, 10, 20, "a")];
        let mut b = vec![make_tx("chr1", Strand::Plus, 19, 25, "b")];
        sort_by_coord(&mut a);
        sort_by_coord(&mut b);

        let opts = IntersectOpts {
            strand_mode: StrandMode::Ignore,
            min_overlap_bp: Some(2),
        };
        assert_eq!(
            sweep_intersect_pairs(&a, &b, &opts),
            Vec::<(usize, usize)>::new()
        );
    }

    fn make_tx_generated(
        chrom: String,
        strand: Strand,
        start: u32,
        len: u32,
        name: String,
    ) -> Transcript {
        let end = start.saturating_add(len);
        Transcript::new(
            chrom,
            strand,
            Coord::new(start),
            Coord::new(end),
            name,
            vec![Interval::new(Coord::new(start), Coord::new(end)).unwrap()],
            Bed12Attrs {
                score: 0,
                thick_start: Coord::new(start),
                thick_end: Coord::new(end),
                item_rgb: "0".to_owned(),
                extra_fields: Vec::new(),
            },
        )
        .unwrap()
    }

    proptest! {
        #[test]
        fn sweep_matches_naive_random(
            a_specs in prop::collection::vec((prop_oneof![Just("chr1".to_owned()), Just("chr2".to_owned())],
                                              prop_oneof![Just(Strand::Plus), Just(Strand::Minus), Just(Strand::Unknown)],
                                              0u32..200,
                                              0u32..50), 0..12),
            b_specs in prop::collection::vec((prop_oneof![Just("chr1".to_owned()), Just("chr2".to_owned())],
                                              prop_oneof![Just(Strand::Plus), Just(Strand::Minus), Just(Strand::Unknown)],
                                              0u32..200,
                                              0u32..50), 0..12),
            strand_mode in prop_oneof![Just(StrandMode::Ignore), Just(StrandMode::Match)],
            min_overlap_bp in prop_oneof![Just(None), (1u32..25).prop_map(Some)],
        ) {
            let mut a: Vec<Transcript> = a_specs.into_iter().enumerate().map(|(i, (chrom, strand, start, len))| {
                make_tx_generated(chrom, strand, start, len, format!("a{i}"))
            }).collect();
            let mut b: Vec<Transcript> = b_specs.into_iter().enumerate().map(|(i, (chrom, strand, start, len))| {
                make_tx_generated(chrom, strand, start, len, format!("b{i}"))
            }).collect();

            sort_by_coord(&mut a);
            sort_by_coord(&mut b);

            let opts = IntersectOpts {
                strand_mode,
                min_overlap_bp,
            };

            let sweep = sweep_intersect_pairs(&a, &b, &opts);
            let naive = naive_pairs(&a, &b, &opts);
            prop_assert_eq!(sweep, naive);
        }
    }
}
