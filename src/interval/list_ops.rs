use crate::model::Interval;

fn sorted_non_overlapping(xs: &[Interval]) -> bool {
    xs.windows(2).all(|window| {
        let left = window[0];
        let right = window[1];
        left.start <= right.start && left.end <= right.start
    })
}

pub fn total_len(xs: &[Interval]) -> u64 {
    debug_assert!(sorted_non_overlapping(xs));
    xs.iter().map(|interval| interval.len() as u64).sum()
}

pub fn intersection_len(a: &[Interval], b: &[Interval]) -> u64 {
    debug_assert!(sorted_non_overlapping(a));
    debug_assert!(sorted_non_overlapping(b));

    let mut total: u64 = 0;
    let mut ai: usize = 0;
    let mut bi: usize = 0;

    while ai < a.len() && bi < b.len() {
        let a_interval = a[ai];
        let b_interval = b[bi];

        if a_interval.end <= b_interval.start {
            ai += 1;
            continue;
        }
        if b_interval.end <= a_interval.start {
            bi += 1;
            continue;
        }

        total += a_interval.overlap_len(b_interval) as u64;

        if a_interval.end <= b_interval.end {
            ai += 1;
        } else {
            bi += 1;
        }
    }

    total
}

pub fn union_len(a: &[Interval], b: &[Interval]) -> u64 {
    total_len(a) + total_len(b) - intersection_len(a, b)
}

pub fn merge_overlaps(mut xs: Vec<Interval>) -> Vec<Interval> {
    xs.sort_unstable_by_key(|interval| (interval.start, interval.end));
    let mut merged: Vec<Interval> = Vec::new();
    for interval in xs {
        match merged.last_mut() {
            Some(last) if interval.start < last.end => {
                if interval.end > last.end {
                    last.end = interval.end;
                }
            }
            _ => merged.push(interval),
        }
    }
    merged
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::model::Coord;
    use proptest::prelude::*;

    fn interval(start: u32, end: u32) -> Interval {
        Interval::new(Coord::new(start), Coord::new(end)).unwrap()
    }

    fn coverage(xs: &[Interval], max: usize) -> Vec<bool> {
        let mut cov = vec![false; max];
        for interval in xs {
            let start = interval.start.get() as usize;
            let end = interval.end.get() as usize;
            let start = start.min(max);
            let end = end.min(max);
            cov[start..end].fill(true);
        }
        cov
    }

    #[test]
    fn intersection_len_is_zero_for_touching_intervals() {
        let a = [interval(0, 10)];
        let b = [interval(10, 20)];
        assert_eq!(intersection_len(&a, &b), 0);
        assert_eq!(union_len(&a, &b), 20);
    }

    #[test]
    fn union_len_matches_expected_on_simple_overlap() {
        let a = [interval(0, 10)];
        let b = [interval(5, 15)];
        assert_eq!(intersection_len(&a, &b), 5);
        assert_eq!(union_len(&a, &b), 15);
    }

    proptest! {
        #[test]
        fn intersection_and_union_match_naive(
            a in prop::collection::vec((0u32..100, 0u32..100), 0..20),
            b in prop::collection::vec((0u32..100, 0u32..100), 0..20),
        ) {
            let mut a_intervals: Vec<Interval> = a.into_iter()
                .map(|(s, e)| if s <= e { interval(s, e) } else { interval(e, s) })
                .collect();
            let mut b_intervals: Vec<Interval> = b.into_iter()
                .map(|(s, e)| if s <= e { interval(s, e) } else { interval(e, s) })
                .collect();

            a_intervals = merge_overlaps(a_intervals);
            b_intervals = merge_overlaps(b_intervals);

            prop_assert!(sorted_non_overlapping(&a_intervals));
            prop_assert!(sorted_non_overlapping(&b_intervals));

            let a_cov = coverage(&a_intervals, 100);
            let b_cov = coverage(&b_intervals, 100);
            let naive_intersection: u64 = a_cov.iter().zip(b_cov.iter()).filter(|(a,b)| **a && **b).count() as u64;
            let naive_union: u64 = a_cov.iter().zip(b_cov.iter()).filter(|(a,b)| **a || **b).count() as u64;

            prop_assert_eq!(intersection_len(&a_intervals, &b_intervals), naive_intersection);
            prop_assert_eq!(union_len(&a_intervals, &b_intervals), naive_union);
        }
    }
}
