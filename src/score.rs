use crate::interval::{intersection_len, total_len, union_len};
use crate::model::Interval;

pub fn score1_interval(a: Interval, b: Interval) -> f64 {
    let overlap = a.overlap_len(b) as u64;
    let union = a.len() as u64 + b.len() as u64 - overlap;
    if union == 0 {
        0.0
    } else {
        overlap as f64 / union as f64
    }
}

pub fn score2_interval(a: Interval, b: Interval) -> f64 {
    let overlap = a.overlap_len(b) as u64;
    let max_len = (a.len() as u64).max(b.len() as u64);
    if max_len == 0 {
        0.0
    } else {
        overlap as f64 / max_len as f64
    }
}

pub fn score1_intervals(a: &[Interval], b: &[Interval]) -> f64 {
    let overlap = intersection_len(a, b);
    let union = union_len(a, b);
    if union == 0 {
        0.0
    } else {
        overlap as f64 / union as f64
    }
}

pub fn score2_intervals(a: &[Interval], b: &[Interval]) -> f64 {
    let overlap = intersection_len(a, b);
    let max_len = total_len(a).max(total_len(b));
    if max_len == 0 {
        0.0
    } else {
        overlap as f64 / max_len as f64
    }
}

#[cfg(test)]
mod tests {
    use crate::model::Coord;

    use super::*;

    fn interval(start: u32, end: u32) -> Interval {
        Interval::new(Coord::new(start), Coord::new(end)).unwrap()
    }

    #[test]
    fn score1_interval_is_one_for_identical() {
        let a = interval(0, 10);
        assert_eq!(score1_interval(a, a), 1.0);
    }

    #[test]
    fn score1_interval_is_zero_for_disjoint() {
        let a = interval(0, 10);
        let b = interval(10, 20);
        assert_eq!(score1_interval(a, b), 0.0);
        assert_eq!(score2_interval(a, b), 0.0);
    }

    #[test]
    fn score2_interval_penalizes_containment_by_length_ratio() {
        let a = interval(0, 10);
        let b = interval(2, 8);
        assert_eq!(score2_interval(a, b), 0.6);
        assert_eq!(score2_interval(b, a), 0.6);
        assert_eq!(score1_interval(a, b), 0.6);
    }

    #[test]
    fn scores_return_zero_when_denominator_is_zero() {
        let empty = interval(5, 5);
        assert_eq!(score1_interval(empty, empty), 0.0);
        assert_eq!(score2_interval(empty, empty), 0.0);
    }
}
