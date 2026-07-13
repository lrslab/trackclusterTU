//! Symmetric overlap scores for single spans and exon interval lists.

use crate::interval::{intersection_len, total_len, union_len};
use crate::model::Interval;

/// Compute span-Jaccard similarity: intersection length divided by union length.
///
/// This legacy `score1` API is symmetric and returns a value in `[0, 1]`. Two
/// empty intervals return `0.0` because their union has zero length.
///
/// # Examples
///
/// ```
/// use trackclustertu::model::{Coord, Interval};
/// use trackclustertu::score::score1_interval;
///
/// let left = Interval::new(Coord::new(0), Coord::new(10))?;
/// let right = Interval::new(Coord::new(5), Coord::new(15))?;
/// assert!((score1_interval(left, right) - (5.0 / 15.0)).abs() < 1e-12);
/// assert_eq!(score1_interval(left, right), score1_interval(right, left));
/// # Ok::<(), trackclustertu::model::IntervalError>(())
/// ```
pub fn score1_interval(a: Interval, b: Interval) -> f64 {
    let overlap = a.overlap_len(b) as u64;
    let union = a.len() as u64 + b.len() as u64 - overlap;
    if union == 0 {
        0.0
    } else {
        overlap as f64 / union as f64
    }
}

/// Compute overlap-over-longer similarity for two spans.
///
/// This legacy `score2` API divides the shared bases by the longer span length,
/// penalizing short/long containment. It is symmetric, ranges from `0.0` to
/// `1.0`, and returns `0.0` when both spans are empty.
pub fn score2_interval(a: Interval, b: Interval) -> f64 {
    let overlap = a.overlap_len(b) as u64;
    let max_len = (a.len() as u64).max(b.len() as u64);
    if max_len == 0 {
        0.0
    } else {
        overlap as f64 / max_len as f64
    }
}

/// Compute span-Jaccard similarity across sorted, non-overlapping interval lists.
///
/// Each input must be ordered by start coordinate and internally non-overlapping;
/// debug builds assert this invariant. The score is `0.0` when the combined
/// union is empty.
///
/// # Examples
///
/// ```
/// use trackclustertu::model::{Coord, Interval};
/// use trackclustertu::score::{score1_intervals, score2_intervals};
///
/// let interval = |start, end| Interval::new(Coord::new(start), Coord::new(end));
/// let left = [interval(0, 10)?, interval(20, 30)?];
/// let right = [interval(5, 10)?, interval(20, 25)?];
/// assert!((score1_intervals(&left, &right) - (10.0 / 20.0)).abs() < 1e-12);
/// assert!((score2_intervals(&left, &right) - (10.0 / 20.0)).abs() < 1e-12);
/// # Ok::<(), trackclustertu::model::IntervalError>(())
/// ```
pub fn score1_intervals(a: &[Interval], b: &[Interval]) -> f64 {
    let overlap = intersection_len(a, b);
    let union = union_len(a, b);
    if union == 0 {
        0.0
    } else {
        overlap as f64 / union as f64
    }
}

/// Compute overlap-over-longer similarity across interval lists.
///
/// Each input must be ordered by start coordinate and internally non-overlapping;
/// debug builds assert this invariant. The score is `0.0` when both lists have
/// zero total length.
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
