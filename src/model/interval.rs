use thiserror::Error;

use super::Coord;

#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash, PartialOrd, Ord)]
pub struct Interval {
    pub start: Coord,
    pub end: Coord,
}

#[derive(Error, Debug)]
pub enum IntervalError {
    #[error("invalid interval: start {start} > end {end}")]
    StartAfterEnd { start: Coord, end: Coord },
}

impl Interval {
    pub fn new(start: Coord, end: Coord) -> Result<Self, IntervalError> {
        if start > end {
            return Err(IntervalError::StartAfterEnd { start, end });
        }
        Ok(Self { start, end })
    }

    pub fn len(self) -> u32 {
        self.end.get() - self.start.get()
    }

    pub fn is_empty(self) -> bool {
        self.start == self.end
    }

    pub fn overlap_len(self, other: Self) -> u32 {
        let start = std::cmp::max(self.start, other.start);
        let end = std::cmp::min(self.end, other.end);
        if start < end {
            end.get() - start.get()
        } else {
            0
        }
    }

    pub fn intersection(self, other: Self) -> Option<Self> {
        let start = std::cmp::max(self.start, other.start);
        let end = std::cmp::min(self.end, other.end);
        if start < end {
            Some(Self { start, end })
        } else {
            None
        }
    }

    pub fn overlaps(self, other: Self) -> bool {
        self.overlap_len(other) > 0
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use proptest::prelude::*;

    #[test]
    fn overlaps_half_open_truth_table() {
        let a = Interval::new(Coord::new(0), Coord::new(1)).unwrap();
        let b = Interval::new(Coord::new(1), Coord::new(2)).unwrap();
        assert!(!a.overlaps(b));
        assert!(!b.overlaps(a));

        let a = Interval::new(Coord::new(0), Coord::new(2)).unwrap();
        let b = Interval::new(Coord::new(1), Coord::new(2)).unwrap();
        assert!(a.overlaps(b));
        assert!(b.overlaps(a));

        let a = Interval::new(Coord::new(5), Coord::new(5)).unwrap();
        let b = Interval::new(Coord::new(5), Coord::new(6)).unwrap();
        assert!(!a.overlaps(b));
        assert!(!b.overlaps(a));

        let a = Interval::new(Coord::new(5), Coord::new(5)).unwrap();
        let b = Interval::new(Coord::new(4), Coord::new(6)).unwrap();
        assert!(!a.overlaps(b));
        assert!(!b.overlaps(a));
    }

    proptest! {
        #[test]
        fn overlap_is_symmetric(a0 in 0u32..1000, a1 in 0u32..1000, b0 in 0u32..1000, b1 in 0u32..1000) {
            let (a_start, a_end) = if a0 <= a1 { (a0, a1) } else { (a1, a0) };
            let (b_start, b_end) = if b0 <= b1 { (b0, b1) } else { (b1, b0) };

            let a = Interval::new(Coord::new(a_start), Coord::new(a_end)).unwrap();
            let b = Interval::new(Coord::new(b_start), Coord::new(b_end)).unwrap();

            prop_assert_eq!(a.overlaps(b), b.overlaps(a));
            prop_assert_eq!(a.overlap_len(b) > 0, a.overlaps(b));
        }

        #[test]
        fn empty_intervals_never_overlap(x in 0u32..1000, b0 in 0u32..1000, b1 in 0u32..1000) {
            let (b_start, b_end) = if b0 <= b1 { (b0, b1) } else { (b1, b0) };
            let empty = Interval::new(Coord::new(x), Coord::new(x)).unwrap();
            let b = Interval::new(Coord::new(b_start), Coord::new(b_end)).unwrap();

            prop_assert_eq!(empty.overlap_len(b), 0);
            prop_assert!(!empty.overlaps(b));
            prop_assert!(!b.overlaps(empty));
        }
    }
}
