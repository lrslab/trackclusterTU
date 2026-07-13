//! Deterministic endpoint-mode detection for boundary-aware transcript calling.
//!
//! The implementation repeatedly selects the supported observed coordinate with
//! the most remaining observations in a centered window, emits that mode, and
//! removes its window. A Fenwick tree supplies indexed window sums, giving a
//! worst-case complexity of O(k² log k) for k distinct coordinates rather than
//! cubic repeated rescans.
//!
//! ```
//! use trackclustertu::endpoint::{detect_endpoint_modes, TieDirection};
//!
//! let modes = detect_endpoint_modes([100, 100, 103, 150], 5, TieDirection::PreferLower);
//! assert_eq!(modes[0].coordinate(), 100);
//! assert_eq!(modes[0].support(), 3);
//! ```

use std::collections::BTreeMap;

/// Direction used to resolve endpoint modes with otherwise equal support.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum TieDirection {
    /// Prefer the lower genomic coordinate.
    PreferLower,
    /// Prefer the higher genomic coordinate.
    PreferHigher,
}

/// One non-overlapping endpoint mode returned by [`detect_endpoint_modes`].
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct EndpointMode {
    coordinate: u32,
    support: usize,
}

impl EndpointMode {
    /// The observed genomic coordinate selected as the mode center.
    pub fn coordinate(self) -> u32 {
        self.coordinate
    }

    /// Number of observations assigned to this centered window.
    pub fn support(self) -> usize {
        self.support
    }
}

/// Detect deterministic, non-overlapping modes in genomic endpoint coordinates.
///
/// Each emitted center is one of the observed coordinates. Support includes
/// remaining coordinates within `window_bp` on either side. After selection,
/// all observations in that window are removed before finding the next mode.
pub fn detect_endpoint_modes<I>(
    observations: I,
    window_bp: u32,
    tie_direction: TieDirection,
) -> Vec<EndpointMode>
where
    I: IntoIterator<Item = u32>,
{
    let mut counts: BTreeMap<u32, usize> = BTreeMap::new();
    for coordinate in observations {
        *counts.entry(coordinate).or_default() += 1;
    }

    let coordinates: Vec<u32> = counts.keys().copied().collect();
    let coordinate_counts: Vec<usize> = counts.values().copied().collect();
    let mut active = vec![true; coordinates.len()];
    let mut index = FenwickCounts::from_counts(&coordinate_counts);
    let (window_starts, window_ends) = centered_window_bounds(&coordinates, window_bp);
    let mut modes = Vec::new();

    while let Some(best_idx) = (0..coordinates.len())
        .filter(|&candidate_idx| active[candidate_idx])
        .max_by(|&left_idx, &right_idx| {
            let left_support = index.range_sum(window_starts[left_idx], window_ends[left_idx] + 1);
            let right_support =
                index.range_sum(window_starts[right_idx], window_ends[right_idx] + 1);
            left_support
                .cmp(&right_support)
                .then_with(|| coordinate_counts[left_idx].cmp(&coordinate_counts[right_idx]))
                .then_with(|| match tie_direction {
                    TieDirection::PreferLower => coordinates[right_idx].cmp(&coordinates[left_idx]),
                    TieDirection::PreferHigher => {
                        coordinates[left_idx].cmp(&coordinates[right_idx])
                    }
                })
        })
    {
        let start = window_starts[best_idx];
        let end = window_ends[best_idx];
        let support = index.range_sum(start, end + 1);
        modes.push(EndpointMode {
            coordinate: coordinates[best_idx],
            support,
        });
        for coordinate_idx in start..=end {
            if active[coordinate_idx] {
                active[coordinate_idx] = false;
                index.remove(coordinate_idx, coordinate_counts[coordinate_idx]);
            }
        }
    }

    modes
}

fn centered_window_bounds(coordinates: &[u32], window_bp: u32) -> (Vec<usize>, Vec<usize>) {
    let mut starts = vec![0; coordinates.len()];
    let mut ends = vec![0; coordinates.len()];

    let mut left = 0usize;
    for (center_idx, &center) in coordinates.iter().enumerate() {
        while center.abs_diff(coordinates[left]) > window_bp {
            left += 1;
        }
        starts[center_idx] = left;
    }

    let mut right = 0usize;
    for (center_idx, &center) in coordinates.iter().enumerate() {
        if right < center_idx {
            right = center_idx;
        }
        while right + 1 < coordinates.len() && center.abs_diff(coordinates[right + 1]) <= window_bp
        {
            right += 1;
        }
        ends[center_idx] = right;
    }

    (starts, ends)
}

#[derive(Clone, Debug)]
struct FenwickCounts {
    tree: Vec<usize>,
}

impl FenwickCounts {
    fn from_counts(counts: &[usize]) -> Self {
        let mut result = Self {
            tree: vec![0; counts.len() + 1],
        };
        for (idx, &count) in counts.iter().enumerate() {
            result.add(idx, count);
        }
        result
    }

    fn add(&mut self, idx: usize, value: usize) {
        let mut tree_idx = idx + 1;
        while tree_idx < self.tree.len() {
            self.tree[tree_idx] += value;
            tree_idx += tree_idx & tree_idx.wrapping_neg();
        }
    }

    fn remove(&mut self, idx: usize, value: usize) {
        let mut tree_idx = idx + 1;
        while tree_idx < self.tree.len() {
            self.tree[tree_idx] -= value;
            tree_idx += tree_idx & tree_idx.wrapping_neg();
        }
    }

    fn prefix_sum(&self, end_exclusive: usize) -> usize {
        let mut tree_idx = end_exclusive;
        let mut total = 0usize;
        while tree_idx > 0 {
            total += self.tree[tree_idx];
            tree_idx &= tree_idx - 1;
        }
        total
    }

    fn range_sum(&self, start: usize, end_exclusive: usize) -> usize {
        self.prefix_sum(end_exclusive) - self.prefix_sum(start)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn ties_follow_the_explicit_direction() {
        let coordinates = [100, 110];
        assert_eq!(
            detect_endpoint_modes(coordinates, 0, TieDirection::PreferLower)[0].coordinate(),
            100
        );
        assert_eq!(
            detect_endpoint_modes(coordinates, 0, TieDirection::PreferHigher)[0].coordinate(),
            110
        );
    }
}
