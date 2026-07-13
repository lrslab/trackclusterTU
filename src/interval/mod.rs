//! Interval-list, partitioning, intersection, and transcript-span algorithms.

mod cluster_span;
mod intersect_sweep;
mod list_ops;
mod partition;
mod refine;
mod sort;

use crate::model::Strand;

#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
/// Whether an operation groups or matches records independently by strand.
pub enum StrandMode {
    /// Ignore strand when constructing partition keys.
    #[default]
    Ignore,
    /// Require matching strands.
    Match,
}

impl StrandMode {
    /// Convert a record strand into the optional strand component of a key.
    pub fn key_strand(self, strand: Strand) -> Option<Strand> {
        match self {
            Self::Ignore => None,
            Self::Match => Some(strand),
        }
    }
}

pub use cluster_span::{cluster_by_span, RangeCluster};
pub use intersect_sweep::{sweep_intersect_pairs, IntersectOpts};
pub use list_ops::{intersection_len, merge_overlaps, total_len, union_len};
pub use partition::{partition, PartitionKey};
pub use refine::{exonic_overlap_bp, junctions_equal, junctions_subset};
pub use sort::sort_by_coord;
