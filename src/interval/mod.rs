pub mod cluster_span;
pub mod intersect_sweep;
pub mod list_ops;
pub mod partition;
pub mod refine;
pub mod sort;

use crate::model::Strand;

#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub enum StrandMode {
    #[default]
    Ignore,
    Match,
}

impl StrandMode {
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
