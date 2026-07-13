//! Checked genomic coordinate and transcript domain models.

/// Unsigned genomic coordinates.
mod coord;
/// Validated half-open intervals.
mod interval;
/// Reference-strand values and parsing.
mod strand;
/// BED12-compatible transcript models.
mod transcript;

pub use coord::Coord;
pub use interval::{Interval, IntervalError};
pub use strand::{Strand, StrandParseError};
pub use transcript::{Bed12Attrs, JunctionSignature, Transcript, TranscriptError};
