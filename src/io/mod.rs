//! Readers and writers for genomic interchange formats.

pub(crate) mod delimited;

/// BED6 and BED12 streaming I/O.
pub mod bed;
/// GFF3 gene extraction and BED6 conversion.
pub mod gff;
