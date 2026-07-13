//! Core data structures and algorithms for calling bacterial transcription units.
//!
//! The `trackclustertu` binary combines the public interval, alignment, clustering,
//! and counting APIs exposed here. Command implementation details remain private;
//! applications embedding the library should use the typed domain APIs instead.
#![deny(rustdoc::broken_intra_doc_links)]
#![warn(missing_docs)]

/// BAM-to-BED conversion and evidence-aware filtering.
pub mod bam;
/// Legacy transcript representative and output helpers.
pub mod cluster;
/// Legacy subread-counting helpers.
pub mod count;
/// Indexed endpoint support-mode detection.
pub mod endpoint;
/// Interval and transcript-span algorithms.
pub mod interval;
/// BED and GFF readers and writers.
pub mod io;
/// Checked coordinate, interval, strand, and transcript models.
pub mod model;
/// Span and exon-list similarity scores.
pub mod score;
/// Transcription-unit clustering and assignment.
pub mod tu;

mod tools;

/// Run the `trackclustertu` command-line application using the process arguments.
///
/// This is the deliberately small facade used by the packaged binary. Library
/// integrations should prefer the typed APIs in modules such as [`tu`] and [`bam`].
pub fn cli_entrypoint() -> anyhow::Result<()> {
    tools::trackclustertu::entrypoint()
}
