//! Output helpers for legacy transcript clusters.

use std::path::Path;

use crate::io::bed::BedError;
use crate::io::delimited::{DelimitedWriter, Delimiter};
use crate::model::Transcript;

/// Create a headerless TSV containing `(read_id, isoform_id)` pairs.
///
/// Rows retain input order and fields use standards-compliant TSV quoting.
///
/// # Errors
///
/// Returns an I/O error if `path` cannot be created or a row cannot be written.
pub fn write_read_to_isoform_tsv<P: AsRef<Path>>(
    path: P,
    pairs: &[(String, String)],
) -> Result<(), std::io::Error> {
    let mut writer = DelimitedWriter::create(path.as_ref(), Delimiter::Tab, &[])
        .map_err(std::io::Error::other)?;
    for (read, isoform) in pairs {
        writer
            .write_record([read.as_str(), isoform.as_str()])
            .map_err(std::io::Error::other)?;
    }
    writer.flush().map_err(std::io::Error::other)?;
    Ok(())
}

/// Create a BED12 file containing representative isoforms in input order.
///
/// # Errors
///
/// Returns [`BedError`] if the file cannot be created or written.
pub fn write_isoforms_bed<P: AsRef<Path>>(
    path: P,
    isoforms: &[Transcript],
) -> Result<(), BedError> {
    crate::io::bed::write_bed12(path, isoforms.iter())
}
