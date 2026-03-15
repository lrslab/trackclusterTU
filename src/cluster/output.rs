use std::io::Write;
use std::path::Path;

use crate::io::bed::BedError;
use crate::model::Transcript;

pub fn write_read_to_isoform_tsv<P: AsRef<Path>>(
    path: P,
    pairs: &[(String, String)],
) -> Result<(), std::io::Error> {
    let mut writer = std::io::BufWriter::new(std::fs::File::create(path)?);
    for (read, isoform) in pairs {
        writeln!(writer, "{read}\t{isoform}")?;
    }
    Ok(())
}

pub fn write_isoforms_bed<P: AsRef<Path>>(
    path: P,
    isoforms: &[Transcript],
) -> Result<(), BedError> {
    crate::io::bed::write_bed12(path, isoforms.iter())
}
