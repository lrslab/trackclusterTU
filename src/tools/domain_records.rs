//! Shared BED6-to-domain loaders used by clustering, diagnosis, and rescue.

use std::collections::HashMap;
use std::io::BufRead;
use std::path::Path;

use anyhow::Context;

use crate::io::gff::GeneRecord;
use crate::model::{Interval, Strand};

#[derive(Clone, Copy, Debug)]
pub(crate) struct Bed6Validation {
    pub(crate) id_label: &'static str,
    pub(crate) require_nonempty_id: bool,
    pub(crate) require_unique_ids: bool,
    pub(crate) unknown_strand_context: Option<&'static str>,
    pub(crate) invalid_interval_label: &'static str,
}

#[derive(Clone, Debug)]
pub(crate) struct NamedIntervalRecord {
    pub(crate) id: String,
    pub(crate) contig: String,
    pub(crate) strand: Strand,
    pub(crate) interval: Interval,
}

pub(crate) fn read_named_bed6(
    path: &Path,
    validation: Bed6Validation,
) -> anyhow::Result<Vec<NamedIntervalRecord>> {
    let line_numbers = bed6_data_line_numbers(path)?;
    let reader = crate::io::bed::read_bed6(path)
        .with_context(|| format!("failed to open {} BED6 {path:?}", validation.id_label))?;
    let mut records = Vec::new();
    let mut first_line_by_id = HashMap::new();
    for (record_index, record) in reader.enumerate() {
        let line_number = line_numbers
            .get(record_index)
            .copied()
            .unwrap_or(record_index + 1);
        let record = record
            .with_context(|| format!("failed to parse {} BED6 {path:?}", validation.id_label))?;
        if validation.require_nonempty_id && record.name.is_empty() {
            anyhow::bail!(
                "{path:?}:{line_number}: {} ID must not be empty",
                validation.id_label
            );
        }
        if validation.require_unique_ids {
            if let Some(first_line) = first_line_by_id.insert(record.name.clone(), line_number) {
                anyhow::bail!(
                    "{path:?}:{line_number}: duplicate {} ID {:?}; first defined at line {first_line}",
                    validation.id_label,
                    record.name
                );
            }
        }
        if let Some(context) = validation.unknown_strand_context {
            if record.strand == Strand::Unknown {
                anyhow::bail!(
                    "{path:?}:{line_number}: {} {:?} has unknown strand; {context}",
                    validation.id_label,
                    record.name
                );
            }
        }
        let interval = Interval::new(record.start, record.end).with_context(|| {
            format!(
                "{path:?}:{line_number}: invalid {} interval",
                validation.invalid_interval_label
            )
        })?;
        records.push(NamedIntervalRecord {
            id: record.name,
            contig: record.chrom,
            strand: record.strand,
            interval,
        });
    }
    Ok(records)
}

pub(crate) fn read_bed6_as_genes(
    path: &Path,
    validation: Bed6Validation,
) -> anyhow::Result<Vec<GeneRecord>> {
    read_named_bed6(path, validation).map(|records| {
        records
            .into_iter()
            .map(|record| GeneRecord {
                contig: record.contig,
                strand: record.strand,
                interval: record.interval,
                id: record.id,
                feature_kind: "gene".to_owned(),
            })
            .collect()
    })
}

fn bed6_data_line_numbers(path: &Path) -> anyhow::Result<Vec<usize>> {
    let file = std::fs::File::open(path)
        .with_context(|| format!("failed to read BED6 line numbers from {path:?}"))?;
    let reader = std::io::BufReader::new(file);
    let mut line_numbers = Vec::new();
    for (line_index, line) in reader.lines().enumerate() {
        let line = line.with_context(|| {
            format!("failed to read BED6 line {} from {path:?}", line_index + 1)
        })?;
        let trimmed = line.trim();
        if !trimmed.is_empty() && !trimmed.starts_with('#') {
            line_numbers.push(line_index + 1);
        }
    }
    Ok(line_numbers)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn configurable_loader_preserves_duplicate_line_context() {
        let root = std::env::temp_dir().join(format!(
            "trackclustertu_shared_bed_loader_{}",
            std::time::SystemTime::now()
                .duration_since(std::time::UNIX_EPOCH)
                .unwrap()
                .as_nanos()
        ));
        std::fs::create_dir_all(&root).unwrap();
        let path = root.join("reads.bed");
        std::fs::write(
            &path,
            "# comment\nchr1\t0\t10\tr1\t0\t+\nchr1\t20\t30\tr1\t0\t+\n",
        )
        .unwrap();
        let error = read_named_bed6(
            &path,
            Bed6Validation {
                id_label: "read",
                require_nonempty_id: true,
                require_unique_ids: true,
                unknown_strand_context: Some("known strand required"),
                invalid_interval_label: "BED6",
            },
        )
        .unwrap_err();
        let message = error.to_string();
        assert!(message.contains(":3:"), "{message}");
        assert!(message.contains("first defined at line 2"), "{message}");
        let _ = std::fs::remove_dir_all(root);
    }
}
