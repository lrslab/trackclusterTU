use std::collections::BTreeMap;
use std::fs::{self, File};
use std::io::{BufWriter, Write};
use std::path::{Path, PathBuf};

use anyhow::Context;

use crate::tu::multi::{parse_manifest, SampleManifestRecord};

pub(crate) fn sanitize_sample_name(sample: &str) -> String {
    sample
        .chars()
        .map(|ch| match ch {
            'A'..='Z' | 'a'..='z' | '0'..='9' | '.' | '_' | '-' => ch,
            _ => '_',
        })
        .collect()
}

pub(crate) fn validate_unique_sanitized_sample_names(
    records: &[SampleManifestRecord],
) -> anyhow::Result<()> {
    let mut collisions: BTreeMap<String, Vec<String>> = BTreeMap::new();
    for record in records {
        collisions
            .entry(sanitize_sample_name(&record.sample))
            .or_default()
            .push(record.sample.clone());
    }

    let duplicates: Vec<(String, Vec<String>)> = collisions
        .into_iter()
        .filter_map(|(sanitized, mut samples)| {
            if samples.len() < 2 {
                return None;
            }
            samples.sort();
            samples.dedup();
            Some((sanitized, samples))
        })
        .collect();

    if duplicates.is_empty() {
        return Ok(());
    }

    let details = duplicates
        .into_iter()
        .map(|(sanitized, samples)| format!("{sanitized:?} <- {}", samples.join(", ")))
        .collect::<Vec<_>>()
        .join("; ");

    anyhow::bail!(
        "sample names must remain unique after filename sanitization to avoid overwriting outputs: {details}"
    )
}

pub(crate) fn convert_single_bam_to_bed(input_bam: &Path, out_bed: &Path) -> anyhow::Result<()> {
    crate::bam::bam_to_bed6(input_bam, out_bed)
}

pub(crate) fn convert_records_to_bed_manifest(
    records: &[SampleManifestRecord],
    out_dir: &Path,
) -> anyhow::Result<PathBuf> {
    validate_unique_sanitized_sample_names(records)?;

    let bed_dir = out_dir.join("bed");
    fs::create_dir_all(&bed_dir)
        .with_context(|| format!("failed to create BED output directory {:?}", bed_dir))?;

    let manifest_path = out_dir.join("samples.bed.tsv");
    let mut writer = BufWriter::new(
        File::create(&manifest_path)
            .with_context(|| format!("failed to create BED manifest {:?}", manifest_path))?,
    );
    writeln!(writer, "sample\tgroup\treads")?;

    for record in records {
        let sample_stem = sanitize_sample_name(&record.sample);
        let bed_path = bed_dir.join(format!("{sample_stem}.bed"));
        convert_single_bam_to_bed(&record.reads, &bed_path).with_context(|| {
            format!(
                "failed to convert BAM {:?} for sample {:?}",
                record.reads, record.sample
            )
        })?;
        let bed_path = bed_path
            .canonicalize()
            .with_context(|| format!("failed to canonicalize BED path {:?}", bed_path))?;

        writeln!(
            writer,
            "{}\t{}\t{}",
            record.sample,
            record.group.as_deref().unwrap_or(""),
            bed_path.display()
        )?;
    }

    Ok(manifest_path)
}

pub(crate) fn convert_manifest_to_bed_manifest(
    manifest: &Path,
    out_dir: &Path,
) -> anyhow::Result<PathBuf> {
    let records = parse_manifest(manifest)
        .with_context(|| format!("failed to parse BAM manifest {:?}", manifest))?;
    convert_records_to_bed_manifest(&records, out_dir)
}
