use std::collections::BTreeMap;
use std::path::{Path, PathBuf};

use anyhow::Context;

use crate::bam::{BamConversionConfig, BamConversionSummary, LibraryProfile};
use crate::io::delimited::{DelimitedReader, DelimitedWriter, Delimiter};
use crate::tools::output_transaction::OutputTransaction;
use crate::tu::multi::{parse_manifest, SampleManifestRecord};

#[derive(Clone, Debug, Default)]
struct SampleEvidenceSettings {
    input_evidence: Option<PathBuf>,
    library_profile: Option<LibraryProfile>,
}

#[derive(Clone, Debug)]
pub(crate) struct ManifestConversionResult {
    pub bed_manifest: PathBuf,
    pub summary: BamConversionSummary,
}

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

pub(crate) fn convert_single_bam_to_bed_with_config(
    input_bam: &Path,
    out_bed: &Path,
    out_evidence: Option<&Path>,
    input_evidence: Option<&Path>,
    config: &BamConversionConfig,
) -> anyhow::Result<BamConversionSummary> {
    Ok(crate::bam::bam_to_bed6_with_evidence(
        input_bam,
        out_bed,
        out_evidence,
        input_evidence,
        config,
    )?)
}

pub(crate) fn convert_records_to_bed_manifest(
    records: &[SampleManifestRecord],
    out_dir: &Path,
) -> anyhow::Result<PathBuf> {
    Ok(convert_records_to_bed_manifest_with_options(
        records,
        out_dir,
        &[],
        &BTreeMap::new(),
        &BamConversionConfig::default(),
        false,
    )?
    .bed_manifest)
}

fn convert_records_to_bed_manifest_with_options(
    records: &[SampleManifestRecord],
    out_dir: &Path,
    additional_inputs: &[PathBuf],
    per_sample_settings: &BTreeMap<String, SampleEvidenceSettings>,
    config: &BamConversionConfig,
    emit_evidence: bool,
) -> anyhow::Result<ManifestConversionResult> {
    validate_unique_sanitized_sample_names(records)?;

    let bed_dir = out_dir.join("bed");
    let manifest_path = out_dir.join("samples.bed.tsv");
    let bed_paths: Vec<PathBuf> = records
        .iter()
        .map(|record| {
            let sample_stem = sanitize_sample_name(&record.sample);
            bed_dir.join(format!("{sample_stem}.bed"))
        })
        .collect();
    let evidence_paths: Vec<Option<PathBuf>> = bed_paths
        .iter()
        .map(|bed_path| {
            emit_evidence.then(|| {
                let stem = bed_path
                    .file_stem()
                    .and_then(|value| value.to_str())
                    .unwrap_or("sample");
                bed_path.with_file_name(format!("{stem}.evidence.tsv"))
            })
        })
        .collect();
    let mut inputs: Vec<PathBuf> = records.iter().map(|record| record.reads.clone()).collect();
    inputs.extend_from_slice(additional_inputs);
    inputs.extend(
        per_sample_settings
            .values()
            .filter_map(|settings| settings.input_evidence.clone()),
    );
    let mut outputs = bed_paths.clone();
    outputs.extend(evidence_paths.iter().flatten().cloned());
    outputs.push(manifest_path.clone());
    let transaction = OutputTransaction::new(inputs, outputs)?;

    let staged_manifest = transaction.staged_path(&manifest_path)?;
    let mut writer = DelimitedWriter::create(&staged_manifest, Delimiter::Tab, &[])?;
    if emit_evidence {
        writer.write_record(["sample", "group", "reads", "evidence", "library_profile"])?;
    } else {
        writer.write_record(["sample", "group", "reads"])?;
    }

    let mut summary = BamConversionSummary::default();
    for ((record, bed_path), evidence_path) in records.iter().zip(&bed_paths).zip(&evidence_paths) {
        let staged_bed = transaction.staged_path(bed_path)?;
        let settings = per_sample_settings.get(&record.sample);
        let mut sample_config = config.clone();
        if let Some(profile) = settings.and_then(|value| value.library_profile) {
            sample_config.library_profile = profile;
        }
        let staged_evidence = evidence_path
            .as_ref()
            .map(|path| transaction.staged_path(path))
            .transpose()?;
        let sample_summary = convert_single_bam_to_bed_with_config(
            &record.reads,
            &staged_bed,
            staged_evidence.as_deref(),
            settings.and_then(|value| value.input_evidence.as_deref()),
            &sample_config,
        )
        .with_context(|| {
            format!(
                "failed to convert BAM {:?} for sample {:?}",
                record.reads, record.sample
            )
        })?;
        summary.add_assign(&sample_summary);
        let parent = bed_path.parent().unwrap_or_else(|| Path::new("."));
        let final_bed_path = parent
            .canonicalize()
            .with_context(|| format!("failed to canonicalize BED directory {parent:?}"))?
            .join(
                bed_path
                    .file_name()
                    .expect("generated BED output always has a file name"),
            );

        if let Some(evidence_path) = evidence_path {
            let evidence_parent = evidence_path.parent().unwrap_or_else(|| Path::new("."));
            let final_evidence_path = evidence_parent
                .canonicalize()
                .with_context(|| {
                    format!("failed to canonicalize evidence directory {evidence_parent:?}")
                })?
                .join(
                    evidence_path
                        .file_name()
                        .expect("generated evidence output always has a file name"),
                );
            writer.write_record([
                record.sample.clone(),
                record.group.clone().unwrap_or_default(),
                final_bed_path.display().to_string(),
                final_evidence_path.display().to_string(),
                sample_config.library_profile.to_string(),
            ])?;
        } else {
            writer.write_record([
                record.sample.clone(),
                record.group.clone().unwrap_or_default(),
                final_bed_path.display().to_string(),
            ])?;
        }
    }
    writer.flush()?;
    drop(writer);
    transaction.commit()?;

    Ok(ManifestConversionResult {
        bed_manifest: manifest_path,
        summary,
    })
}

pub(crate) fn convert_manifest_to_bed_manifest_with_config(
    manifest: &Path,
    out_dir: &Path,
    config: &BamConversionConfig,
    emit_evidence: bool,
) -> anyhow::Result<ManifestConversionResult> {
    let records = parse_manifest(manifest)
        .with_context(|| format!("failed to parse BAM manifest {:?}", manifest))?;
    let settings = parse_sample_evidence_settings(manifest)?;
    convert_records_to_bed_manifest_with_options(
        &records,
        out_dir,
        std::slice::from_ref(&manifest.to_path_buf()),
        &settings,
        config,
        emit_evidence,
    )
}

fn parse_sample_evidence_settings(
    manifest: &Path,
) -> anyhow::Result<BTreeMap<String, SampleEvidenceSettings>> {
    let mut reader = DelimitedReader::open(manifest, Delimiter::Tab)?;
    let base_dir = manifest.parent().unwrap_or_else(|| Path::new("."));
    let mut columns: Option<Vec<String>> = None;
    let mut settings = BTreeMap::new();

    for result in reader.records() {
        let record = result?;
        let line_number = record.line_number() as usize;
        let fields = record.fields();
        if columns.is_none() {
            columns = Some(
                fields
                    .iter()
                    .map(|value| value.trim().to_ascii_lowercase())
                    .collect(),
            );
            continue;
        }

        let header = columns.as_ref().expect("initialized above");
        let get = |name: &str| -> &str {
            header
                .iter()
                .position(|column| column == name)
                .and_then(|index| fields.get(index))
                .unwrap_or("")
                .trim()
        };
        let sample = get("sample");
        if sample.is_empty() {
            continue;
        }
        let evidence_value = {
            let preferred = get("evidence");
            if preferred.is_empty() {
                get("evidence_tsv")
            } else {
                preferred
            }
        };
        let input_evidence = if evidence_value.is_empty() || evidence_value == "." {
            None
        } else {
            let path = PathBuf::from(evidence_value);
            Some(if path.is_absolute() {
                path
            } else {
                base_dir.join(path)
            })
        };
        let library_profile = match get("library_profile") {
            "" | "." => None,
            value => Some(value.parse::<LibraryProfile>().with_context(|| {
                format!(
                    "invalid library_profile in BAM manifest {:?}:{}",
                    manifest, line_number
                )
            })?),
        };
        if input_evidence.is_some() || library_profile.is_some() {
            settings.insert(
                sample.to_owned(),
                SampleEvidenceSettings {
                    input_evidence,
                    library_profile,
                },
            );
        }
    }
    Ok(settings)
}
