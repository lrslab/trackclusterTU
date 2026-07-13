//! Read input detection, parsing, pooling, and filtering.

use std::collections::{BTreeMap, HashMap};
use std::path::Path;

use crate::io::bed::{BedError, BedParseError};
use crate::io::delimited::{DelimitedError, DelimitedReader, Delimiter};
use crate::model::{Coord, Interval, Strand};
use crate::tools::read_rejections::{
    LocatedRead, ReadLoadOutcome, ReadLocation, ReadRejection, ReadRejectionReason,
};
use crate::tu::multi::{pooled_read_id, SampleManifestRecord};
use crate::tu::ReadRecord;

use super::config::InputFormat;

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(super) enum ResolvedInputFormat {
    Bed6,
    Bed12,
    Tsv,
    Fastq,
}

#[derive(Clone, Debug, PartialEq, Eq)]
struct FullLengthEvidenceRecord {
    contig: String,
    strand: Strand,
    interval: Interval,
    effective_full_length: bool,
}

pub(super) fn cmp_read(a: &ReadRecord, b: &ReadRecord) -> std::cmp::Ordering {
    a.contig
        .cmp(&b.contig)
        .then_with(|| a.strand.cmp(&b.strand))
        .then_with(|| a.interval.start().cmp(&b.interval.start()))
        .then_with(|| a.interval.end().cmp(&b.interval.end()))
        .then_with(|| a.id.cmp(&b.id))
}

impl ResolvedInputFormat {
    const fn as_str(self) -> &'static str {
        match self {
            Self::Bed6 => "bed6",
            Self::Bed12 => "bed12",
            Self::Tsv => "tsv",
            Self::Fastq => "fastq",
        }
    }
}

fn location(
    path: &Path,
    format: ResolvedInputFormat,
    record: Option<u64>,
    line: Option<u64>,
) -> ReadLocation {
    let mut location = ReadLocation::new(path, format.as_str());
    location.record = record;
    location.line = line;
    location
}

fn reject(
    outcome: &mut ReadLoadOutcome,
    location: ReadLocation,
    read_id: Option<String>,
    stage: &'static str,
    reason: ReadRejectionReason,
    detail: impl Into<String>,
) {
    outcome
        .rejections
        .push(ReadRejection::new(location, read_id, stage, reason, detail));
}

fn retain_validated_read(outcome: &mut ReadLoadOutcome, location: ReadLocation, read: ReadRecord) {
    if read.contig.is_empty() {
        reject(
            outcome,
            location,
            Some(read.id),
            "validation",
            ReadRejectionReason::EmptyContig,
            "read contig must not be empty",
        );
        return;
    }
    if read.id.is_empty() {
        reject(
            outcome,
            location,
            None,
            "validation",
            ReadRejectionReason::EmptyReadId,
            "read ID must not be empty",
        );
        return;
    }
    if read.strand == Strand::Unknown {
        let read_id = read.id.clone();
        reject(
            outcome,
            location,
            Some(read_id.clone()),
            "validation",
            ReadRejectionReason::UnknownStrand,
            format!(
                "read {read_id:?} has unknown strand; boundary-aware TU clustering requires '+' or '-'"
            ),
        );
        return;
    }
    outcome.reads.push(LocatedRead {
        read,
        location,
        full_length_evidence: false,
    });
}

fn bed_rejection_reason(error: &BedParseError) -> ReadRejectionReason {
    match error {
        BedParseError::InvalidUtf8 { .. } => ReadRejectionReason::InvalidUtf8,
        BedParseError::TooFewColumnsBed6 { .. } | BedParseError::TooFewColumns { .. } => {
            ReadRejectionReason::TooFewColumns
        }
        BedParseError::InvalidInt { .. } => ReadRejectionReason::InvalidInteger,
        BedParseError::InvalidBlockCount { .. }
        | BedParseError::BlockListLengthMismatch { .. }
        | BedParseError::BlockOverflow
        | BedParseError::Transcript(_) => ReadRejectionReason::InvalidBed12Blocks,
        BedParseError::Strand(_) => ReadRejectionReason::InvalidStrand,
        BedParseError::Interval(_) => ReadRejectionReason::InvalidInterval,
    }
}

fn read_bed6(path: &Path) -> anyhow::Result<ReadLoadOutcome> {
    let mut reader = crate::io::bed::read_bed6(path)?;
    let mut outcome = ReadLoadOutcome::default();
    while let Some(result) = reader.next() {
        outcome.total_records = outcome.total_records.saturating_add(1);
        let location = location(
            path,
            ResolvedInputFormat::Bed6,
            Some(outcome.total_records),
            Some(reader.line_number() as u64),
        );
        match result {
            Ok(record) => retain_validated_read(
                &mut outcome,
                location,
                ReadRecord {
                    contig: record.chrom,
                    strand: record.strand,
                    interval: Interval::new(record.start, record.end)
                        .expect("BED6 parser already validated the interval"),
                    id: record.name,
                },
            ),
            Err(BedError::Parse { source, .. }) => {
                let reason = bed_rejection_reason(&source);
                reject(
                    &mut outcome,
                    location,
                    None,
                    "parse",
                    reason,
                    source.to_string(),
                );
            }
            Err(error) => return Err(error.into()),
        }
    }
    Ok(outcome)
}

fn read_bed12(path: &Path) -> anyhow::Result<ReadLoadOutcome> {
    let mut reader = crate::io::bed::read_bed12(path)?;
    let mut outcome = ReadLoadOutcome::default();
    while let Some(result) = reader.next() {
        outcome.total_records = outcome.total_records.saturating_add(1);
        let location = location(
            path,
            ResolvedInputFormat::Bed12,
            Some(outcome.total_records),
            Some(reader.line_number() as u64),
        );
        match result {
            Ok(tx) => {
                let interval = Interval::new(tx.tx_start(), tx.tx_end())
                    .expect("BED12 parser already validated the transcript span");
                retain_validated_read(
                    &mut outcome,
                    location,
                    ReadRecord {
                        contig: tx.chrom().to_owned(),
                        strand: tx.strand(),
                        interval,
                        id: tx.name().to_owned(),
                    },
                );
            }
            Err(BedError::Parse { source, .. }) => {
                let reason = bed_rejection_reason(&source);
                reject(
                    &mut outcome,
                    location,
                    None,
                    "parse",
                    reason,
                    source.to_string(),
                );
            }
            Err(error) => return Err(error.into()),
        }
    }
    Ok(outcome)
}

fn read_tsv(path: &Path) -> anyhow::Result<ReadLoadOutcome> {
    let mut reader = DelimitedReader::open(path, Delimiter::Tab)?;
    let mut outcome = ReadLoadOutcome::default();

    for table_record in reader.records() {
        let table_record = match table_record {
            Ok(record) => record,
            Err(DelimitedError::Read {
                path: error_path,
                record,
                line,
                source,
            }) if matches!(source.kind(), &csv::ErrorKind::Utf8 { .. }) => {
                outcome.total_records = outcome.total_records.saturating_add(1);
                let location = location(
                    &error_path,
                    ResolvedInputFormat::Tsv,
                    record.or(Some(outcome.total_records)),
                    line,
                );
                reject(
                    &mut outcome,
                    location,
                    None,
                    "parse",
                    ReadRejectionReason::InvalidUtf8,
                    source.to_string(),
                );
                continue;
            }
            Err(error) => return Err(error.into()),
        };
        outcome.total_records = outcome.total_records.saturating_add(1);
        let line_number = table_record.line_number();
        let fields = table_record.fields();
        let location = location(
            path,
            ResolvedInputFormat::Tsv,
            Some(outcome.total_records),
            Some(line_number),
        );
        let read_id = fields
            .get(3)
            .filter(|value| !value.is_empty())
            .map(str::to_owned);
        if fields.len() < 5 {
            reject(
                &mut outcome,
                location,
                read_id,
                "parse",
                ReadRejectionReason::TooFewColumns,
                format!(
                    "expected at least 5 columns (contig, start, end, id, strand), got {}",
                    fields.len()
                ),
            );
            continue;
        }
        let start_text = fields.get(1).unwrap_or("");
        let start: u32 = match start_text.parse() {
            Ok(value) => value,
            Err(_) => {
                reject(
                    &mut outcome,
                    location,
                    read_id,
                    "parse",
                    ReadRejectionReason::InvalidInteger,
                    format!("invalid integer for start: {start_text:?}"),
                );
                continue;
            }
        };
        let end_text = fields.get(2).unwrap_or("");
        let end: u32 = match end_text.parse() {
            Ok(value) => value,
            Err(_) => {
                reject(
                    &mut outcome,
                    location,
                    read_id,
                    "parse",
                    ReadRejectionReason::InvalidInteger,
                    format!("invalid integer for end: {end_text:?}"),
                );
                continue;
            }
        };
        let id = fields.get(3).unwrap_or("").to_owned();
        let strand_text = fields.get(4).unwrap_or("");
        let strand = match Strand::try_from(strand_text) {
            Ok(value) => value,
            Err(error) => {
                reject(
                    &mut outcome,
                    location,
                    Some(id),
                    "parse",
                    ReadRejectionReason::InvalidStrand,
                    format!("invalid strand {strand_text:?}: {error}"),
                );
                continue;
            }
        };
        let interval = match Interval::new(Coord::new(start), Coord::new(end)) {
            Ok(value) => value,
            Err(error) => {
                reject(
                    &mut outcome,
                    location,
                    Some(id),
                    "validation",
                    ReadRejectionReason::InvalidInterval,
                    error.to_string(),
                );
                continue;
            }
        };
        retain_validated_read(
            &mut outcome,
            location,
            ReadRecord {
                contig: fields.get(0).unwrap_or("").to_owned(),
                strand,
                interval,
                id,
            },
        );
    }
    Ok(outcome)
}

fn resolve_explicit_input_format(format: InputFormat) -> Option<ResolvedInputFormat> {
    match format {
        InputFormat::Auto => None,
        InputFormat::Bed6 => Some(ResolvedInputFormat::Bed6),
        InputFormat::Bed12 => Some(ResolvedInputFormat::Bed12),
        InputFormat::Tsv => Some(ResolvedInputFormat::Tsv),
    }
}

pub(super) fn read_reads(
    path: &Path,
    format: ResolvedInputFormat,
) -> anyhow::Result<ReadLoadOutcome> {
    let mut outcome = match format {
        ResolvedInputFormat::Bed6 => read_bed6(path),
        ResolvedInputFormat::Bed12 => read_bed12(path),
        ResolvedInputFormat::Tsv => read_tsv(path),
        ResolvedInputFormat::Fastq => anyhow::bail!(
            "FASTQ inputs are not supported directly; use `trackclustertu map` or `trackclustertu run`, or convert sorted BAMs with `trackclustertu bam-to-bed` before clustering"
        ),
    }?;
    quarantine_duplicate_ids(&mut outcome);
    Ok(outcome)
}

fn parse_evidence_flag(path: &Path, line: u64, column: &str, value: &str) -> anyhow::Result<bool> {
    match value.trim().to_ascii_lowercase().as_str() {
        "1" | "true" => Ok(true),
        "0" | "false" => Ok(false),
        _ => anyhow::bail!(
            "{path:?}:{line}: invalid {column} value {value:?}; expected 0, 1, false, or true"
        ),
    }
}

fn read_full_length_evidence(
    path: &Path,
) -> anyhow::Result<BTreeMap<String, FullLengthEvidenceRecord>> {
    const REQUIRED_COLUMNS: [&str; 8] = [
        "schema_version",
        "read_name",
        "chrom",
        "start",
        "end",
        "strand",
        "retained",
        "effective_full_length",
    ];

    let mut reader = DelimitedReader::open(path, Delimiter::Tab)?;
    let mut records = reader.records();
    let header = records
        .next()
        .transpose()?
        .ok_or_else(|| anyhow::anyhow!("evidence sidecar {path:?} is empty"))?;
    let mut columns = HashMap::new();
    for (index, raw_column) in header.fields().iter().enumerate() {
        let column = raw_column.trim().to_ascii_lowercase();
        if columns.insert(column.clone(), index).is_some() {
            anyhow::bail!(
                "{path:?}:{}: duplicate evidence column {column:?}",
                header.line_number()
            );
        }
    }
    for required in REQUIRED_COLUMNS {
        if !columns.contains_key(required) {
            anyhow::bail!(
                "{path:?}:{}: evidence sidecar is missing required column {required:?}",
                header.line_number()
            );
        }
    }

    let column = |name: &str| -> usize {
        *columns
            .get(name)
            .expect("required evidence columns checked above")
    };
    let required_width = REQUIRED_COLUMNS
        .iter()
        .map(|name| column(name))
        .max()
        .expect("required evidence columns are non-empty")
        + 1;
    let mut evidence = BTreeMap::new();

    for result in records {
        let record = result?;
        let line = record.line_number();
        let fields = record.fields();
        if fields.len() < required_width {
            anyhow::bail!(
                "{path:?}:{line}: expected at least {required_width} evidence columns, got {}",
                fields.len()
            );
        }
        let get = |name: &str| fields.get(column(name)).unwrap_or("").trim();
        let schema = get("schema_version");
        if schema != crate::bam::BAM_EVIDENCE_SCHEMA_VERSION {
            anyhow::bail!(
                "{path:?}:{line}: unsupported evidence schema {schema:?}; expected {:?}",
                crate::bam::BAM_EVIDENCE_SCHEMA_VERSION
            );
        }
        if !parse_evidence_flag(path, line, "retained", get("retained"))? {
            continue;
        }

        let read_name = get("read_name");
        if read_name.is_empty() || read_name == "." {
            anyhow::bail!("{path:?}:{line}: retained evidence row has no read_name");
        }
        let start = get("start").parse::<u32>().map_err(|_| {
            anyhow::anyhow!(
                "{path:?}:{line}: invalid evidence start value {:?}",
                get("start")
            )
        })?;
        let end = get("end").parse::<u32>().map_err(|_| {
            anyhow::anyhow!(
                "{path:?}:{line}: invalid evidence end value {:?}",
                get("end")
            )
        })?;
        let interval = Interval::new(Coord::new(start), Coord::new(end)).map_err(|error| {
            anyhow::anyhow!("{path:?}:{line}: invalid retained evidence interval: {error}")
        })?;
        let strand = Strand::try_from(get("strand")).map_err(|error| {
            anyhow::anyhow!("{path:?}:{line}: invalid evidence strand: {error}")
        })?;
        if strand == Strand::Unknown {
            anyhow::bail!("{path:?}:{line}: retained evidence row has unknown strand");
        }
        let evidence_record = FullLengthEvidenceRecord {
            contig: get("chrom").to_owned(),
            strand,
            interval,
            effective_full_length: parse_evidence_flag(
                path,
                line,
                "effective_full_length",
                get("effective_full_length"),
            )?,
        };
        if evidence
            .insert(read_name.to_owned(), evidence_record)
            .is_some()
        {
            anyhow::bail!("{path:?}:{line}: duplicate retained evidence for read {read_name:?}");
        }
    }

    Ok(evidence)
}

fn attach_sample_evidence(
    sample: &SampleManifestRecord,
    outcome: &mut ReadLoadOutcome,
) -> anyhow::Result<()> {
    let Some(path) = sample.evidence.as_deref() else {
        return Ok(());
    };
    let mut evidence = read_full_length_evidence(path)?;
    for located in &mut outcome.reads {
        let read_id = located.read.id.as_str();
        let evidence_record = evidence.remove(read_id).ok_or_else(|| {
            anyhow::anyhow!(
                "evidence sidecar {path:?} has no retained row for sample {:?} read {read_id:?}",
                sample.sample
            )
        })?;
        if evidence_record.contig != located.read.contig
            || evidence_record.strand != located.read.strand
            || evidence_record.interval != located.read.interval
        {
            anyhow::bail!(
                "evidence sidecar {path:?} does not match sample {:?} read {read_id:?}: evidence={}:{}-{}({}), reads={}:{}-{}({})",
                sample.sample,
                evidence_record.contig,
                evidence_record.interval.start(),
                evidence_record.interval.end(),
                evidence_record.strand.as_char(),
                located.read.contig,
                located.read.interval.start(),
                located.read.interval.end(),
                located.read.strand.as_char(),
            );
        }
        located.full_length_evidence = evidence_record.effective_full_length;
    }
    if let Some((read_id, _)) = evidence.first_key_value() {
        anyhow::bail!(
            "evidence sidecar {path:?} has retained read {read_id:?} that is absent from sample {:?} reads",
            sample.sample
        );
    }
    Ok(())
}

fn quarantine_duplicate_ids(outcome: &mut ReadLoadOutcome) {
    let mut counts: HashMap<String, usize> = HashMap::new();
    for located in &outcome.reads {
        *counts.entry(located.read.id.clone()).or_default() += 1;
    }

    let mut retained = Vec::with_capacity(outcome.reads.len());
    for located in outcome.reads.drain(..) {
        let duplicate_count = counts.get(&located.read.id).copied().unwrap_or(0);
        if duplicate_count > 1 {
            let read_id = located.read.id;
            outcome.rejections.push(ReadRejection::new(
                located.location,
                Some(read_id.clone()),
                "deduplicate",
                ReadRejectionReason::DuplicateReadId,
                format!(
                    "read ID {read_id:?} occurs {duplicate_count} times in this input; all occurrences were skipped"
                ),
            ));
        } else {
            retained.push(located);
        }
    }
    outcome.reads = retained;
}

fn infer_format_from_path(path: &Path) -> Option<ResolvedInputFormat> {
    let file_name = path.file_name()?.to_string_lossy().to_ascii_lowercase();
    if file_name.ends_with(".fastq")
        || file_name.ends_with(".fq")
        || file_name.ends_with(".fastq.gz")
        || file_name.ends_with(".fq.gz")
    {
        return Some(ResolvedInputFormat::Fastq);
    }
    if file_name.ends_with(".bed12") {
        return Some(ResolvedInputFormat::Bed12);
    }
    if file_name.ends_with(".bed") {
        return Some(ResolvedInputFormat::Bed6);
    }
    if file_name.ends_with(".tsv") || file_name.ends_with(".txt") {
        return Some(ResolvedInputFormat::Tsv);
    }
    None
}

pub(super) fn resolve_single_input_format(
    format: InputFormat,
    input: &Path,
) -> ResolvedInputFormat {
    resolve_explicit_input_format(format)
        .or_else(|| infer_format_from_path(input))
        .unwrap_or(ResolvedInputFormat::Bed6)
}

pub(super) fn resolve_manifest_input_format(
    format: InputFormat,
    samples: &[SampleManifestRecord],
) -> anyhow::Result<ResolvedInputFormat> {
    if let Some(format) = resolve_explicit_input_format(format) {
        return Ok(format);
    }
    let mut inferred = None;
    for sample in samples {
        let Some(sample_format) = infer_format_from_path(&sample.reads) else {
            continue;
        };
        if let Some(existing) = inferred {
            if existing != sample_format {
                anyhow::bail!(
                    "auto-detected mixed manifest input formats ({existing:?} and {sample_format:?}); pass --format explicitly"
                );
            }
        } else {
            inferred = Some(sample_format);
        }
    }
    Ok(inferred.unwrap_or(ResolvedInputFormat::Bed6))
}

pub(super) fn load_manifest_reads(
    samples: &[SampleManifestRecord],
    format: ResolvedInputFormat,
) -> anyhow::Result<ReadLoadOutcome> {
    let mut pooled = ReadLoadOutcome::default();
    for sample in samples {
        let mut sample_outcome = read_reads(&sample.reads, format)?;
        attach_sample_evidence(sample, &mut sample_outcome)?;
        sample_outcome.set_sample(&sample.sample);
        for located in &mut sample_outcome.reads {
            located.read.id = pooled_read_id(&sample.sample, &located.read.id);
        }
        pooled.append(sample_outcome);
    }
    Ok(pooled)
}

pub(super) fn filter_reads(outcome: &mut ReadLoadOutcome, min_read_len: Option<u32>) {
    let mut retained = Vec::with_capacity(outcome.reads.len());
    for located in outcome.reads.drain(..) {
        let reason = if located.read.interval.is_empty() {
            Some((
                ReadRejectionReason::ZeroLength,
                "read interval has zero length".to_owned(),
            ))
        } else if let Some(minimum) = min_read_len {
            (located.read.interval.len() < minimum).then(|| {
                (
                    ReadRejectionReason::BelowMinReadLen,
                    format!(
                        "read length {} is below --min-read-len {minimum}",
                        located.read.interval.len()
                    ),
                )
            })
        } else {
            None
        };

        if let Some((reason, detail)) = reason {
            outcome.rejections.push(ReadRejection::new(
                located.location,
                Some(located.read.id),
                "filter",
                reason,
                detail,
            ));
        } else {
            retained.push(located);
        }
    }
    outcome.reads = retained;
}

#[cfg(test)]
mod tests {
    use std::sync::atomic::{AtomicU64, Ordering};

    use super::*;

    static NEXT_TMP_ID: AtomicU64 = AtomicU64::new(0);

    fn write_temp_bytes(label: &str, extension: &str, contents: &[u8]) -> std::path::PathBuf {
        let path = std::env::temp_dir().join(format!(
            "trackclustertu_input_{label}_{}_{}.{}",
            std::process::id(),
            NEXT_TMP_ID.fetch_add(1, Ordering::Relaxed),
            extension
        ));
        std::fs::write(&path, contents).unwrap();
        path
    }

    fn write_temp_tsv(label: &str, contents: &str) -> std::path::PathBuf {
        write_temp_bytes(label, "tsv", contents.as_bytes())
    }

    #[test]
    fn tsv_unknown_strand_is_quarantined_with_path_and_physical_line() {
        let path = write_temp_tsv(
            "unknown_strand",
            "# metadata\n\nchr1\t0\t10\tknown\t+\nchr1\t20\t30\tunknown\t.\n",
        );

        let outcome = read_tsv(&path).unwrap();
        assert_eq!(outcome.total_records, 2);
        assert_eq!(outcome.reads.len(), 1);
        assert_eq!(outcome.reads[0].read.id, "known");
        assert_eq!(outcome.rejections.len(), 1);
        let rejection = &outcome.rejections[0];
        assert_eq!(rejection.location.source_path, path);
        assert_eq!(rejection.location.line, Some(4));
        assert_eq!(rejection.read_id.as_deref(), Some("unknown"));
        assert_eq!(rejection.stage, "validation");
        assert_eq!(rejection.reason, ReadRejectionReason::UnknownStrand);

        let _ = std::fs::remove_file(path);
    }

    #[test]
    fn tsv_invalid_interval_is_quarantined_with_physical_line() {
        let path = write_temp_tsv(
            "invalid_interval",
            "# metadata\nchr1\t0\t10\t\"multi\nline\"\t+\nchr1\t30\t20\tinvalid\t+\n",
        );

        let outcome = read_tsv(&path).unwrap();
        assert_eq!(outcome.total_records, 2);
        assert_eq!(outcome.reads.len(), 1);
        assert_eq!(outcome.reads[0].read.id, "multi\nline");
        assert_eq!(outcome.rejections.len(), 1);
        let rejection = &outcome.rejections[0];
        assert_eq!(rejection.location.source_path, path);
        assert_eq!(rejection.location.line, Some(4));
        assert_eq!(rejection.read_id.as_deref(), Some("invalid"));
        assert_eq!(rejection.reason, ReadRejectionReason::InvalidInterval);

        let _ = std::fs::remove_file(path);
    }

    #[test]
    fn bed6_good_bad_good_rows_continue_and_keep_locations() {
        let path = write_temp_bytes(
            "bed6_continue",
            "bed",
            b"chr1\t0\t10\tgood1\t0\t+\nmalformed\nchr1\t20\t30\tgood2\t0\t+\n",
        );

        let outcome = read_reads(&path, ResolvedInputFormat::Bed6).unwrap();
        assert_eq!(outcome.total_records, 3);
        assert_eq!(
            outcome
                .reads
                .iter()
                .map(|located| located.read.id.as_str())
                .collect::<Vec<_>>(),
            ["good1", "good2"]
        );
        assert_eq!(outcome.rejections.len(), 1);
        assert_eq!(outcome.rejections[0].location.line, Some(2));
        assert_eq!(
            outcome.rejections[0].reason,
            ReadRejectionReason::TooFewColumns
        );

        let _ = std::fs::remove_file(path);
    }

    #[test]
    fn bed12_good_bad_good_rows_continue() {
        let path = write_temp_bytes(
            "bed12_continue",
            "bed12",
            concat!(
                "chr1\t0\t10\tgood1\t0\t+\t0\t10\t0\t1\t10,\t0,\n",
                "chr1\t20\t30\tbad\t0\t+\t20\t30\t0\t1\t0,\t0,\n",
                "chr1\t40\t50\tgood2\t0\t+\t40\t50\t0\t1\t10,\t0,\n",
            )
            .as_bytes(),
        );

        let outcome = read_reads(&path, ResolvedInputFormat::Bed12).unwrap();
        assert_eq!(outcome.total_records, 3);
        assert_eq!(
            outcome
                .reads
                .iter()
                .map(|located| located.read.id.as_str())
                .collect::<Vec<_>>(),
            ["good1", "good2"]
        );
        assert_eq!(outcome.rejections.len(), 1);
        assert_eq!(outcome.rejections[0].location.line, Some(2));
        assert_eq!(
            outcome.rejections[0].reason,
            ReadRejectionReason::InvalidBed12Blocks
        );

        let _ = std::fs::remove_file(path);
    }

    #[test]
    fn tsv_invalid_utf8_record_does_not_hide_later_good_record() {
        let mut input = b"chr1\t0\t10\tgood1\t+\nchr1\t20\t30\t".to_vec();
        input.push(0xff);
        input.extend_from_slice(b"\t+\nchr1\t40\t50\tgood2\t+\n");
        let path = write_temp_bytes("tsv_utf8", "tsv", &input);

        let outcome = read_reads(&path, ResolvedInputFormat::Tsv).unwrap();
        assert_eq!(outcome.total_records, 3);
        assert_eq!(
            outcome
                .reads
                .iter()
                .map(|located| located.read.id.as_str())
                .collect::<Vec<_>>(),
            ["good1", "good2"]
        );
        assert_eq!(outcome.rejections.len(), 1);
        assert_eq!(outcome.rejections[0].location.line, Some(2));
        assert_eq!(
            outcome.rejections[0].reason,
            ReadRejectionReason::InvalidUtf8
        );

        let _ = std::fs::remove_file(path);
    }

    #[test]
    fn duplicate_ids_quarantine_every_occurrence() {
        let path = write_temp_bytes(
            "duplicate_ids",
            "bed",
            b"chr1\t0\t10\tdup\t0\t+\nchr1\t20\t30\tkeep\t0\t+\nchr1\t40\t50\tdup\t0\t+\n",
        );

        let outcome = read_reads(&path, ResolvedInputFormat::Bed6).unwrap();
        assert_eq!(outcome.reads.len(), 1);
        assert_eq!(outcome.reads[0].read.id, "keep");
        assert_eq!(outcome.rejections.len(), 2);
        assert!(outcome.rejections.iter().all(|rejection| {
            rejection.reason == ReadRejectionReason::DuplicateReadId
                && rejection.read_id.as_deref() == Some("dup")
        }));
        assert_eq!(
            outcome
                .rejections
                .iter()
                .filter_map(|rejection| rejection.location.line)
                .collect::<Vec<_>>(),
            [1, 3]
        );

        let _ = std::fs::remove_file(path);
    }

    #[test]
    fn empty_contig_and_read_id_are_quarantined_without_stopping() {
        let path = write_temp_bytes(
            "empty_fields",
            "bed",
            b"\t0\t10\tempty_contig\t0\t+\nchr1\t20\t30\t\t0\t+\nchr1\t40\t50\tkeep\t0\t+\n",
        );

        let outcome = read_reads(&path, ResolvedInputFormat::Bed6).unwrap();
        assert_eq!(outcome.total_records, 3);
        assert_eq!(outcome.reads.len(), 1);
        assert_eq!(outcome.reads[0].read.id, "keep");
        assert_eq!(outcome.rejections.len(), 2);
        assert_eq!(
            outcome.rejections[0].reason,
            ReadRejectionReason::EmptyContig
        );
        assert_eq!(
            outcome.rejections[1].reason,
            ReadRejectionReason::EmptyReadId
        );

        let _ = std::fs::remove_file(path);
    }

    #[test]
    fn zero_length_and_short_reads_are_recorded_by_filter() {
        let path = write_temp_bytes(
            "filters",
            "bed",
            b"chr1\t0\t0\tzero\t0\t+\nchr1\t10\t15\tshort\t0\t+\nchr1\t20\t40\tkeep\t0\t+\n",
        );

        let mut outcome = read_reads(&path, ResolvedInputFormat::Bed6).unwrap();
        filter_reads(&mut outcome, Some(10));
        assert_eq!(outcome.reads.len(), 1);
        assert_eq!(outcome.reads[0].read.id, "keep");
        assert_eq!(outcome.rejections.len(), 2);
        assert_eq!(
            outcome.rejections[0].reason,
            ReadRejectionReason::ZeroLength
        );
        assert_eq!(
            outcome.rejections[1].reason,
            ReadRejectionReason::BelowMinReadLen
        );

        let _ = std::fs::remove_file(path);
    }

    #[test]
    fn identical_raw_ids_in_different_samples_remain_valid() {
        let first = write_temp_bytes("sample_a", "bed", b"chr1\t0\t10\tsame\t0\t+\n");
        let second = write_temp_bytes("sample_b", "bed", b"chr1\t20\t30\tsame\t0\t+\n");
        let samples = vec![
            SampleManifestRecord {
                sample: "sampleA".to_owned(),
                reads: first.clone(),
                group: None,
                evidence: None,
            },
            SampleManifestRecord {
                sample: "sampleB".to_owned(),
                reads: second.clone(),
                group: None,
                evidence: None,
            },
        ];

        let outcome = load_manifest_reads(&samples, ResolvedInputFormat::Bed6).unwrap();
        assert!(outcome.rejections.is_empty());
        assert_eq!(
            outcome
                .reads
                .iter()
                .map(|located| located.read.id.as_str())
                .collect::<Vec<_>>(),
            ["sampleA::same", "sampleB::same"]
        );
        assert_eq!(
            outcome
                .reads
                .iter()
                .filter_map(|located| located.location.sample.as_deref())
                .collect::<Vec<_>>(),
            ["sampleA", "sampleB"]
        );

        let _ = std::fs::remove_file(first);
        let _ = std::fs::remove_file(second);
    }

    #[test]
    fn manifest_evidence_sidecar_is_attached_to_exact_retained_reads() {
        let reads = write_temp_bytes(
            "evidence_reads",
            "bed",
            b"chr1\t0\t10\tr1\t0\t+\nchr1\t20\t30\tr2\t0\t-\n",
        );
        let evidence = write_temp_bytes(
            "evidence_sidecar",
            "tsv",
            concat!(
                "schema_version\tread_name\tchrom\tstart\tend\tstrand\tretained\teffective_full_length\n",
                "trackclustertu.bam-evidence.v1\tr1\tchr1\t0\t10\t+\t1\t1\n",
                "trackclustertu.bam-evidence.v1\tr2\tchr1\t20\t30\t-\t1\t0\n",
                "trackclustertu.bam-evidence.v1\tfiltered\t.\t.\t.\t.\t0\t0\n",
            )
            .as_bytes(),
        );
        let samples = vec![SampleManifestRecord {
            sample: "sampleA".to_owned(),
            reads: reads.clone(),
            group: None,
            evidence: Some(evidence.clone()),
        }];

        let outcome = load_manifest_reads(&samples, ResolvedInputFormat::Bed6).unwrap();
        assert_eq!(outcome.reads.len(), 2);
        assert_eq!(outcome.reads[0].read.id, "sampleA::r1");
        assert!(outcome.reads[0].full_length_evidence);
        assert_eq!(outcome.reads[1].read.id, "sampleA::r2");
        assert!(!outcome.reads[1].full_length_evidence);

        let _ = std::fs::remove_file(reads);
        let _ = std::fs::remove_file(evidence);
    }
}
