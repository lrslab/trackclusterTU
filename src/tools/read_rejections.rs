//! Structured, transaction-friendly diagnostics for reads excluded from a run.

use std::collections::BTreeMap;
use std::path::{Path, PathBuf};

use crate::io::delimited::{DelimitedWriter, Delimiter};
use crate::tu::ReadRecord;

pub(crate) const READ_REJECTIONS_SCHEMA_VERSION: &str = "trackclustertu.read-rejections.v1";

#[derive(Clone, Copy, Debug, Eq, Ord, PartialEq, PartialOrd)]
pub(crate) enum ReadRejectionReason {
    TooFewColumns,
    InvalidInteger,
    InvalidStrand,
    UnknownStrand,
    InvalidInterval,
    InvalidBed12Blocks,
    InvalidUtf8,
    EmptyContig,
    EmptyReadId,
    ZeroLength,
    BelowMinReadLen,
    DuplicateReadId,
}

impl ReadRejectionReason {
    pub(crate) const fn as_str(self) -> &'static str {
        match self {
            Self::TooFewColumns => "too_few_columns",
            Self::InvalidInteger => "invalid_integer",
            Self::InvalidStrand => "invalid_strand",
            Self::UnknownStrand => "unknown_strand",
            Self::InvalidInterval => "invalid_interval",
            Self::InvalidBed12Blocks => "invalid_bed12_blocks",
            Self::InvalidUtf8 => "invalid_utf8",
            Self::EmptyContig => "empty_contig",
            Self::EmptyReadId => "empty_read_id",
            Self::ZeroLength => "zero_length",
            Self::BelowMinReadLen => "below_min_read_len",
            Self::DuplicateReadId => "duplicate_read_id",
        }
    }
}

#[derive(Clone, Debug, Eq, PartialEq)]
pub(crate) struct ReadLocation {
    pub(crate) source_path: PathBuf,
    pub(crate) sample: Option<String>,
    pub(crate) input_format: &'static str,
    pub(crate) record: Option<u64>,
    pub(crate) line: Option<u64>,
}

impl ReadLocation {
    pub(crate) fn new(source_path: &Path, input_format: &'static str) -> Self {
        Self {
            source_path: source_path.to_path_buf(),
            sample: None,
            input_format,
            record: None,
            line: None,
        }
    }
}

#[derive(Clone, Debug, Eq, PartialEq)]
pub(crate) struct LocatedRead {
    pub(crate) read: ReadRecord,
    pub(crate) location: ReadLocation,
    pub(crate) full_length_evidence: bool,
}

#[derive(Clone, Debug, Eq, PartialEq)]
pub(crate) struct ReadRejection {
    pub(crate) location: ReadLocation,
    pub(crate) read_id: Option<String>,
    pub(crate) stage: &'static str,
    pub(crate) reason: ReadRejectionReason,
    pub(crate) detail: String,
}

impl ReadRejection {
    pub(crate) fn new(
        location: ReadLocation,
        read_id: Option<String>,
        stage: &'static str,
        reason: ReadRejectionReason,
        detail: impl Into<String>,
    ) -> Self {
        Self {
            location,
            read_id,
            stage,
            reason,
            detail: detail.into(),
        }
    }

    pub(crate) fn context(&self) -> String {
        let mut context = self.location.source_path.display().to_string();
        if let Some(line) = self.location.line {
            context.push(':');
            context.push_str(&line.to_string());
        } else if let Some(record) = self.location.record {
            context.push_str(":record-");
            context.push_str(&record.to_string());
        }
        if let Some(read_id) = self.read_id.as_deref() {
            context.push_str(" read ");
            context.push_str(&format!("{read_id:?}"));
        }
        format!("{context}: {} ({})", self.detail, self.reason.as_str())
    }
}

#[derive(Clone, Debug, Default, Eq, PartialEq)]
pub(crate) struct ReadLoadOutcome {
    pub(crate) reads: Vec<LocatedRead>,
    pub(crate) rejections: Vec<ReadRejection>,
    pub(crate) total_records: u64,
}

impl ReadLoadOutcome {
    pub(crate) fn append(&mut self, mut other: Self) {
        self.reads.append(&mut other.reads);
        self.rejections.append(&mut other.rejections);
        self.total_records = self.total_records.saturating_add(other.total_records);
    }

    pub(crate) fn set_sample(&mut self, sample: &str) {
        for located in &mut self.reads {
            located.location.sample = Some(sample.to_owned());
        }
        for rejection in &mut self.rejections {
            rejection.location.sample = Some(sample.to_owned());
        }
    }
}

pub(crate) fn write_read_rejections(
    path: &Path,
    rejections: &[ReadRejection],
) -> anyhow::Result<()> {
    let metadata = vec![format!(
        "#trackclustertu_read_rejections_schema={READ_REJECTIONS_SCHEMA_VERSION}"
    )];
    let mut writer = DelimitedWriter::create(path, Delimiter::Tab, &metadata)?;
    writer.write_record([
        "sample",
        "source_path",
        "input_format",
        "record",
        "line",
        "read_id",
        "stage",
        "reason",
        "detail",
    ])?;
    for rejection in rejections {
        writer.write_record([
            rejection
                .location
                .sample
                .as_deref()
                .unwrap_or(".")
                .to_owned(),
            rejection.location.source_path.display().to_string(),
            rejection.location.input_format.to_owned(),
            rejection
                .location
                .record
                .map_or_else(|| ".".to_owned(), |value| value.to_string()),
            rejection
                .location
                .line
                .map_or_else(|| ".".to_owned(), |value| value.to_string()),
            rejection.read_id.as_deref().unwrap_or(".").to_owned(),
            rejection.stage.to_owned(),
            rejection.reason.as_str().to_owned(),
            rejection.detail.clone(),
        ])?;
    }
    writer.flush()?;
    Ok(())
}

pub(crate) fn emit_read_input_summary(
    total_records: u64,
    retained: usize,
    rejections: &[ReadRejection],
) {
    let mut counts: BTreeMap<&'static str, u64> = BTreeMap::new();
    for rejection in rejections {
        *counts.entry(rejection.reason.as_str()).or_default() += 1;
    }
    eprint!(
        "read_input_counts\ttotal={total_records}\tretained={retained}\trejected={}",
        rejections.len()
    );
    for (reason, count) in counts {
        eprint!("\t{reason}={count}");
    }
    eprintln!();
}
