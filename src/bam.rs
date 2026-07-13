//! Evidence-aware conversion of BAM alignments into boundary-preserving BED6 records.

use std::collections::{BTreeMap, HashMap};
use std::fmt;
use std::fs::File;
use std::io::{self, BufWriter, Write};
use std::path::{Path, PathBuf};
use std::str::FromStr;

use noodles::{
    bam as noodles_bam, bgzf as noodles_bgzf,
    sam::{
        self,
        alignment::{record::cigar::op::Kind as CigarKind, RecordBuf},
        header::record::value::map::header::{sort_order, tag as header_tag},
    },
};
use sha2::{Digest, Sha256};
use thiserror::Error;

use crate::io::bed::{self, Bed6Record};
use crate::io::delimited::{DelimitedError, DelimitedReader, DelimitedWriter, Delimiter};
use crate::model::{Coord, Strand};

/// Schema identifier written in every BAM evidence sidecar row.
pub const BAM_EVIDENCE_SCHEMA_VERSION: &str = "trackclustertu.bam-evidence.v1";

/// Typed failure returned by BAM-to-BED conversion APIs.
///
/// I/O variants retain their original [`std::io::Error`] as a source. Logical
/// validation variants carry the input path and the most specific available
/// line or BAM-record position so library callers do not need to parse display
/// strings to handle failures.
#[derive(Debug, Error)]
pub enum BamConversionError {
    /// Exact-boundary support was configured as zero.
    #[error("minimum boundary support must be at least 1, got {value}")]
    InvalidMinimumBoundarySupport {
        /// Rejected support value.
        value: usize,
    },

    /// The optional input evidence table could not be opened.
    #[error("failed to open input evidence TSV {path:?}: {source}")]
    EvidenceOpen {
        /// Evidence-table path.
        path: PathBuf,
        /// Underlying filesystem error.
        source: std::io::Error,
    },

    /// A physical evidence-table row could not be read.
    #[error("failed to read input evidence TSV {path:?}:{line}: {source}")]
    EvidenceRead {
        /// Evidence-table path.
        path: PathBuf,
        /// One-based source line number.
        line: usize,
        /// Underlying read error.
        source: std::io::Error,
    },

    /// The evidence-table header omits the required `read_name` column.
    #[error("input evidence TSV {path:?}:{line} is missing required read_name column")]
    EvidenceMissingReadNameColumn {
        /// Evidence-table path.
        path: PathBuf,
        /// One-based header line number.
        line: usize,
    },

    /// The evidence-table header repeats a column name.
    #[error("input evidence TSV {path:?}:{line} contains duplicate column {column:?}")]
    EvidenceDuplicateColumn {
        /// Evidence-table path.
        path: PathBuf,
        /// One-based header line number.
        line: usize,
        /// Normalized duplicate column name.
        column: String,
    },

    /// An evidence row has no usable read name.
    #[error("input evidence TSV {path:?}:{line} has an empty read_name")]
    EvidenceEmptyReadName {
        /// Evidence-table path.
        path: PathBuf,
        /// One-based source line number.
        line: usize,
    },

    /// A boolean-like evidence cell contains an unsupported value.
    #[error(
        "invalid evidence value {value:?} for column {column:?} at {path:?}:{line}; expected present/absent/unknown"
    )]
    EvidenceInvalidValue {
        /// Evidence-table path.
        path: PathBuf,
        /// One-based source line number.
        line: usize,
        /// Evidence column name.
        column: &'static str,
        /// Rejected cell text.
        value: String,
    },

    /// More than one evidence row uses the same read name.
    #[error("input evidence TSV {path:?}:{line} contains duplicate read_name {read_name:?}")]
    EvidenceDuplicateReadName {
        /// Evidence-table path.
        path: PathBuf,
        /// One-based source line number of the duplicate.
        line: usize,
        /// Repeated read name.
        read_name: String,
    },

    /// The evidence table contains no header or data row.
    #[error("input evidence TSV {path:?} is empty")]
    EvidenceEmpty {
        /// Evidence-table path.
        path: PathBuf,
    },

    /// A BAM file could not be opened.
    #[error("failed to open BAM {path:?} while {context}: {source}")]
    BamOpen {
        /// BAM path.
        path: PathBuf,
        /// Conversion phase that attempted the open.
        context: &'static str,
        /// Underlying filesystem error.
        source: std::io::Error,
    },

    /// A BAM header could not be decoded.
    #[error("failed to read BAM header from {path:?} while {context}: {source}")]
    BamHeaderRead {
        /// BAM path.
        path: PathBuf,
        /// Conversion phase reading the header.
        context: &'static str,
        /// Underlying BAM/header error.
        source: std::io::Error,
    },

    /// A BAM record could not be decoded.
    #[error("failed to read BAM record {record_ordinal} from {path:?} while {context}: {source}")]
    BamRecordRead {
        /// BAM path.
        path: PathBuf,
        /// One-based ordinal of the record being decoded.
        record_ordinal: usize,
        /// Conversion phase reading the record.
        context: &'static str,
        /// Underlying BAM record error.
        source: std::io::Error,
    },

    /// A mapped record appears after the unmapped tail in a declared coordinate-sorted BAM.
    #[error(
        "BAM {path:?} declares @HD SO:coordinate but mapped record {record_ordinal} follows the unmapped tail"
    )]
    MappedAfterUnmappedTail {
        /// BAM path.
        path: PathBuf,
        /// One-based offending record ordinal.
        record_ordinal: usize,
    },

    /// Coordinates decrease in a BAM declared as coordinate sorted.
    #[error(
        "BAM {path:?} declares @HD SO:coordinate but record {record_ordinal} is out of order: ({current_reference_id}, {current_position}) follows ({previous_reference_id}, {previous_position})"
    )]
    OutOfOrder {
        /// BAM path.
        path: PathBuf,
        /// One-based offending record ordinal.
        record_ordinal: usize,
        /// Current zero-based reference-sequence index.
        current_reference_id: usize,
        /// Current one-based alignment position.
        current_position: usize,
        /// Previous zero-based reference-sequence index.
        previous_reference_id: usize,
        /// Previous one-based alignment position.
        previous_position: usize,
    },

    /// A retained boundary appears only during the second support-counting pass.
    #[error(
        "BAM {path:?} changed between boundary-support passes: retained record {record_ordinal} ({read_name:?}) has a new boundary"
    )]
    BoundaryChangedBetweenPasses {
        /// BAM path.
        path: PathBuf,
        /// One-based record ordinal in the second pass.
        record_ordinal: usize,
        /// Read name reported by the second pass.
        read_name: String,
    },

    /// The number of records differs between boundary-support passes.
    #[error(
        "BAM {path:?} changed between boundary-support passes: first pass read {first_pass_records} records, second pass read {second_pass_records}"
    )]
    RecordCountChangedBetweenPasses {
        /// BAM path.
        path: PathBuf,
        /// Records observed in the first pass.
        first_pass_records: u64,
        /// Records observed in the second pass.
        second_pass_records: u64,
    },

    /// The number of recoverably malformed records differs between support-counting passes.
    #[error(
        "BAM {path:?} changed between boundary-support passes: first pass skipped {first_pass_records} malformed records, second pass skipped {second_pass_records}"
    )]
    MalformedRecordCountChangedBetweenPasses {
        /// BAM path.
        path: PathBuf,
        /// Recoverably malformed records in the first pass.
        first_pass_records: u64,
        /// Recoverably malformed records in the second pass.
        second_pass_records: u64,
    },

    /// A BAM header could not be encoded for the two-pass consistency fingerprint.
    #[error("failed to fingerprint BAM header from {path:?} while {context}: {source}")]
    BamHeaderFingerprint {
        /// BAM path.
        path: PathBuf,
        /// Conversion phase computing the fingerprint.
        context: &'static str,
        /// Underlying SAM header encoding error.
        source: std::io::Error,
    },

    /// Critical BAM content differs between boundary-support passes.
    #[error(
        "BAM {path:?} changed between boundary-support passes: first content SHA-256 {first_digest}, second content SHA-256 {second_digest}"
    )]
    ContentChangedBetweenPasses {
        /// BAM path.
        path: PathBuf,
        /// Stable content digest from the first pass.
        first_digest: String,
        /// Stable content digest from the second pass.
        second_digest: String,
    },

    /// A BAM record's CIGAR could not be rendered for the evidence table.
    #[error(
        "failed to encode CIGAR for BAM record {record_ordinal} ({read_name:?}) from {path:?}: {source}"
    )]
    CigarWrite {
        /// BAM path.
        path: PathBuf,
        /// One-based BAM record ordinal.
        record_ordinal: usize,
        /// Read name, or `.` when absent.
        read_name: String,
        /// Underlying CIGAR writer error.
        source: std::io::Error,
    },

    /// The SAM CIGAR writer unexpectedly produced non-UTF-8 bytes.
    #[error(
        "CIGAR for BAM record {record_ordinal} ({read_name:?}) from {path:?} is not UTF-8: {source}"
    )]
    CigarUtf8 {
        /// BAM path.
        path: PathBuf,
        /// One-based BAM record ordinal.
        record_ordinal: usize,
        /// Read name, or `.` when absent.
        read_name: String,
        /// UTF-8 decoding failure.
        source: std::string::FromUtf8Error,
    },

    /// A BED or evidence output could not be created.
    #[error("failed to create {output} output {path:?}: {source}")]
    OutputCreate {
        /// Output path.
        path: PathBuf,
        /// Human-readable output kind (`BED` or `BAM evidence TSV`).
        output: &'static str,
        /// Underlying filesystem error.
        source: std::io::Error,
    },

    /// A BED or evidence output could not be written or flushed.
    #[error("failed to write {output} output {path:?} while {context}: {source}")]
    OutputWrite {
        /// Output path.
        path: PathBuf,
        /// Human-readable output kind (`BED` or `BAM evidence TSV`).
        output: &'static str,
        /// Record/header/flush operation that failed.
        context: String,
        /// Underlying writer error.
        source: std::io::Error,
    },
}

type BamResult<T> = std::result::Result<T, BamConversionError>;

/// Library chemistry used to interpret optional full-length evidence.
///
/// [`FromStr`] accepts the displayed kebab-case names plus underscore and short
/// compatibility spellings (`drna`, `dcdna`, and `pcrcdna`). [`fmt::Display`]
/// always emits the canonical kebab-case name.
#[derive(Clone, Copy, Debug, Eq, Hash, Ord, PartialEq, PartialOrd)]
pub enum LibraryProfile {
    /// Oxford Nanopore direct-RNA sequencing.
    DirectRna,
    /// Oxford Nanopore direct-cDNA sequencing.
    DirectCdna,
    /// PCR-amplified cDNA sequencing.
    PcrCdna,
}

impl Default for LibraryProfile {
    fn default() -> Self {
        Self::DirectRna
    }
}

impl LibraryProfile {
    /// Returns the evidence rule used when no explicit `full_length_evidence` value is supplied.
    ///
    /// These rules are intentionally conservative. Direct RNA requires both poly(A) and 5-prime
    /// adapter evidence, while cDNA profiles require evidence at both ends. The default conversion
    /// does not require this rule to pass; it is activated by `require_full_length`.
    ///
    /// # Examples
    ///
    /// ```
    /// use trackclustertu::bam::LibraryProfile;
    ///
    /// assert_eq!(
    ///     LibraryProfile::DirectRna.inferred_full_length_rule(),
    ///     "poly_a_and_five_prime_adapter",
    /// );
    /// assert_eq!(
    ///     "pcr-cdna".parse::<LibraryProfile>()?,
    ///     LibraryProfile::PcrCdna,
    /// );
    /// # Ok::<(), Box<dyn std::error::Error>>(())
    /// ```
    pub const fn inferred_full_length_rule(self) -> &'static str {
        match self {
            Self::DirectRna => "poly_a_and_five_prime_adapter",
            Self::DirectCdna | Self::PcrCdna => "five_prime_adapter_and_three_prime_adapter",
        }
    }
}

impl fmt::Display for LibraryProfile {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str(match self {
            Self::DirectRna => "direct-rna",
            Self::DirectCdna => "direct-cdna",
            Self::PcrCdna => "pcr-cdna",
        })
    }
}

#[derive(Clone, Debug, Eq, PartialEq)]
/// Error returned when a library-profile name is not recognized.
pub struct ParseLibraryProfileError {
    value: String,
}

impl fmt::Display for ParseLibraryProfileError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(
            f,
            "unknown library profile {:?}; expected direct-rna, direct-cdna, or pcr-cdna",
            self.value
        )
    }
}

impl std::error::Error for ParseLibraryProfileError {}

impl FromStr for LibraryProfile {
    type Err = ParseLibraryProfileError;

    fn from_str(value: &str) -> std::result::Result<Self, Self::Err> {
        match value.trim().to_ascii_lowercase().as_str() {
            "direct-rna" | "direct_rna" | "drna" => Ok(Self::DirectRna),
            "direct-cdna" | "direct_cdna" | "dcdna" => Ok(Self::DirectCdna),
            "pcr-cdna" | "pcr_cdna" | "pcrcdna" => Ok(Self::PcrCdna),
            _ => Err(ParseLibraryProfileError {
                value: value.to_owned(),
            }),
        }
    }
}

/// Filters and evidence interpretation used during BAM conversion.
///
/// The default retains every otherwise-eligible primary unspliced alignment:
/// minimum MAPQ zero, no full-length requirement, exact-boundary support one,
/// and the [`LibraryProfile::DirectRna`] evidence rule.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct BamConversionConfig {
    /// Minimum mapping quality for a record to reach BED output.
    pub min_mapq: u8,
    /// Require explicit or profile-inferred full-length evidence.
    pub require_full_length: bool,
    /// Minimum number of otherwise-retained reads with exactly matching boundaries and strand.
    pub min_boundary_support: usize,
    /// Chemistry/library rule used for full-length inference and recorded in the sidecar.
    pub library_profile: LibraryProfile,
}

impl Default for BamConversionConfig {
    fn default() -> Self {
        Self {
            min_mapq: 0,
            require_full_length: false,
            min_boundary_support: 1,
            library_profile: LibraryProfile::default(),
        }
    }
}

/// Deterministic BAM conversion counters. Category counters are not mutually exclusive; e.g. a
/// supplementary record containing soft clipping contributes to both counters.
#[derive(Clone, Debug, Default, Eq, PartialEq)]
pub struct BamConversionSummary {
    /// BAM records inspected.
    pub total: u64,
    /// Primary records retained in BED output.
    pub retained: u64,
    /// Unmapped records.
    pub unmapped: u64,
    /// Secondary alignments.
    pub secondary: u64,
    /// Supplementary alignments.
    pub supplementary: u64,
    /// Records containing skipped/spliced reference operations.
    pub spliced_or_skipped: u64,
    /// Records with terminal soft clipping, regardless of final retention.
    pub soft_clipped: u64,
    /// Records without a usable read name.
    pub missing_name: u64,
    /// Complete BAM record blocks whose payload could not be decoded and were skipped.
    pub malformed_record: u64,
    /// Records below the configured mapping-quality threshold.
    pub below_min_mapq: u64,
    /// Records rejected for missing required full-length evidence.
    pub missing_full_length_evidence: u64,
    /// Otherwise eligible records below exact-boundary support.
    pub below_min_boundary_support: u64,
    /// Records missing a usable reference or alignment boundary.
    pub invalid_alignment: u64,
}

impl BamConversionSummary {
    /// Add every counter from `other` into this summary.
    pub fn add_assign(&mut self, other: &Self) {
        self.total += other.total;
        self.retained += other.retained;
        self.unmapped += other.unmapped;
        self.secondary += other.secondary;
        self.supplementary += other.supplementary;
        self.spliced_or_skipped += other.spliced_or_skipped;
        self.soft_clipped += other.soft_clipped;
        self.missing_name += other.missing_name;
        self.malformed_record += other.malformed_record;
        self.below_min_mapq += other.below_min_mapq;
        self.missing_full_length_evidence += other.missing_full_length_evidence;
        self.below_min_boundary_support += other.below_min_boundary_support;
        self.invalid_alignment += other.invalid_alignment;
    }

    /// Return all counters in the stable, machine-readable order used by the CLI.
    ///
    /// This includes every public counter exactly once and is suitable for
    /// deterministic text or structured reporting.
    pub fn fields(&self) -> [(&'static str, u64); 13] {
        [
            ("total", self.total),
            ("retained", self.retained),
            ("unmapped", self.unmapped),
            ("secondary", self.secondary),
            ("supplementary", self.supplementary),
            ("spliced_or_skipped", self.spliced_or_skipped),
            ("soft_clipped", self.soft_clipped),
            ("missing_name", self.missing_name),
            ("malformed_record", self.malformed_record),
            ("below_min_mapq", self.below_min_mapq),
            (
                "missing_full_length_evidence",
                self.missing_full_length_evidence,
            ),
            (
                "below_min_boundary_support",
                self.below_min_boundary_support,
            ),
            ("invalid_alignment", self.invalid_alignment),
        ]
    }
}

#[derive(Clone, Copy, Debug, Default, Eq, PartialEq)]
enum EvidenceState {
    Present,
    Absent,
    #[default]
    Unknown,
}

impl EvidenceState {
    fn parse(value: &str, path: &Path, line: usize, column: &'static str) -> BamResult<Self> {
        match value.trim().to_ascii_lowercase().as_str() {
            "1" | "true" | "yes" | "present" | "pass" => Ok(Self::Present),
            "0" | "false" | "no" | "absent" | "fail" => Ok(Self::Absent),
            "" | "." | "na" | "n/a" | "unknown" => Ok(Self::Unknown),
            _ => Err(BamConversionError::EvidenceInvalidValue {
                path: path.to_path_buf(),
                line,
                column,
                value: value.to_owned(),
            }),
        }
    }

    const fn is_present(self) -> bool {
        matches!(self, Self::Present)
    }
}

impl fmt::Display for EvidenceState {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str(match self {
            Self::Present => "present",
            Self::Absent => "absent",
            Self::Unknown => "unknown",
        })
    }
}

#[derive(Clone, Debug, Default, Eq, PartialEq)]
struct InputEvidence {
    poly_a: EvidenceState,
    five_prime_adapter: EvidenceState,
    three_prime_adapter: EvidenceState,
    full_length: EvidenceState,
    library_preparation: Option<String>,
}

impl InputEvidence {
    fn effective_full_length(&self, profile: LibraryProfile) -> bool {
        match self.full_length {
            EvidenceState::Present => true,
            EvidenceState::Absent => false,
            EvidenceState::Unknown => match profile {
                LibraryProfile::DirectRna => {
                    self.poly_a.is_present() && self.five_prime_adapter.is_present()
                }
                LibraryProfile::DirectCdna | LibraryProfile::PcrCdna => {
                    self.five_prime_adapter.is_present() && self.three_prime_adapter.is_present()
                }
            },
        }
    }
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
enum FilterReason {
    Retained,
    Unmapped,
    MalformedRecord,
    MissingName,
    Secondary,
    Supplementary,
    SplicedOrSkipped,
    MissingReference,
    MissingAlignmentStart,
    MissingAlignmentEnd,
    BelowMinMapq,
    MissingFullLengthEvidence,
    BelowMinBoundarySupport,
}

impl FilterReason {
    const fn as_str(self) -> &'static str {
        match self {
            Self::Retained => "retained",
            Self::Unmapped => "unmapped",
            Self::MalformedRecord => "malformed_record",
            Self::MissingName => "missing_name",
            Self::Secondary => "secondary",
            Self::Supplementary => "supplementary",
            Self::SplicedOrSkipped => "spliced_or_skipped",
            Self::MissingReference => "missing_reference",
            Self::MissingAlignmentStart => "missing_alignment_start",
            Self::MissingAlignmentEnd => "missing_alignment_end",
            Self::BelowMinMapq => "below_min_mapq",
            Self::MissingFullLengthEvidence => "missing_full_length_evidence",
            Self::BelowMinBoundarySupport => "below_min_boundary_support",
        }
    }
}

#[derive(Clone, Debug, Eq, Hash, PartialEq)]
struct BoundaryKey {
    chrom: String,
    start: Coord,
    end: Coord,
    strand: Strand,
}

#[derive(Debug)]
struct ConvertedRecord {
    ordinal: usize,
    name: String,
    chrom: Option<String>,
    start: Option<Coord>,
    end: Option<Coord>,
    strand: Strand,
    mapq: Option<u8>,
    five_prime_soft_clip: usize,
    three_prime_soft_clip: usize,
    cigar: String,
    evidence: InputEvidence,
    effective_full_length: bool,
    filter_reason: FilterReason,
    bed_record: Option<Bed6Record>,
    boundary_key: Option<BoundaryKey>,
}

impl ConvertedRecord {
    fn malformed(ordinal: usize) -> Self {
        Self {
            ordinal,
            name: ".".to_owned(),
            chrom: None,
            start: None,
            end: None,
            strand: Strand::Unknown,
            mapq: None,
            five_prime_soft_clip: 0,
            three_prime_soft_clip: 0,
            cigar: ".".to_owned(),
            evidence: InputEvidence::default(),
            effective_full_length: false,
            filter_reason: FilterReason::MalformedRecord,
            bed_record: None,
            boundary_key: None,
        }
    }
}

struct RecordConversionContext<'a> {
    bam_path: &'a Path,
    header: &'a sam::Header,
    input_evidence: &'a HashMap<String, InputEvidence>,
    config: &'a BamConversionConfig,
    capture_evidence: bool,
}

const BED_OUTPUT: &str = "BED";
const EVIDENCE_OUTPUT: &str = "BAM evidence TSV";

type FileBamReader = noodles_bam::io::Reader<noodles_bgzf::io::Reader<File>>;

enum BamRecordReadOutcome {
    End,
    Decoded,
    Malformed { message: String },
}

struct BamPassDigest {
    hasher: Sha256,
}

impl BamPassDigest {
    fn new(header: &sam::Header) -> io::Result<Self> {
        let mut header_writer = sam::io::Writer::new(Vec::new());
        header_writer.write_header(header)?;

        let mut digest = Self {
            hasher: Sha256::new(),
        };
        digest.update_field(b"digest_schema", b"trackclustertu.bam-pass-digest.v1");
        digest.update_field(b"header", &header_writer.into_inner());
        Ok(digest)
    }

    fn update_decoded(&mut self, raw: &RecordBuf, converted: &ConvertedRecord) {
        self.update_field(b"record_kind", b"decoded");
        self.update_usize(b"ordinal", converted.ordinal);
        self.update_optional_field(b"raw_name", raw.name().map(|name| name.as_ref()));
        self.update_optional_usize(b"raw_reference_id", raw.reference_sequence_id());
        self.update_optional_usize(
            b"raw_alignment_start",
            raw.alignment_start().map(usize::from),
        );
        self.update_u64(b"raw_flags", u64::from(u16::from(raw.flags())));
        self.update_optional_u64(
            b"raw_mapq",
            raw.mapping_quality()
                .map(|value| u64::from(u8::from(value))),
        );

        self.update_field(b"name", converted.name.as_bytes());
        self.update_optional_field(b"chrom", converted.chrom.as_deref().map(str::as_bytes));
        self.update_optional_u64(
            b"start",
            converted.start.map(|value| u64::from(value.get())),
        );
        self.update_optional_u64(b"end", converted.end.map(|value| u64::from(value.get())));
        self.update_u64(b"strand", u64::from(converted.strand.as_char() as u32));
        self.update_optional_u64(b"mapq", converted.mapq.map(u64::from));
        self.update_usize(b"five_prime_soft_clip", converted.five_prime_soft_clip);
        self.update_usize(b"three_prime_soft_clip", converted.three_prime_soft_clip);
        self.update_field(b"cigar", converted.cigar.as_bytes());
        self.update_u64(
            b"poly_a_evidence",
            u64::from(evidence_state_code(converted.evidence.poly_a)),
        );
        self.update_u64(
            b"five_prime_adapter_evidence",
            u64::from(evidence_state_code(converted.evidence.five_prime_adapter)),
        );
        self.update_u64(
            b"three_prime_adapter_evidence",
            u64::from(evidence_state_code(converted.evidence.three_prime_adapter)),
        );
        self.update_u64(
            b"full_length_evidence",
            u64::from(evidence_state_code(converted.evidence.full_length)),
        );
        self.update_optional_field(
            b"library_preparation",
            converted
                .evidence
                .library_preparation
                .as_deref()
                .map(str::as_bytes),
        );
        self.update_u64(
            b"effective_full_length",
            u64::from(converted.effective_full_length),
        );
        self.update_field(
            b"filter_reason",
            converted.filter_reason.as_str().as_bytes(),
        );

        self.update_u64(
            b"boundary_present",
            u64::from(converted.boundary_key.is_some()),
        );
        if let Some(boundary) = &converted.boundary_key {
            self.update_field(b"boundary_chrom", boundary.chrom.as_bytes());
            self.update_u64(b"boundary_start", u64::from(boundary.start.get()));
            self.update_u64(b"boundary_end", u64::from(boundary.end.get()));
            self.update_u64(
                b"boundary_strand",
                u64::from(boundary.strand.as_char() as u32),
            );
        }

        self.update_u64(b"bed_present", u64::from(converted.bed_record.is_some()));
        if let Some(bed) = &converted.bed_record {
            self.update_field(b"bed_chrom", bed.chrom.as_bytes());
            self.update_u64(b"bed_start", u64::from(bed.start.get()));
            self.update_u64(b"bed_end", u64::from(bed.end.get()));
            self.update_field(b"bed_name", bed.name.as_bytes());
            self.update_u64(b"bed_score", u64::from(bed.score));
            self.update_u64(b"bed_strand", u64::from(bed.strand.as_char() as u32));
            self.update_usize(b"bed_extra_field_count", bed.extra_fields.len());
            for field in &bed.extra_fields {
                self.update_field(b"bed_extra_field", field.as_bytes());
            }
        }
    }

    fn update_malformed(&mut self, ordinal: usize, message: &str) {
        self.update_field(b"record_kind", b"malformed");
        self.update_usize(b"ordinal", ordinal);
        self.update_field(b"decoder_error", message.as_bytes());
    }

    fn update_field(&mut self, label: &[u8], value: &[u8]) {
        self.hasher
            .update(u64::try_from(label.len()).unwrap_or(u64::MAX).to_le_bytes());
        self.hasher.update(label);
        self.hasher
            .update(u64::try_from(value.len()).unwrap_or(u64::MAX).to_le_bytes());
        self.hasher.update(value);
    }

    fn update_usize(&mut self, label: &[u8], value: usize) {
        self.update_u64(label, u64::try_from(value).unwrap_or(u64::MAX));
    }

    fn update_u64(&mut self, label: &[u8], value: u64) {
        self.update_field(label, &value.to_le_bytes());
    }

    fn update_optional_u64(&mut self, label: &[u8], value: Option<u64>) {
        self.update_u64(label, u64::from(value.is_some()));
        if let Some(value) = value {
            self.update_field(b"optional_u64_value", &value.to_le_bytes());
        }
    }

    fn update_optional_usize(&mut self, label: &[u8], value: Option<usize>) {
        self.update_optional_u64(
            label,
            value.map(|value| u64::try_from(value).unwrap_or(u64::MAX)),
        );
    }

    fn update_optional_field(&mut self, label: &[u8], value: Option<&[u8]>) {
        self.update_u64(label, u64::from(value.is_some()));
        if let Some(value) = value {
            self.update_field(b"optional_bytes_value", value);
        }
    }

    fn finish(self) -> [u8; 32] {
        self.hasher.finalize().into()
    }
}

const fn evidence_state_code(state: EvidenceState) -> u8 {
    match state {
        EvidenceState::Present => 1,
        EvidenceState::Absent => 2,
        EvidenceState::Unknown => 0,
    }
}

fn digest_hex(digest: &[u8; 32]) -> String {
    use std::fmt::Write as _;

    let mut output = String::with_capacity(digest.len() * 2);
    for byte in digest {
        write!(&mut output, "{byte:02x}").expect("writing to a String cannot fail");
    }
    output
}

fn ensure_pass_digests_match(
    path: &Path,
    first_digest: [u8; 32],
    second_digest: [u8; 32],
) -> BamResult<()> {
    if first_digest == second_digest {
        Ok(())
    } else {
        Err(BamConversionError::ContentChangedBetweenPasses {
            path: path.to_path_buf(),
            first_digest: digest_hex(&first_digest),
            second_digest: digest_hex(&second_digest),
        })
    }
}

fn open_bam(path: &Path, context: &'static str) -> BamResult<FileBamReader> {
    File::open(path)
        .map_err(|source| BamConversionError::BamOpen {
            path: path.to_path_buf(),
            context,
            source,
        })
        .map(noodles_bam::io::Reader::new)
}

fn read_bam_header(
    reader: &mut FileBamReader,
    path: &Path,
    context: &'static str,
) -> BamResult<sam::Header> {
    reader
        .read_header()
        .map_err(|source| BamConversionError::BamHeaderRead {
            path: path.to_path_buf(),
            context,
            source,
        })
}

fn read_bam_record(
    reader: &mut FileBamReader,
    header: &sam::Header,
    record: &mut RecordBuf,
    path: &Path,
    record_ordinal: usize,
    context: &'static str,
) -> BamResult<BamRecordReadOutcome> {
    match reader.read_record_buf(header, record) {
        Ok(0) => Ok(BamRecordReadOutcome::End),
        Ok(_) => Ok(BamRecordReadOutcome::Decoded),
        Err(source) => {
            let decode_message = if source.kind() == io::ErrorKind::InvalidData {
                source
                    .get_ref()
                    .and_then(|error| {
                        error.downcast_ref::<noodles_bam::record::codec::decoder::DecodeError>()
                    })
                    .map(ToString::to_string)
            } else {
                None
            };

            if let Some(message) = decode_message {
                // noodles reads the declared record block completely before decoding it, so a
                // decoder error leaves the stream at the next record. The reusable buffer may be
                // partially mutated and must not be reused as-is.
                *record = RecordBuf::default();
                Ok(BamRecordReadOutcome::Malformed { message })
            } else {
                Err(BamConversionError::BamRecordRead {
                    path: path.to_path_buf(),
                    record_ordinal,
                    context,
                    source,
                })
            }
        }
    }
}

fn report_malformed_record(path: &Path, ordinal: usize, message: &str) {
    eprintln!(
        "bam_record_skipped\tpath={path:?}\trecord_ordinal={ordinal}\treason=malformed_record\terror={message:?}"
    );
}

fn create_output(path: &Path, output: &'static str) -> BamResult<BufWriter<File>> {
    File::create(path)
        .map(BufWriter::new)
        .map_err(|source| BamConversionError::OutputCreate {
            path: path.to_path_buf(),
            output,
            source,
        })
}

fn output_write_error(
    path: &Path,
    output: &'static str,
    context: impl Into<String>,
    source: io::Error,
) -> BamConversionError {
    BamConversionError::OutputWrite {
        path: path.to_path_buf(),
        output,
        context: context.into(),
        source,
    }
}

fn flush_output(writer: &mut BufWriter<File>, path: &Path, output: &'static str) -> BamResult<()> {
    writer
        .flush()
        .map_err(|source| output_write_error(path, output, "flushing", source))
}

fn write_bed_record(
    writer: &mut BufWriter<File>,
    path: &Path,
    record: &Bed6Record,
    context: impl Into<String>,
) -> BamResult<()> {
    bed::write_bed6_to_writer(writer, [record])
        .map_err(|source| output_write_error(path, BED_OUTPUT, context, source))
}

/// Convert primary, unspliced BAM alignments to BED6 with backward-compatible defaults.
///
/// This is equivalent to calling [`bam_to_bed6_with_evidence`] with no evidence
/// input/output and [`BamConversionConfig::default()`]. Coordinate-sorted BAMs
/// stream directly; BAMs without a coordinate-sort declaration use the
/// documented in-memory fallback.
///
/// # Errors
///
/// Returns an error if the BAM/header cannot be read, an eligible alignment has
/// invalid coordinates, or the BED output cannot be created or written.
pub fn bam_to_bed6(
    bam_path: &Path,
    out_bed_path: &Path,
) -> std::result::Result<(), BamConversionError> {
    bam_to_bed6_with_evidence(
        bam_path,
        out_bed_path,
        None,
        None,
        &BamConversionConfig::default(),
    )?;
    Ok(())
}

/// Convert BAM alignments to BED6 and optionally write a versioned per-record evidence TSV.
///
/// `input_evidence_path` is a tab-separated file keyed by `read_name`. It may contain
/// `poly_a_evidence`, `five_prime_adapter_evidence`, `three_prime_adapter_evidence`,
/// `full_length_evidence`, and `library_preparation`. Boolean evidence accepts
/// present/absent/unknown and common true/false spellings. The evidence output contains every BAM
/// record, including records filtered from BED6, in BAM order.
///
/// BAMs declared `SO:coordinate` are processed with memory independent of BAM record count. A
/// minimum-boundary-support filter uses two passes, stores only distinct boundary keys, and checks
/// a stable SHA-256 digest of the header and every critical per-record output field before
/// publishing second-pass output. Input evidence is keyed by read name and therefore uses memory
/// proportional to input evidence rows. BAMs without `SO:coordinate` retain the documented
/// in-memory sort fallback.
/// A length-delimited record block that was read completely but fails noodles' payload decoder is
/// skipped, counted as `malformed_record`, and represented by a placeholder evidence row. BAM
/// framing, BGZF, truncation, and I/O failures remain fatal because the next record boundary is not
/// known to be recoverable.
///
/// `out_evidence_path` receives one audit row for every BAM record, including
/// filtered records. `input_evidence_path` may supply optional evidence keyed by
/// read name. The returned [`BamConversionSummary`] reports deterministic filter
/// and retention counters.
///
/// # Errors
///
/// Returns an error for invalid configuration (including zero boundary support), malformed
/// evidence input, structurally unreadable BAM data, false coordinate-sort declarations,
/// inconsistent two-pass input, or output I/O failures. Complete length-delimited record blocks
/// with malformed payloads are skipped instead of returned as errors.
pub fn bam_to_bed6_with_evidence(
    bam_path: &Path,
    out_bed_path: &Path,
    out_evidence_path: Option<&Path>,
    input_evidence_path: Option<&Path>,
    config: &BamConversionConfig,
) -> std::result::Result<BamConversionSummary, BamConversionError> {
    if config.min_boundary_support == 0 {
        return Err(BamConversionError::InvalidMinimumBoundarySupport {
            value: config.min_boundary_support,
        });
    }

    let input_evidence = match input_evidence_path {
        Some(path) => read_input_evidence(path)?,
        None => HashMap::new(),
    };

    let mut header_reader = open_bam(bam_path, "inspecting sort order")?;
    let header = read_bam_header(&mut header_reader, bam_path, "inspecting sort order")?;
    let coordinate_sorted = header
        .header()
        .and_then(|header| header.other_fields().get(&header_tag::SORT_ORDER))
        .is_some_and(|value| value.as_slice() == sort_order::COORDINATE);
    drop(header_reader);

    if coordinate_sorted && config.min_boundary_support > 1 {
        convert_coordinate_sorted_two_pass(
            bam_path,
            out_bed_path,
            out_evidence_path,
            &input_evidence,
            config,
        )
    } else if coordinate_sorted {
        convert_coordinate_sorted_one_pass(
            bam_path,
            out_bed_path,
            out_evidence_path,
            &input_evidence,
            config,
        )
    } else {
        eprintln!(
            "warning: BAM {:?} does not declare @HD SO:coordinate; falling back to in-memory record buffering and BED sorting; coordinate-sort the BAM for record-count-independent memory use",
            bam_path
        );
        convert_unsorted_buffered(
            bam_path,
            out_bed_path,
            out_evidence_path,
            &input_evidence,
            config,
        )
    }
}

fn convert_coordinate_sorted_one_pass(
    bam_path: &Path,
    out_bed_path: &Path,
    out_evidence_path: Option<&Path>,
    input_evidence: &HashMap<String, InputEvidence>,
    config: &BamConversionConfig,
) -> BamResult<BamConversionSummary> {
    const CONTEXT: &str = "streaming coordinate-sorted conversion";
    let mut reader = open_bam(bam_path, CONTEXT)?;
    let header = read_bam_header(&mut reader, bam_path, CONTEXT)?;
    let mut bed_writer = create_output(out_bed_path, BED_OUTPUT)?;
    let mut evidence_writer = out_evidence_path.map(create_evidence_writer).transpose()?;
    let mut summary = BamConversionSummary::default();
    let mut previous_coordinate: Option<(usize, usize)> = None;
    let mut observed_unmapped_tail = false;
    let mut record = RecordBuf::default();
    let conversion = RecordConversionContext {
        bam_path,
        header: &header,
        input_evidence,
        config,
        capture_evidence: evidence_writer.is_some(),
    };
    loop {
        let ordinal = usize::try_from(summary.total.saturating_add(1)).unwrap_or(usize::MAX);
        let read_outcome = read_bam_record(
            &mut reader,
            &header,
            &mut record,
            bam_path,
            ordinal,
            CONTEXT,
        )?;
        let converted = match read_outcome {
            BamRecordReadOutcome::End => break,
            BamRecordReadOutcome::Decoded => {
                summary.total += 1;
                validate_next_coordinate(
                    bam_path,
                    ordinal,
                    &record,
                    &mut previous_coordinate,
                    &mut observed_unmapped_tail,
                )?;
                convert_record(&conversion, ordinal, &record, &mut summary)?
            }
            BamRecordReadOutcome::Malformed { message } => {
                summary.total += 1;
                summary.malformed_record += 1;
                report_malformed_record(bam_path, ordinal, &message);
                ConvertedRecord::malformed(ordinal)
            }
        };
        if converted.filter_reason == FilterReason::Retained {
            summary.retained += 1;
            if let Some(bed_record) = converted.bed_record.as_ref() {
                write_bed_record(
                    &mut bed_writer,
                    out_bed_path,
                    bed_record,
                    format!("writing BAM record {ordinal}"),
                )?;
            }
        }
        if let (Some(path), Some(writer)) = (out_evidence_path, evidence_writer.as_mut()) {
            write_evidence_record(writer, &converted, config).map_err(|source| {
                output_write_error(
                    path,
                    EVIDENCE_OUTPUT,
                    format!("writing BAM record {ordinal}"),
                    io::Error::other(source),
                )
            })?;
        }
    }
    flush_output(&mut bed_writer, out_bed_path, BED_OUTPUT)?;
    if let Some(mut writer) = evidence_writer {
        flush_evidence_writer(
            &mut writer,
            out_evidence_path.expect("writer exists only when path is present"),
        )?;
    }
    Ok(summary)
}

fn convert_coordinate_sorted_two_pass(
    bam_path: &Path,
    out_bed_path: &Path,
    out_evidence_path: Option<&Path>,
    input_evidence: &HashMap<String, InputEvidence>,
    config: &BamConversionConfig,
) -> BamResult<BamConversionSummary> {
    // Pass 1 validates actual coordinate order, records deterministic counters once, and stores
    // only distinct boundary keys. Memory is O(distinct boundaries + input evidence rows), not
    // O(BAM records).
    const PASS_ONE_CONTEXT: &str = "reading boundary-support pass 1";
    let mut reader = open_bam(bam_path, PASS_ONE_CONTEXT)?;
    let header = read_bam_header(&mut reader, bam_path, PASS_ONE_CONTEXT)?;
    let mut pass_one_digest =
        BamPassDigest::new(&header).map_err(|source| BamConversionError::BamHeaderFingerprint {
            path: bam_path.to_path_buf(),
            context: PASS_ONE_CONTEXT,
            source,
        })?;
    let mut summary = BamConversionSummary::default();
    let mut support: HashMap<BoundaryKey, usize> = HashMap::new();
    let mut previous_coordinate: Option<(usize, usize)> = None;
    let mut observed_unmapped_tail = false;
    let mut record = RecordBuf::default();
    let pass_one_conversion = RecordConversionContext {
        bam_path,
        header: &header,
        input_evidence,
        config,
        capture_evidence: true,
    };
    loop {
        let ordinal = usize::try_from(summary.total.saturating_add(1)).unwrap_or(usize::MAX);
        match read_bam_record(
            &mut reader,
            &header,
            &mut record,
            bam_path,
            ordinal,
            PASS_ONE_CONTEXT,
        )? {
            BamRecordReadOutcome::End => break,
            BamRecordReadOutcome::Malformed { message } => {
                summary.total += 1;
                summary.malformed_record += 1;
                pass_one_digest.update_malformed(ordinal, &message);
                report_malformed_record(bam_path, ordinal, &message);
                continue;
            }
            BamRecordReadOutcome::Decoded => {}
        }
        summary.total += 1;
        validate_next_coordinate(
            bam_path,
            ordinal,
            &record,
            &mut previous_coordinate,
            &mut observed_unmapped_tail,
        )?;
        let converted = convert_record(&pass_one_conversion, ordinal, &record, &mut summary)?;
        pass_one_digest.update_decoded(&record, &converted);
        if converted.filter_reason == FilterReason::Retained {
            if let Some(key) = converted.boundary_key {
                let count = support.entry(key).or_default();
                *count = count.saturating_add(1);
            }
        }
    }
    let pass_one_digest = pass_one_digest.finish();
    drop(reader);

    // Pass 2 reconstructs one record at a time, applies the known support, and streams BED and
    // evidence in BAM order. Its conversion counters are deliberately discarded so the public
    // summary describes records, not passes.
    const PASS_TWO_CONTEXT: &str = "reading boundary-support pass 2";
    let mut reader = open_bam(bam_path, PASS_TWO_CONTEXT)?;
    let header = read_bam_header(&mut reader, bam_path, PASS_TWO_CONTEXT)?;
    let mut pass_two_digest =
        BamPassDigest::new(&header).map_err(|source| BamConversionError::BamHeaderFingerprint {
            path: bam_path.to_path_buf(),
            context: PASS_TWO_CONTEXT,
            source,
        })?;
    let mut bed_writer = create_output(out_bed_path, BED_OUTPUT)?;
    let mut evidence_writer = out_evidence_path.map(create_evidence_writer).transpose()?;
    let mut discarded_summary = BamConversionSummary::default();
    let mut pass_two_total = 0u64;
    let mut previous_coordinate: Option<(usize, usize)> = None;
    let mut observed_unmapped_tail = false;
    let mut record = RecordBuf::default();
    let pass_two_conversion = RecordConversionContext {
        bam_path,
        header: &header,
        input_evidence,
        config,
        capture_evidence: true,
    };
    loop {
        let ordinal = usize::try_from(pass_two_total.saturating_add(1)).unwrap_or(usize::MAX);
        let read_outcome = read_bam_record(
            &mut reader,
            &header,
            &mut record,
            bam_path,
            ordinal,
            PASS_TWO_CONTEXT,
        )?;
        let mut converted = match read_outcome {
            BamRecordReadOutcome::End => break,
            BamRecordReadOutcome::Malformed { message } => {
                pass_two_total += 1;
                discarded_summary.malformed_record += 1;
                pass_two_digest.update_malformed(ordinal, &message);
                ConvertedRecord::malformed(ordinal)
            }
            BamRecordReadOutcome::Decoded => {
                pass_two_total += 1;
                validate_next_coordinate(
                    bam_path,
                    ordinal,
                    &record,
                    &mut previous_coordinate,
                    &mut observed_unmapped_tail,
                )?;
                let converted = convert_record(
                    &pass_two_conversion,
                    ordinal,
                    &record,
                    &mut discarded_summary,
                )?;
                pass_two_digest.update_decoded(&record, &converted);
                converted
            }
        };
        let boundary_support = converted
            .boundary_key
            .as_ref()
            .and_then(|key| support.get(key).copied());
        if converted.filter_reason == FilterReason::Retained && boundary_support.is_none() {
            return Err(BamConversionError::BoundaryChangedBetweenPasses {
                path: bam_path.to_path_buf(),
                record_ordinal: ordinal,
                read_name: converted.name,
            });
        }
        let below_support = converted.filter_reason == FilterReason::Retained
            && boundary_support.is_some_and(|count| count < config.min_boundary_support);
        if below_support {
            converted.filter_reason = FilterReason::BelowMinBoundarySupport;
            converted.bed_record = None;
            summary.below_min_boundary_support += 1;
        } else if converted.filter_reason == FilterReason::Retained {
            summary.retained += 1;
            if let Some(bed_record) = converted.bed_record.as_ref() {
                write_bed_record(
                    &mut bed_writer,
                    out_bed_path,
                    bed_record,
                    format!("writing BAM record {ordinal}"),
                )?;
            }
        }
        if let (Some(path), Some(writer)) = (out_evidence_path, evidence_writer.as_mut()) {
            write_evidence_record(writer, &converted, config).map_err(|source| {
                output_write_error(
                    path,
                    EVIDENCE_OUTPUT,
                    format!("writing BAM record {ordinal}"),
                    io::Error::other(source),
                )
            })?;
        }
    }
    let pass_two_digest = pass_two_digest.finish();
    if pass_two_total != summary.total {
        return Err(BamConversionError::RecordCountChangedBetweenPasses {
            path: bam_path.to_path_buf(),
            first_pass_records: summary.total,
            second_pass_records: pass_two_total,
        });
    }
    if discarded_summary.malformed_record != summary.malformed_record {
        return Err(
            BamConversionError::MalformedRecordCountChangedBetweenPasses {
                path: bam_path.to_path_buf(),
                first_pass_records: summary.malformed_record,
                second_pass_records: discarded_summary.malformed_record,
            },
        );
    }
    ensure_pass_digests_match(bam_path, pass_one_digest, pass_two_digest)?;
    flush_output(&mut bed_writer, out_bed_path, BED_OUTPUT)?;
    if let Some(mut writer) = evidence_writer {
        flush_evidence_writer(
            &mut writer,
            out_evidence_path.expect("writer exists only when path is present"),
        )?;
    }
    Ok(summary)
}

fn convert_unsorted_buffered(
    bam_path: &Path,
    out_bed_path: &Path,
    out_evidence_path: Option<&Path>,
    input_evidence: &HashMap<String, InputEvidence>,
    config: &BamConversionConfig,
) -> BamResult<BamConversionSummary> {
    const CONTEXT: &str = "buffering a BAM without a coordinate-sort declaration";
    let mut reader = open_bam(bam_path, CONTEXT)?;
    let header = read_bam_header(&mut reader, bam_path, CONTEXT)?;
    let capture_evidence = out_evidence_path.is_some();
    let mut records = Vec::new();
    let mut summary = BamConversionSummary::default();
    let mut record = RecordBuf::default();
    let conversion = RecordConversionContext {
        bam_path,
        header: &header,
        input_evidence,
        config,
        capture_evidence,
    };
    loop {
        let ordinal = usize::try_from(summary.total.saturating_add(1)).unwrap_or(usize::MAX);
        match read_bam_record(
            &mut reader,
            &header,
            &mut record,
            bam_path,
            ordinal,
            CONTEXT,
        )? {
            BamRecordReadOutcome::End => break,
            BamRecordReadOutcome::Malformed { message } => {
                summary.total += 1;
                summary.malformed_record += 1;
                report_malformed_record(bam_path, ordinal, &message);
                records.push(ConvertedRecord::malformed(ordinal));
            }
            BamRecordReadOutcome::Decoded => {
                summary.total += 1;
                records.push(convert_record(&conversion, ordinal, &record, &mut summary)?);
            }
        }
    }

    if config.min_boundary_support > 1 {
        apply_min_boundary_support(&mut records, config.min_boundary_support, &mut summary);
    }

    let mut bed_records: Vec<Bed6Record> = records
        .iter_mut()
        .filter_map(|record| {
            if record.filter_reason == FilterReason::Retained {
                summary.retained += 1;
                record.bed_record.take()
            } else {
                None
            }
        })
        .collect();
    bed_records.sort_by(|left, right| {
        left.chrom
            .cmp(&right.chrom)
            .then_with(|| left.start.cmp(&right.start))
            .then_with(|| left.end.cmp(&right.end))
            .then_with(|| left.strand.cmp(&right.strand))
            .then_with(|| left.name.cmp(&right.name))
    });
    let mut bed_writer = create_output(out_bed_path, BED_OUTPUT)?;
    for (index, bed_record) in bed_records.iter().enumerate() {
        write_bed_record(
            &mut bed_writer,
            out_bed_path,
            bed_record,
            format!("writing sorted BED row {}", index + 1),
        )?;
    }
    flush_output(&mut bed_writer, out_bed_path, BED_OUTPUT)?;

    if let Some(path) = out_evidence_path {
        write_evidence_sidecar(path, &records, config)?;
    }

    Ok(summary)
}

fn validate_next_coordinate(
    bam_path: &Path,
    ordinal: usize,
    record: &RecordBuf,
    previous: &mut Option<(usize, usize)>,
    observed_unmapped_tail: &mut bool,
) -> BamResult<()> {
    if record.flags().is_unmapped() {
        *observed_unmapped_tail = true;
        return Ok(());
    }

    let (Some(reference_id), Some(start)) =
        (record.reference_sequence_id(), record.alignment_start())
    else {
        return Ok(());
    };
    if *observed_unmapped_tail {
        return Err(BamConversionError::MappedAfterUnmappedTail {
            path: bam_path.to_path_buf(),
            record_ordinal: ordinal,
        });
    }

    let coordinate = (reference_id, usize::from(start));
    if let Some(&(previous_reference_id, previous_position)) = previous.as_ref() {
        if coordinate < (previous_reference_id, previous_position) {
            return Err(BamConversionError::OutOfOrder {
                path: bam_path.to_path_buf(),
                record_ordinal: ordinal,
                current_reference_id: coordinate.0,
                current_position: coordinate.1,
                previous_reference_id,
                previous_position,
            });
        }
    }
    *previous = Some(coordinate);
    Ok(())
}

fn convert_record(
    context: &RecordConversionContext<'_>,
    ordinal: usize,
    record: &RecordBuf,
    summary: &mut BamConversionSummary,
) -> BamResult<ConvertedRecord> {
    let RecordConversionContext {
        bam_path,
        header,
        input_evidence,
        config,
        capture_evidence,
    } = context;
    let flags = record.flags();
    if flags.is_unmapped() {
        summary.unmapped += 1;
    }
    if flags.is_secondary() {
        summary.secondary += 1;
    }
    if flags.is_supplementary() {
        summary.supplementary += 1;
    }

    let cigar_ops = record.cigar().as_ref();
    let is_spliced = cigar_ops.iter().any(|op| op.kind() == CigarKind::Skip);
    if is_spliced {
        summary.spliced_or_skipped += 1;
    }
    let (leading_soft_clip, trailing_soft_clip) = terminal_soft_clips(cigar_ops);
    if leading_soft_clip > 0 || trailing_soft_clip > 0 {
        summary.soft_clipped += 1;
    }

    let strand = if flags.is_reverse_complemented() {
        Strand::Minus
    } else {
        Strand::Plus
    };
    let (five_prime_soft_clip, three_prime_soft_clip) = match strand {
        Strand::Plus | Strand::Unknown => (leading_soft_clip, trailing_soft_clip),
        Strand::Minus => (trailing_soft_clip, leading_soft_clip),
    };

    let record_name = record.name().map(|name| name.to_string());
    let missing_name = record_name.is_none();
    if missing_name {
        summary.missing_name += 1;
    }
    let name = record_name.unwrap_or_else(|| ".".to_owned());
    let evidence = if missing_name {
        InputEvidence::default()
    } else {
        input_evidence.get(&name).cloned().unwrap_or_default()
    };
    let effective_full_length = evidence.effective_full_length(config.library_profile);
    let mapq = record.mapping_quality().map(u8::from);

    let chrom = record
        .reference_sequence_id()
        .and_then(|id| header.reference_sequences().get_index(id))
        .map(|(name, _)| name.to_string());
    let start = record
        .alignment_start()
        .map(|position| Coord::new((usize::from(position) - 1) as u32));
    let end = record
        .alignment_end()
        .map(|position| Coord::new(usize::from(position) as u32));

    let filter_reason = if flags.is_unmapped() {
        FilterReason::Unmapped
    } else if missing_name {
        FilterReason::MissingName
    } else if flags.is_secondary() {
        FilterReason::Secondary
    } else if flags.is_supplementary() {
        FilterReason::Supplementary
    } else if is_spliced {
        FilterReason::SplicedOrSkipped
    } else if chrom.is_none() {
        summary.invalid_alignment += 1;
        FilterReason::MissingReference
    } else if start.is_none() {
        summary.invalid_alignment += 1;
        FilterReason::MissingAlignmentStart
    } else if end.is_none() {
        summary.invalid_alignment += 1;
        FilterReason::MissingAlignmentEnd
    } else if mapq.unwrap_or(0) < config.min_mapq {
        summary.below_min_mapq += 1;
        FilterReason::BelowMinMapq
    } else if config.require_full_length && !effective_full_length {
        summary.missing_full_length_evidence += 1;
        FilterReason::MissingFullLengthEvidence
    } else {
        FilterReason::Retained
    };

    let bed_record = match (&chrom, start, end, filter_reason) {
        (Some(chrom), Some(start), Some(end), FilterReason::Retained) => Some(Bed6Record {
            chrom: chrom.clone(),
            start,
            end,
            name: name.clone(),
            score: mapq.map(u32::from).unwrap_or(0),
            strand,
            extra_fields: Vec::new(),
        }),
        _ => None,
    };
    let boundary_key = bed_record.as_ref().map(|bed| BoundaryKey {
        chrom: bed.chrom.clone(),
        start: bed.start,
        end: bed.end,
        strand: bed.strand,
    });
    let cigar = if *capture_evidence {
        cigar_text(record, bam_path, ordinal, &name)?
    } else {
        String::new()
    };

    Ok(ConvertedRecord {
        ordinal,
        name,
        chrom,
        start,
        end,
        strand,
        mapq,
        five_prime_soft_clip,
        three_prime_soft_clip,
        cigar,
        evidence,
        effective_full_length,
        filter_reason,
        bed_record,
        boundary_key,
    })
}

fn terminal_soft_clips(ops: &[sam::alignment::record::cigar::Op]) -> (usize, usize) {
    let first_non_hard = ops.iter().find(|op| op.kind() != CigarKind::HardClip);
    let last_non_hard = ops.iter().rev().find(|op| op.kind() != CigarKind::HardClip);
    let leading = first_non_hard
        .filter(|op| op.kind() == CigarKind::SoftClip)
        .map(|op| op.len())
        .unwrap_or(0);
    let trailing = last_non_hard
        .filter(|op| op.kind() == CigarKind::SoftClip)
        .map(|op| op.len())
        .unwrap_or(0);
    (leading, trailing)
}

fn cigar_text(
    record: &RecordBuf,
    bam_path: &Path,
    record_ordinal: usize,
    read_name: &str,
) -> BamResult<String> {
    let mut buf = Vec::new();
    sam::io::writer::record::write_cigar(&mut buf, record.cigar()).map_err(|source| {
        BamConversionError::CigarWrite {
            path: bam_path.to_path_buf(),
            record_ordinal,
            read_name: read_name.to_owned(),
            source,
        }
    })?;
    String::from_utf8(buf).map_err(|source| BamConversionError::CigarUtf8 {
        path: bam_path.to_path_buf(),
        record_ordinal,
        read_name: read_name.to_owned(),
        source,
    })
}

fn apply_min_boundary_support(
    records: &mut [ConvertedRecord],
    minimum: usize,
    summary: &mut BamConversionSummary,
) {
    let mut support: HashMap<BoundaryKey, usize> = HashMap::new();
    for record in records.iter() {
        if record.filter_reason == FilterReason::Retained {
            if let Some(key) = &record.boundary_key {
                *support.entry(key.clone()).or_default() += 1;
            }
        }
    }

    for record in records {
        let is_below = record
            .boundary_key
            .as_ref()
            .is_some_and(|key| support.get(key).copied().unwrap_or(0) < minimum);
        if record.filter_reason == FilterReason::Retained && is_below {
            record.filter_reason = FilterReason::BelowMinBoundarySupport;
            record.bed_record = None;
            summary.below_min_boundary_support += 1;
        }
    }
}

fn read_input_evidence(path: &Path) -> BamResult<HashMap<String, InputEvidence>> {
    let mut reader = DelimitedReader::open(path, Delimiter::Tab).map_err(|error| match error {
        DelimitedError::Open { source, .. } => BamConversionError::EvidenceOpen {
            path: path.to_path_buf(),
            source,
        },
        other => BamConversionError::EvidenceRead {
            path: path.to_path_buf(),
            line: 1,
            source: io::Error::new(io::ErrorKind::InvalidData, other),
        },
    })?;
    let mut header: Option<Vec<String>> = None;
    let mut evidence = HashMap::new();

    for result in reader.records() {
        let record = result.map_err(|error| {
            let line = match &error {
                DelimitedError::Read { line, .. } => line.unwrap_or(1) as usize,
                _ => 1,
            };
            BamConversionError::EvidenceRead {
                path: path.to_path_buf(),
                line,
                source: io::Error::new(io::ErrorKind::InvalidData, error),
            }
        })?;
        let line_number = record.line_number() as usize;
        let fields = record.fields();
        if header.is_none() {
            let columns: Vec<String> = fields
                .iter()
                .map(|value| value.trim().to_ascii_lowercase())
                .collect();
            if !columns.iter().any(|column| column == "read_name") {
                return Err(BamConversionError::EvidenceMissingReadNameColumn {
                    path: path.to_path_buf(),
                    line: line_number,
                });
            }
            let mut counts = BTreeMap::new();
            for column in &columns {
                *counts.entry(column).or_insert(0usize) += 1;
            }
            if let Some((duplicate, _)) = counts.into_iter().find(|(_, count)| *count > 1) {
                return Err(BamConversionError::EvidenceDuplicateColumn {
                    path: path.to_path_buf(),
                    line: line_number,
                    column: (*duplicate).clone(),
                });
            }
            header = Some(columns);
            continue;
        }

        let columns = header.as_ref().expect("header initialized above");
        let get = |name: &str| -> &str {
            columns
                .iter()
                .position(|column| column == name)
                .and_then(|index| fields.get(index))
                .unwrap_or("")
        };
        let read_name = get("read_name").trim();
        if read_name.is_empty() || read_name == "." {
            return Err(BamConversionError::EvidenceEmptyReadName {
                path: path.to_path_buf(),
                line: line_number,
            });
        }
        let row = InputEvidence {
            poly_a: EvidenceState::parse(
                get("poly_a_evidence"),
                path,
                line_number,
                "poly_a_evidence",
            )?,
            five_prime_adapter: EvidenceState::parse(
                get("five_prime_adapter_evidence"),
                path,
                line_number,
                "five_prime_adapter_evidence",
            )?,
            three_prime_adapter: EvidenceState::parse(
                get("three_prime_adapter_evidence"),
                path,
                line_number,
                "three_prime_adapter_evidence",
            )?,
            full_length: EvidenceState::parse(
                get("full_length_evidence"),
                path,
                line_number,
                "full_length_evidence",
            )?,
            library_preparation: normalize_optional_text(get("library_preparation")),
        };
        if evidence.insert(read_name.to_owned(), row).is_some() {
            return Err(BamConversionError::EvidenceDuplicateReadName {
                path: path.to_path_buf(),
                line: line_number,
                read_name: read_name.to_owned(),
            });
        }
    }

    if header.is_none() {
        return Err(BamConversionError::EvidenceEmpty {
            path: path.to_path_buf(),
        });
    }
    Ok(evidence)
}

fn normalize_optional_text(value: &str) -> Option<String> {
    let value = value.trim();
    if value.is_empty() || value == "." {
        None
    } else {
        Some(value.to_owned())
    }
}

fn write_evidence_sidecar(
    path: &Path,
    records: &[ConvertedRecord],
    config: &BamConversionConfig,
) -> BamResult<()> {
    let mut writer = create_evidence_writer(path)?;
    for record in records {
        write_evidence_record(&mut writer, record, config).map_err(|source| {
            output_write_error(
                path,
                EVIDENCE_OUTPUT,
                format!("writing BAM record {}", record.ordinal),
                io::Error::other(source),
            )
        })?;
    }
    flush_evidence_writer(&mut writer, path)
}

fn create_evidence_writer(path: &Path) -> BamResult<DelimitedWriter<BufWriter<File>>> {
    let mut writer = DelimitedWriter::create(
        path,
        Delimiter::Tab,
        &[format!("#schema_version={BAM_EVIDENCE_SCHEMA_VERSION}")],
    )
    .map_err(|error| match error {
        DelimitedError::Create { source, .. } => BamConversionError::OutputCreate {
            path: path.to_path_buf(),
            output: EVIDENCE_OUTPUT,
            source,
        },
        other => output_write_error(
            path,
            EVIDENCE_OUTPUT,
            "writing schema header",
            io::Error::other(other),
        ),
    })?;
    writer
        .write_record([
            "schema_version",
            "record_ordinal",
            "read_name",
            "chrom",
            "start",
            "end",
            "strand",
            "mapq",
            "five_prime_soft_clip",
            "three_prime_soft_clip",
            "cigar",
            "retained",
            "filter_reason",
            "library_profile",
            "full_length_rule",
            "poly_a_evidence",
            "five_prime_adapter_evidence",
            "three_prime_adapter_evidence",
            "full_length_evidence",
            "effective_full_length",
            "library_preparation",
            "min_mapq_filter",
            "require_full_length_filter",
            "min_boundary_support_filter",
        ])
        .map_err(|source| {
            output_write_error(
                path,
                EVIDENCE_OUTPUT,
                "writing column header",
                io::Error::other(source),
            )
        })?;
    Ok(writer)
}

fn write_evidence_record<W: Write>(
    writer: &mut DelimitedWriter<W>,
    record: &ConvertedRecord,
    config: &BamConversionConfig,
) -> Result<(), DelimitedError> {
    let retained = record.filter_reason == FilterReason::Retained;
    writer.write_record([
        BAM_EVIDENCE_SCHEMA_VERSION.to_owned(),
        record.ordinal.to_string(),
        record.name.clone(),
        record.chrom.clone().unwrap_or_else(|| ".".to_owned()),
        record
            .start
            .map(|value| value.get().to_string())
            .unwrap_or_else(|| ".".to_owned()),
        record
            .end
            .map(|value| value.get().to_string())
            .unwrap_or_else(|| ".".to_owned()),
        record.strand.as_char().to_string(),
        record
            .mapq
            .map(|value| value.to_string())
            .unwrap_or_else(|| ".".to_owned()),
        record.five_prime_soft_clip.to_string(),
        record.three_prime_soft_clip.to_string(),
        record.cigar.clone(),
        if retained { "1" } else { "0" }.to_owned(),
        record.filter_reason.as_str().to_owned(),
        config.library_profile.to_string(),
        config
            .library_profile
            .inferred_full_length_rule()
            .to_owned(),
        record.evidence.poly_a.to_string(),
        record.evidence.five_prime_adapter.to_string(),
        record.evidence.three_prime_adapter.to_string(),
        record.evidence.full_length.to_string(),
        if record.effective_full_length {
            "1"
        } else {
            "0"
        }
        .to_owned(),
        record
            .evidence
            .library_preparation
            .clone()
            .unwrap_or_else(|| ".".to_owned()),
        config.min_mapq.to_string(),
        if config.require_full_length { "1" } else { "0" }.to_owned(),
        config.min_boundary_support.to_string(),
    ])
}

fn flush_evidence_writer(
    writer: &mut DelimitedWriter<BufWriter<File>>,
    path: &Path,
) -> BamResult<()> {
    writer.flush().map_err(|source| {
        output_write_error(path, EVIDENCE_OUTPUT, "flushing", io::Error::other(source))
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use noodles::sam::alignment::record::cigar::Op;

    #[test]
    fn library_profile_rules_are_explicit_and_stable() {
        assert_eq!(
            LibraryProfile::DirectRna.inferred_full_length_rule(),
            "poly_a_and_five_prime_adapter"
        );
        assert_eq!(
            LibraryProfile::DirectCdna.inferred_full_length_rule(),
            "five_prime_adapter_and_three_prime_adapter"
        );
        assert_eq!(
            "pcr-cdna".parse::<LibraryProfile>().unwrap(),
            LibraryProfile::PcrCdna
        );
    }

    #[test]
    fn terminal_soft_clips_ignore_terminal_hard_clips() {
        let ops = [
            Op::new(CigarKind::HardClip, 3),
            Op::new(CigarKind::SoftClip, 5),
            Op::new(CigarKind::Match, 20),
            Op::new(CigarKind::SoftClip, 2),
            Op::new(CigarKind::HardClip, 4),
        ];
        assert_eq!(terminal_soft_clips(&ops), (5, 2));
    }

    #[test]
    fn profile_inference_is_conservative() {
        let direct_rna = InputEvidence {
            poly_a: EvidenceState::Present,
            five_prime_adapter: EvidenceState::Present,
            ..InputEvidence::default()
        };
        assert!(direct_rna.effective_full_length(LibraryProfile::DirectRna));
        assert!(!direct_rna.effective_full_length(LibraryProfile::DirectCdna));

        let explicit_absence = InputEvidence {
            full_length: EvidenceState::Absent,
            poly_a: EvidenceState::Present,
            five_prime_adapter: EvidenceState::Present,
            ..InputEvidence::default()
        };
        assert!(!explicit_absence.effective_full_length(LibraryProfile::DirectRna));
    }

    #[test]
    fn pass_digest_rejects_same_count_and_key_set_with_changed_content() {
        fn converted_record(ordinal: usize, start: u32, mapq: u8) -> ConvertedRecord {
            let start = Coord::new(start);
            let end = Coord::new(start.get() + 10);
            let name = format!("r{ordinal}");
            let boundary_key = BoundaryKey {
                chrom: "chr1".to_owned(),
                start,
                end,
                strand: Strand::Plus,
            };
            let bed_record = Bed6Record {
                chrom: "chr1".to_owned(),
                start,
                end,
                name: name.clone(),
                score: u32::from(mapq),
                strand: Strand::Plus,
                extra_fields: Vec::new(),
            };

            ConvertedRecord {
                ordinal,
                name,
                chrom: Some("chr1".to_owned()),
                start: Some(start),
                end: Some(end),
                strand: Strand::Plus,
                mapq: Some(mapq),
                five_prime_soft_clip: 0,
                three_prime_soft_clip: 0,
                cigar: "10M".to_owned(),
                evidence: InputEvidence::default(),
                effective_full_length: false,
                filter_reason: FilterReason::Retained,
                bed_record: Some(bed_record),
                boundary_key: Some(boundary_key),
            }
        }

        fn digest(layout: &[(u32, u8)]) -> [u8; 32] {
            let mut digest = BamPassDigest::new(&sam::Header::default()).unwrap();
            let raw = RecordBuf::default();
            for (index, &(start, mapq)) in layout.iter().enumerate() {
                digest.update_decoded(&raw, &converted_record(index + 1, start, mapq));
            }
            digest.finish()
        }

        let first = digest(&[(10, 30), (10, 30), (20, 30)]);
        let changed_multiplicity = digest(&[(10, 30), (20, 30), (20, 30)]);
        let changed_output_field = digest(&[(10, 30), (10, 31), (20, 30)]);

        assert_ne!(first, changed_multiplicity);
        assert_ne!(first, changed_output_field);
        ensure_pass_digests_match(Path::new("reads.bam"), first, first).unwrap();

        for second in [changed_multiplicity, changed_output_field] {
            match ensure_pass_digests_match(Path::new("reads.bam"), first, second).unwrap_err() {
                BamConversionError::ContentChangedBetweenPasses {
                    path,
                    first_digest,
                    second_digest,
                } => {
                    assert_eq!(path, Path::new("reads.bam"));
                    assert_eq!(first_digest.len(), 64);
                    assert_eq!(second_digest.len(), 64);
                    assert_ne!(first_digest, second_digest);
                }
                other => panic!("unexpected error: {other:?}"),
            }
        }
    }
}
