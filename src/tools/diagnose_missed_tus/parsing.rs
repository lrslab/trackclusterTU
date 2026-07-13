use std::collections::{HashMap, HashSet};
use std::io::BufRead;
use std::path::Path;

use crate::io::bed::{BedError, BedParseError};
pub(super) use crate::io::gff::GeneRecord;
use crate::model::{Interval, Strand};
use crate::tools::domain_records::{read_bed6_as_genes, read_named_bed6, Bed6Validation};
use crate::tools::membership::{read_memberships, MembershipRecord};
use crate::tools::read_rejections::{
    LocatedRead, ReadLoadOutcome, ReadLocation, ReadRejection, ReadRejectionReason,
};
use crate::tu::ReadRecord;
use anyhow::{Context, Result};

#[derive(Clone, Debug)]
pub(super) struct ExistingTu {
    pub(super) id: String,
    pub(super) contig: String,
    pub(super) strand: Strand,
    pub(super) interval: Interval,
}

pub(super) type ExistingMembershipRow = MembershipRecord;

#[derive(Debug)]
pub(super) struct BedReadLoad {
    pub(super) outcome: ReadLoadOutcome,
    /// Every syntactically decoded, non-empty read ID, including IDs later
    /// rejected by validation/filtering. Rescue membership foreign-key checks
    /// use this set so a membership row for a deliberately skipped read does
    /// not turn a recoverable read problem into a global failure.
    pub(super) input_read_ids: HashSet<String>,
}

pub(super) fn load_detection_inputs(
    input: &Path,
    existing_tu: &Path,
    annotation_bed: Option<&Path>,
    min_read_len: Option<u32>,
) -> Result<(BedReadLoad, Vec<ExistingTu>, Vec<GeneRecord>)> {
    let reads = read_reads_bed6(input, min_read_len)?;
    let existing_tus = read_existing_tus_bed6(existing_tu)?;
    let genes = match annotation_bed {
        Some(path) => read_gene_bed6(path)?,
        None => Vec::new(),
    };
    Ok((reads, existing_tus, genes))
}

pub(super) fn read_existing_membership(path: &Path) -> Result<Vec<ExistingMembershipRow>> {
    read_memberships(path).map_err(anyhow::Error::new)
}

pub(super) fn read_reads_bed6(path: &Path, min_read_len: Option<u32>) -> Result<BedReadLoad> {
    // Capture IDs independently of full BED decoding. This lets rescue regard a
    // membership row for a malformed-but-identifiable BED record as referring
    // to an input read, while the malformed read itself remains quarantined.
    let input_read_ids = bed_input_read_ids(path)?;
    let mut reader = crate::io::bed::read_bed6(path)?;
    let mut parsed = Vec::new();
    let mut rejections = Vec::new();
    let mut total_records = 0u64;

    while let Some(result) = reader.next() {
        total_records = total_records.saturating_add(1);
        let line = reader.line_number() as u64;
        let mut location = ReadLocation::new(path, "bed6");
        location.record = Some(total_records);
        location.line = Some(line);

        match result {
            Ok(record) => {
                parsed.push(LocatedRead {
                    read: ReadRecord {
                        contig: record.chrom,
                        strand: record.strand,
                        interval: Interval::new(record.start, record.end)
                            .expect("BED6 reader already validates interval ordering"),
                        id: record.name,
                    },
                    location,
                    full_length_evidence: false,
                });
            }
            Err(BedError::Parse { source, .. }) => {
                let reason = bed_parse_rejection_reason(&source);
                rejections.push(ReadRejection::new(
                    location,
                    None,
                    "parse",
                    reason,
                    source.to_string(),
                ));
            }
            Err(error @ BedError::IoRead { .. }) => return Err(anyhow::Error::new(error)),
            Err(error @ BedError::IoWrite { .. }) => return Err(anyhow::Error::new(error)),
        }
    }

    // Match the main cluster loader's ordering exactly: validation precedes
    // duplicate detection, so an invalid row cannot poison a valid row that
    // happens to carry the same raw ID.
    let mut validated = Vec::with_capacity(parsed.len());
    for located in parsed {
        let read = &located.read;
        let rejection = if read.id.is_empty() {
            Some((
                "validate",
                ReadRejectionReason::EmptyReadId,
                "read ID must not be empty".to_owned(),
            ))
        } else if read.contig.is_empty() {
            Some((
                "validate",
                ReadRejectionReason::EmptyContig,
                "read contig must not be empty".to_owned(),
            ))
        } else if read.strand == Strand::Unknown {
            Some((
                "validate",
                ReadRejectionReason::UnknownStrand,
                "boundary-aware diagnosis/rescue requires '+' or '-' strand".to_owned(),
            ))
        } else {
            None
        };

        if let Some((stage, reason, detail)) = rejection {
            rejections.push(ReadRejection::new(
                located.location,
                (!read.id.is_empty()).then(|| read.id.clone()),
                stage,
                reason,
                detail,
            ));
        } else {
            validated.push(located);
        }
    }

    let duplicate_counts: HashMap<String, usize> =
        validated
            .iter()
            .fold(HashMap::new(), |mut counts, located| {
                *counts.entry(located.read.id.clone()).or_default() += 1;
                counts
            });
    let mut deduplicated = Vec::with_capacity(validated.len());
    for located in validated {
        let read = &located.read;
        let count = duplicate_counts.get(read.id.as_str()).copied().unwrap_or(1);
        if count > 1 {
            rejections.push(ReadRejection::new(
                located.location,
                Some(read.id.clone()),
                "deduplicate",
                ReadRejectionReason::DuplicateReadId,
                format!(
                    "duplicate read ID {:?}; all {count} copies were quarantined",
                    read.id
                ),
            ));
        } else {
            deduplicated.push(located);
        }
    }

    let mut reads = Vec::with_capacity(deduplicated.len());
    for located in deduplicated {
        let read = &located.read;
        let rejection = if read.interval.is_empty() {
            Some((
                "filter",
                ReadRejectionReason::ZeroLength,
                "zero-length read interval".to_owned(),
            ))
        } else if min_read_len.is_some_and(|minimum| read.interval.len() < minimum) {
            let minimum = min_read_len.expect("minimum checked above");
            Some((
                "filter",
                ReadRejectionReason::BelowMinReadLen,
                format!(
                    "read length {} is below minimum {minimum}",
                    read.interval.len()
                ),
            ))
        } else {
            None
        };

        if let Some((stage, reason, detail)) = rejection {
            rejections.push(ReadRejection::new(
                located.location,
                (!read.id.is_empty()).then(|| read.id.clone()),
                stage,
                reason,
                detail,
            ));
        } else {
            reads.push(located);
        }
    }

    rejections.sort_by(|a, b| {
        a.location
            .source_path
            .cmp(&b.location.source_path)
            .then_with(|| a.location.line.cmp(&b.location.line))
            .then_with(|| a.read_id.cmp(&b.read_id))
    });

    Ok(BedReadLoad {
        outcome: ReadLoadOutcome {
            reads,
            rejections,
            total_records,
        },
        input_read_ids,
    })
}

fn bed_input_read_ids(path: &Path) -> Result<HashSet<String>> {
    let file = std::fs::File::open(path)
        .with_context(|| format!("failed to scan read IDs from BED6 {path:?}"))?;
    let mut reader = std::io::BufReader::new(file);
    let mut line = Vec::new();
    let mut read_ids = HashSet::new();

    loop {
        line.clear();
        let bytes_read = reader
            .read_until(b'\n', &mut line)
            .with_context(|| format!("failed to scan read IDs from BED6 {path:?}"))?;
        if bytes_read == 0 {
            break;
        }
        let content = trim_ascii_end(&line);
        let Some(first_non_whitespace) =
            content.iter().position(|byte| !byte.is_ascii_whitespace())
        else {
            continue;
        };
        if content[first_non_whitespace] == b'#' {
            continue;
        }

        let raw_id = if content.contains(&b'\t') {
            content.split(|byte| *byte == b'\t').nth(3)
        } else {
            content
                .split(|byte| byte.is_ascii_whitespace())
                .filter(|field| !field.is_empty())
                .nth(3)
        };
        let Some(read_id) = raw_id
            .and_then(|raw_id| std::str::from_utf8(raw_id).ok())
            .filter(|read_id| !read_id.is_empty())
        else {
            continue;
        };
        read_ids.insert(read_id.to_owned());
    }

    Ok(read_ids)
}

fn trim_ascii_end(mut bytes: &[u8]) -> &[u8] {
    while bytes.last().is_some_and(|byte| byte.is_ascii_whitespace()) {
        bytes = &bytes[..bytes.len() - 1];
    }
    bytes
}

fn bed_parse_rejection_reason(error: &BedParseError) -> ReadRejectionReason {
    match error {
        BedParseError::InvalidUtf8 { .. } => ReadRejectionReason::InvalidUtf8,
        BedParseError::TooFewColumnsBed6 { .. } | BedParseError::TooFewColumns { .. } => {
            ReadRejectionReason::TooFewColumns
        }
        BedParseError::InvalidInt { .. } | BedParseError::InvalidBlockCount { .. } => {
            ReadRejectionReason::InvalidInteger
        }
        BedParseError::Strand(_) => ReadRejectionReason::InvalidStrand,
        BedParseError::Interval(_) => ReadRejectionReason::InvalidInterval,
        BedParseError::BlockListLengthMismatch { .. }
        | BedParseError::BlockOverflow
        | BedParseError::Transcript(_) => ReadRejectionReason::InvalidBed12Blocks,
    }
}

pub(super) fn read_existing_tus_bed6(path: &Path) -> Result<Vec<ExistingTu>> {
    read_named_bed6(
        path,
        Bed6Validation {
            id_label: "existing TU",
            require_nonempty_id: true,
            require_unique_ids: true,
            unknown_strand_context: Some("boundary-aware diagnosis/rescue requires '+' or '-'"),
            invalid_interval_label: "TU",
        },
    )
    .map(|records| {
        records
            .into_iter()
            .map(|record| ExistingTu {
                id: record.id,
                contig: record.contig,
                strand: record.strand,
                interval: record.interval,
            })
            .collect()
    })
}

pub(super) fn read_gene_bed6(path: &Path) -> Result<Vec<GeneRecord>> {
    read_bed6_as_genes(
        path,
        Bed6Validation {
            id_label: "annotation feature",
            require_nonempty_id: false,
            require_unique_ids: false,
            unknown_strand_context: Some("boundary-aware diagnosis/rescue requires '+' or '-'"),
            invalid_interval_label: "gene",
        },
    )
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::io::delimited::{DelimitedReader, Delimiter};
    use crate::tools::membership::{parse_memberships, MembershipError};

    fn temp_bed(label: &str, bytes: &[u8]) -> std::path::PathBuf {
        let path = std::env::temp_dir().join(format!(
            "trackclustertu_diagnose_{label}_{}_{}.bed",
            std::process::id(),
            std::time::SystemTime::now()
                .duration_since(std::time::UNIX_EPOCH)
                .unwrap()
                .as_nanos()
        ));
        std::fs::write(&path, bytes).unwrap();
        path
    }

    #[test]
    fn bed_read_errors_are_quarantined_and_later_good_rows_survive() {
        let path = temp_bed(
            "tolerant",
            concat!(
                "# metadata\n",
                "chr1\t0\t10\tgood1\t0\t+\n",
                "chr1\tbad\t20\tbad_integer\t0\t+\n",
                "\t20\t30\tempty_contig\t0\t+\n",
                "chr1\t30\t40\t\t0\t+\n",
                "chr1\t40\t50\tunknown\t0\t.\n",
                "chr1\t60\t60\tzero\t0\t+\n",
                "chr1\t70\t72\tshort\t0\t+\n",
                "chr1\t80\t90\tdup\t0\t+\n",
                "chr1\t81\t91\tdup\t0\t+\n",
                "chr1\t100\t110\tgood2\t0\t-\n",
            )
            .as_bytes(),
        );

        let load = read_reads_bed6(&path, Some(5)).unwrap();
        assert_eq!(load.outcome.total_records, 10);
        assert_eq!(
            load.outcome
                .reads
                .iter()
                .map(|located| located.read.id.as_str())
                .collect::<Vec<_>>(),
            ["good1", "good2"]
        );
        assert_eq!(load.outcome.rejections.len(), 8);
        assert_eq!(
            load.outcome
                .rejections
                .iter()
                .map(|rejection| rejection.reason)
                .collect::<Vec<_>>(),
            [
                ReadRejectionReason::InvalidInteger,
                ReadRejectionReason::EmptyContig,
                ReadRejectionReason::EmptyReadId,
                ReadRejectionReason::UnknownStrand,
                ReadRejectionReason::ZeroLength,
                ReadRejectionReason::BelowMinReadLen,
                ReadRejectionReason::DuplicateReadId,
                ReadRejectionReason::DuplicateReadId,
            ]
        );
        assert!(load.input_read_ids.contains("dup"));
        assert!(load.input_read_ids.contains("short"));
        assert!(load.input_read_ids.contains("bad_integer"));
        assert_eq!(
            load.outcome.total_records as usize,
            load.outcome.reads.len() + load.outcome.rejections.len()
        );

        let _ = std::fs::remove_file(path);
    }

    #[test]
    fn invalid_utf8_bed_row_is_recorded_without_hiding_following_read() {
        let mut bytes = b"chr1\t0\t10\tgood1\t0\t+\nchr1\t20\t30\t".to_vec();
        bytes.push(0xff);
        bytes.extend_from_slice(b"\t0\t+\nchr1\t40\t50\tgood2\t0\t-\n");
        let path = temp_bed("invalid_utf8", &bytes);

        let load = read_reads_bed6(&path, None).unwrap();
        assert_eq!(load.outcome.total_records, 3);
        assert_eq!(load.outcome.reads.len(), 2);
        assert_eq!(load.outcome.reads[1].read.id, "good2");
        assert_eq!(load.outcome.rejections.len(), 1);
        assert_eq!(
            load.outcome.rejections[0].reason,
            ReadRejectionReason::InvalidUtf8
        );
        assert_eq!(load.outcome.rejections[0].location.line, Some(2));

        let _ = std::fs::remove_file(path);
    }

    #[test]
    fn invalid_same_id_rows_do_not_duplicate_quarantine_a_valid_read() {
        let path = temp_bed(
            "invalid_same_id",
            concat!(
                "chr1\t0\t10\tr1\t0\t+\n",
                "chr1\t20\t30\tr1\t0\t.\n",
                "chr1\t50\t40\tr1\t0\t+\n",
                "chr1\t60\t70\tr2\t0\t-\n",
            )
            .as_bytes(),
        );

        let load = read_reads_bed6(&path, None).unwrap();
        assert_eq!(
            load.outcome
                .reads
                .iter()
                .map(|located| located.read.id.as_str())
                .collect::<Vec<_>>(),
            ["r1", "r2"]
        );
        assert_eq!(load.outcome.rejections.len(), 2);
        assert_eq!(
            load.outcome
                .rejections
                .iter()
                .map(|rejection| rejection.reason)
                .collect::<Vec<_>>(),
            [
                ReadRejectionReason::UnknownStrand,
                ReadRejectionReason::InvalidInterval,
            ]
        );
        assert!(load
            .outcome
            .rejections
            .iter()
            .all(|rejection| rejection.reason != ReadRejectionReason::DuplicateReadId));
        assert!(load.input_read_ids.contains("r1"));

        let _ = std::fs::remove_file(path);
    }

    #[test]
    fn v2_no_hard_assignments_are_typed_and_duplicates_still_fail() {
        let header = concat!(
            "#trackclustertu_membership_schema=v2\n",
            "#columns=read_id\ttu_id\tscore1\tscore2\tschema_version\tassignment_status\n",
        );
        let rows_text = format!(
            "{header}ambiguous\t.\t.\t.\tv2\tambiguous\tTU1\tTU2\t.\t.\t.\t.\t.\t.\t.\t.\t0\t0\t.\tfalse\nunassigned\t.\t.\t.\tv2\tunassigned\t.\t.\t.\t.\t.\t.\t.\t.\t.\t.\t.\t0\t.\tfalse\nhard\tTU1\t1\t1\tv2\tunique\tTU1\t.\t.\t.\t0\t0\t.\t.\t1\t.\t.\t1\tTU1:1\tfalse\n"
        );
        let reader = DelimitedReader::from_reader(
            Path::new("membership.tsv"),
            Delimiter::Tab,
            rows_text.as_bytes(),
        );
        let rows = parse_memberships(reader, true).unwrap();
        assert_eq!(rows[0].line_number, 3);
        assert_eq!(rows[1].line_number, 4);
        assert_eq!(rows[2].line_number, 5);
        assert_eq!(rows[0].hard_tu_id, None);
        assert_eq!(rows[1].hard_tu_id, None);
        assert_eq!(rows[2].hard_tu_id.as_deref(), Some("TU1"));

        let duplicate = format!(
            "{header}same\t.\t.\t.\tv2\tambiguous\tTU1\tTU2\t.\t.\t.\t.\t.\t.\t.\t.\t0\t0\t.\tfalse\nsame\t.\t.\t.\tv2\tunassigned\t.\t.\t.\t.\t.\t.\t.\t.\t.\t.\t.\t0\t.\tfalse\n"
        );
        let reader = DelimitedReader::from_reader(
            Path::new("membership.tsv"),
            Delimiter::Tab,
            duplicate.as_bytes(),
        );
        let error = parse_memberships(reader, true).unwrap_err();
        assert!(matches!(error, MembershipError::DuplicateRead { .. }));
    }
}
