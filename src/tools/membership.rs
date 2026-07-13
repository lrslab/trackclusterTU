//! Shared version-aware membership table model, parser, and writer.

use std::collections::{HashMap, HashSet};
use std::io::{BufRead, BufReader, Read, Write};
use std::path::Path;

use thiserror::Error;

use crate::io::delimited::{DelimitedError, DelimitedReader, DelimitedWriter, Delimiter};

pub(crate) const V2_COLUMNS: [&str; 20] = [
    "read_id",
    "tu_id",
    "score1",
    "score2",
    "schema_version",
    "assignment_status",
    "best_candidate_tu_id",
    "second_candidate_tu_id",
    "second_score1",
    "second_score2",
    "best_five_prime_delta_bp",
    "best_three_prime_delta_bp",
    "second_five_prime_delta_bp",
    "second_three_prime_delta_bp",
    "best_assignment_score",
    "second_assignment_score",
    "score_margin",
    "primary_weight",
    "fractional_assignments",
    "full_length_evidence",
];

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) enum MembershipStatus {
    Legacy,
    Unique,
    Ambiguous,
    Partial,
    Unassigned,
}

impl MembershipStatus {
    pub(crate) const fn as_str(self) -> &'static str {
        match self {
            Self::Legacy => "legacy",
            Self::Unique => "unique",
            Self::Ambiguous => "ambiguous",
            Self::Partial => "partial",
            Self::Unassigned => "unassigned",
        }
    }
}

#[derive(Clone, Debug, PartialEq)]
pub(crate) struct WeightedTu {
    pub(crate) tu_id: String,
    pub(crate) weight: f64,
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub(crate) struct V2MembershipFields {
    /// Whether the source row carried the complete canonical v2 field set.
    ///
    /// Early v2 files were allowed to contain only the six required prefix
    /// columns. Keeping that distinction lets a parse/write cycle preserve
    /// those compatibility rows without manufacturing invalid audit fields.
    pub(crate) complete: bool,
    pub(crate) best_candidate_tu_id: Option<String>,
    pub(crate) second_candidate_tu_id: Option<String>,
    pub(crate) second_score1: String,
    pub(crate) second_score2: String,
    pub(crate) best_five_prime_delta_bp: String,
    pub(crate) best_three_prime_delta_bp: String,
    pub(crate) second_five_prime_delta_bp: String,
    pub(crate) second_three_prime_delta_bp: String,
    pub(crate) best_assignment_score: String,
    pub(crate) second_assignment_score: String,
    pub(crate) score_margin: String,
    pub(crate) primary_weight: String,
    pub(crate) full_length_evidence: bool,
}

#[derive(Clone, Debug, PartialEq)]
pub(crate) struct MembershipRecord {
    pub(crate) read_id: String,
    pub(crate) hard_tu_id: Option<String>,
    pub(crate) score1: String,
    pub(crate) score2: String,
    pub(crate) status: MembershipStatus,
    pub(crate) fractional_assignments: Vec<WeightedTu>,
    pub(crate) v2: Option<V2MembershipFields>,
    pub(crate) line_number: usize,
}

impl MembershipRecord {
    pub(crate) fn legacy(read_id: String, tu_id: String, score1: String, score2: String) -> Self {
        Self {
            read_id,
            hard_tu_id: Some(tu_id),
            score1,
            score2,
            status: MembershipStatus::Legacy,
            fractional_assignments: Vec::new(),
            v2: None,
            line_number: 0,
        }
    }

    pub(crate) fn candidate_tu_ids(&self) -> impl Iterator<Item = &str> {
        let candidates = self.v2.iter().flat_map(|v2| {
            [
                v2.best_candidate_tu_id.as_deref(),
                v2.second_candidate_tu_id.as_deref(),
            ]
            .into_iter()
            .flatten()
        });
        self.hard_tu_id
            .iter()
            .map(String::as_str)
            .chain(candidates)
            .chain(
                self.fractional_assignments
                    .iter()
                    .map(|assignment| assignment.tu_id.as_str()),
            )
    }

    pub(crate) fn contributions(&self) -> Vec<MembershipContribution<'_>> {
        match self.status {
            MembershipStatus::Legacy => self
                .hard_tu_id
                .as_deref()
                .map(|tu_id| {
                    vec![MembershipContribution {
                        tu_id,
                        unique: false,
                        full_length_evidence: false,
                        total: true,
                        fractional_weight: 1.0,
                    }]
                })
                .unwrap_or_default(),
            MembershipStatus::Unique | MembershipStatus::Partial => {
                let Some(tu_id) = self.hard_tu_id.as_deref() else {
                    return Vec::new();
                };
                let fractional_weight = self
                    .fractional_assignments
                    .iter()
                    .find(|assignment| assignment.tu_id == tu_id)
                    .map_or(1.0, |assignment| assignment.weight);
                vec![MembershipContribution {
                    tu_id,
                    unique: self.status == MembershipStatus::Unique,
                    full_length_evidence: self
                        .v2
                        .as_ref()
                        .is_some_and(|v2| v2.full_length_evidence),
                    total: true,
                    fractional_weight,
                }]
            }
            MembershipStatus::Ambiguous => self
                .fractional_assignments
                .iter()
                .map(|assignment| MembershipContribution {
                    tu_id: assignment.tu_id.as_str(),
                    unique: false,
                    full_length_evidence: false,
                    total: false,
                    fractional_weight: assignment.weight,
                })
                .collect(),
            MembershipStatus::Unassigned => Vec::new(),
        }
    }
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub(crate) struct MembershipContribution<'a> {
    pub(crate) tu_id: &'a str,
    pub(crate) unique: bool,
    pub(crate) full_length_evidence: bool,
    pub(crate) total: bool,
    pub(crate) fractional_weight: f64,
}

#[derive(Debug, Error)]
pub(crate) enum MembershipError {
    #[error(transparent)]
    Delimited(#[from] DelimitedError),

    #[error("{path:?}:{line}: expected at least {expected} membership columns, got {got}")]
    TooFewColumns {
        path: std::path::PathBuf,
        line: usize,
        expected: usize,
        got: usize,
    },

    #[error("{path:?}:{line}: membership read ID must not be empty")]
    EmptyReadId {
        path: std::path::PathBuf,
        line: usize,
    },

    #[error("{path:?}:{line}: membership TU ID must not be empty")]
    EmptyTuId {
        path: std::path::PathBuf,
        line: usize,
    },

    #[error(
        "{path:?}:{line}: unsupported membership schema_version {schema:?} for read {read_id:?}; expected v2 or an absent legacy marker"
    )]
    UnsupportedSchema {
        path: std::path::PathBuf,
        line: usize,
        read_id: String,
        schema: String,
    },

    #[error("{path:?}:{line}: invalid assignment_status {status:?} for read {read_id:?}")]
    InvalidStatus {
        path: std::path::PathBuf,
        line: usize,
        read_id: String,
        status: String,
    },

    #[error("{path:?}:{line}: {status} read {read_id:?} requires a hard TU ID")]
    MissingHardAssignment {
        path: std::path::PathBuf,
        line: usize,
        read_id: String,
        status: &'static str,
    },

    #[error("{path:?}:{line}: {status} read {read_id:?} cannot have a hard TU ID")]
    UnexpectedHardAssignment {
        path: std::path::PathBuf,
        line: usize,
        read_id: String,
        status: &'static str,
    },

    #[error("{path:?}:{line}: unassigned read {read_id:?} cannot have fractional assignments")]
    UnassignedFractions {
        path: std::path::PathBuf,
        line: usize,
        read_id: String,
    },

    #[error("{path:?}:{line}: invalid full_length_evidence value {value:?}")]
    InvalidFullLengthEvidence {
        path: std::path::PathBuf,
        line: usize,
        value: String,
    },

    #[error(
        "{path:?}:{line}: invalid fractional assignment {value:?} for read {read_id:?}; expected comma-separated tu_id:weight entries, and TU IDs containing ',' cannot be represented"
    )]
    InvalidFraction {
        path: std::path::PathBuf,
        line: usize,
        read_id: String,
        value: String,
    },

    #[error("{path:?}:{line}: empty or duplicate fractional TU {tu_id:?} for read {read_id:?}")]
    DuplicateFractionalTu {
        path: std::path::PathBuf,
        line: usize,
        read_id: String,
        tu_id: String,
    },

    #[error("{path:?}:{line}: invalid fractional weight {value:?} for read {read_id:?}")]
    InvalidFractionalWeight {
        path: std::path::PathBuf,
        line: usize,
        read_id: String,
        value: String,
    },

    #[error(
        "{path:?}:{line}: fractional TU ID {tu_id:?} for read {read_id:?} contains ',', which the tu_id:weight mini-format cannot represent"
    )]
    ReservedFractionalTuId {
        path: std::path::PathBuf,
        line: usize,
        read_id: String,
        tu_id: String,
    },

    #[error(
        "cannot write membership for read {read_id:?}: fractional TU ID {tu_id:?} contains ',', which the tu_id:weight mini-format cannot represent"
    )]
    ReservedFractionalTuIdWrite { read_id: String, tu_id: String },

    #[error(
        "{path:?}:{line}: fractional TU ID for read {read_id:?} cannot be '.', which is reserved for no assignment"
    )]
    FractionalSentinelTuId {
        path: std::path::PathBuf,
        line: usize,
        read_id: String,
    },

    #[error(
        "cannot write membership for read {read_id:?}: fractional TU ID cannot be '.', which is reserved for no assignment"
    )]
    FractionalSentinelTuIdWrite { read_id: String },

    #[error("{path:?}:{line}: invalid primary_weight {value:?} for read {read_id:?}")]
    InvalidPrimaryWeight {
        path: std::path::PathBuf,
        line: usize,
        read_id: String,
        value: String,
    },

    #[error("{path:?}:{line}: inconsistent v2 membership for read {read_id:?}: {reason}")]
    InconsistentV2 {
        path: std::path::PathBuf,
        line: usize,
        read_id: String,
        reason: String,
    },

    #[error("{path:?}:{line}: fractional weights for read {read_id:?} must sum exactly to 1.0, got {sum}")]
    FractionalWeightSum {
        path: std::path::PathBuf,
        line: usize,
        read_id: String,
        sum: f64,
    },

    #[error("{path:?}:{line}: duplicate membership for read ID {read_id:?}; first assigned at line {first_line}")]
    DuplicateRead {
        path: std::path::PathBuf,
        line: usize,
        read_id: String,
        first_line: usize,
    },

    #[error("failed to scan membership metadata {path:?}: {source}")]
    MetadataRead {
        path: std::path::PathBuf,
        source: std::io::Error,
    },
}

/// Read unquoted membership comment lines while respecting quoted multiline fields.
pub(crate) fn read_membership_metadata(path: &Path) -> Result<Vec<String>, MembershipError> {
    let file = std::fs::File::open(path).map_err(|source| MembershipError::MetadataRead {
        path: path.to_path_buf(),
        source,
    })?;
    let mut reader = BufReader::new(file);
    scan_membership_metadata(&mut reader).map_err(|source| MembershipError::MetadataRead {
        path: path.to_path_buf(),
        source,
    })
}

fn scan_membership_metadata<R: BufRead>(reader: &mut R) -> std::io::Result<Vec<String>> {
    let mut line = Vec::new();
    let mut metadata = Vec::new();
    let mut in_quotes = false;
    let mut at_field_start = true;

    loop {
        line.clear();
        if reader.read_until(b'\n', &mut line)? == 0 {
            break;
        }

        let content_end = line
            .iter()
            .rposition(|byte| !matches!(byte, b'\n' | b'\r'))
            .map_or(0, |index| index + 1);
        let content = &line[..content_end];
        let first_non_whitespace = content.iter().position(|byte| !byte.is_ascii_whitespace());
        if !in_quotes && first_non_whitespace.is_some_and(|index| content[index] == b'#') {
            let start = first_non_whitespace.expect("comment position checked above");
            metadata.push(String::from_utf8_lossy(&content[start..]).into_owned());
            at_field_start = true;
            continue;
        }

        let mut index = 0;
        while index < line.len() {
            let byte = line[index];
            if in_quotes {
                if byte == b'"' {
                    if line.get(index + 1) == Some(&b'"') {
                        index += 2;
                        continue;
                    }
                    in_quotes = false;
                }
            } else if byte == b'"' && at_field_start {
                in_quotes = true;
                at_field_start = false;
            } else {
                at_field_start = byte == b'\t' || byte == b'\n';
            }
            index += 1;
        }
    }

    Ok(metadata)
}

pub(crate) fn read_memberships(path: &Path) -> Result<Vec<MembershipRecord>, MembershipError> {
    let reader = DelimitedReader::open(path, Delimiter::Tab)?;
    parse_memberships(reader, true)
}

pub(crate) fn parse_memberships<R: Read>(
    mut reader: DelimitedReader<R>,
    reject_duplicate_reads: bool,
) -> Result<Vec<MembershipRecord>, MembershipError> {
    let path = reader.path().to_path_buf();
    let mut records = Vec::new();
    let mut first_line_by_read = HashMap::new();
    for table_record in reader.records() {
        let table_record = table_record?;
        let line = table_record.line_number() as usize;
        let fields = table_record.fields();
        if fields.len() < 2 {
            return Err(MembershipError::TooFewColumns {
                path: path.clone(),
                line,
                expected: 2,
                got: fields.len(),
            });
        }
        let read_id = fields.get(0).unwrap_or("").to_owned();
        let raw_tu_id = fields.get(1).unwrap_or("").to_owned();
        if read_id.is_empty() {
            return Err(MembershipError::EmptyReadId {
                path: path.clone(),
                line,
            });
        }
        if raw_tu_id.is_empty() {
            return Err(MembershipError::EmptyTuId {
                path: path.clone(),
                line,
            });
        }
        if reject_duplicate_reads {
            if let Some(first_line) = first_line_by_read.insert(read_id.clone(), line) {
                return Err(MembershipError::DuplicateRead {
                    path: path.clone(),
                    line,
                    read_id,
                    first_line,
                });
            }
        }

        let score1 = fields.get(2).unwrap_or(".").to_owned();
        let score2 = fields.get(3).unwrap_or(".").to_owned();
        let schema = fields.get(4).unwrap_or("");
        let is_v2 = match schema {
            "v2" => true,
            "" | "." => false,
            _ => {
                return Err(MembershipError::UnsupportedSchema {
                    path: path.clone(),
                    line,
                    read_id,
                    schema: schema.to_owned(),
                })
            }
        };
        let mut record = if !is_v2 {
            if raw_tu_id == "." {
                return Err(MembershipError::MissingHardAssignment {
                    path: path.clone(),
                    line,
                    read_id,
                    status: "legacy",
                });
            }
            MembershipRecord::legacy(read_id, raw_tu_id, score1, score2)
        } else {
            if fields.len() < 6 {
                return Err(MembershipError::TooFewColumns {
                    path: path.clone(),
                    line,
                    expected: 6,
                    got: fields.len(),
                });
            }
            if fields.len() > 6 && fields.len() < V2_COLUMNS.len() {
                return Err(MembershipError::TooFewColumns {
                    path: path.clone(),
                    line,
                    expected: V2_COLUMNS.len(),
                    got: fields.len(),
                });
            }
            parse_v2(&path, fields, line, read_id, raw_tu_id, score1, score2)?
        };
        record.line_number = line;
        records.push(record);
    }
    Ok(records)
}

fn parse_v2(
    path: &Path,
    fields: &csv::StringRecord,
    line: usize,
    read_id: String,
    raw_tu_id: String,
    score1: String,
    score2: String,
) -> Result<MembershipRecord, MembershipError> {
    let status = match fields.get(5).unwrap_or("") {
        "unique" => MembershipStatus::Unique,
        "ambiguous" => MembershipStatus::Ambiguous,
        "partial" => MembershipStatus::Partial,
        "unassigned" => MembershipStatus::Unassigned,
        status => {
            return Err(MembershipError::InvalidStatus {
                path: path.to_path_buf(),
                line,
                read_id,
                status: status.to_owned(),
            })
        }
    };
    let hard_tu_id = (raw_tu_id != ".").then_some(raw_tu_id);
    match status {
        MembershipStatus::Unique | MembershipStatus::Partial if hard_tu_id.is_none() => {
            return Err(MembershipError::MissingHardAssignment {
                path: path.to_path_buf(),
                line,
                read_id,
                status: status.as_str(),
            });
        }
        MembershipStatus::Ambiguous | MembershipStatus::Unassigned if hard_tu_id.is_some() => {
            return Err(MembershipError::UnexpectedHardAssignment {
                path: path.to_path_buf(),
                line,
                read_id,
                status: status.as_str(),
            });
        }
        _ => {}
    }
    let best_candidate_tu_id = optional_id(fields.get(6));
    let second_candidate_tu_id = optional_id(fields.get(7));
    let fractional_text = fields.get(18).unwrap_or(".");
    if !fractional_text.is_empty() && fractional_text != "." {
        for tu_id in [
            hard_tu_id.as_deref(),
            best_candidate_tu_id.as_deref(),
            second_candidate_tu_id.as_deref(),
        ]
        .into_iter()
        .flatten()
        {
            validate_fractional_tu_id(path, line, read_id.as_str(), tu_id)?;
        }
    }
    let fractional_assignments =
        parse_fractional_assignments(path, fractional_text, line, read_id.as_str())?;
    let primary_weight_text = fields.get(17).unwrap_or("0").to_owned();
    let primary_weight = parse_primary_weight(path, line, read_id.as_str(), &primary_weight_text)?;
    if fields.len() >= V2_COLUMNS.len() {
        validate_v2_invariants(
            path,
            line,
            read_id.as_str(),
            status,
            hard_tu_id.as_deref(),
            best_candidate_tu_id.as_deref(),
            second_candidate_tu_id.as_deref(),
            &fractional_assignments,
            primary_weight,
        )?;
    }
    let full_length_evidence = match fields.get(19) {
        None => false,
        Some(value) => {
            value
                .parse::<bool>()
                .map_err(|_| MembershipError::InvalidFullLengthEvidence {
                    path: path.to_path_buf(),
                    line,
                    value: value.to_owned(),
                })?
        }
    };

    Ok(MembershipRecord {
        read_id,
        hard_tu_id,
        score1,
        score2,
        status,
        fractional_assignments,
        v2: Some(V2MembershipFields {
            complete: fields.len() >= V2_COLUMNS.len(),
            best_candidate_tu_id,
            second_candidate_tu_id,
            second_score1: fields.get(8).unwrap_or(".").to_owned(),
            second_score2: fields.get(9).unwrap_or(".").to_owned(),
            best_five_prime_delta_bp: fields.get(10).unwrap_or(".").to_owned(),
            best_three_prime_delta_bp: fields.get(11).unwrap_or(".").to_owned(),
            second_five_prime_delta_bp: fields.get(12).unwrap_or(".").to_owned(),
            second_three_prime_delta_bp: fields.get(13).unwrap_or(".").to_owned(),
            best_assignment_score: fields.get(14).unwrap_or(".").to_owned(),
            second_assignment_score: fields.get(15).unwrap_or(".").to_owned(),
            score_margin: fields.get(16).unwrap_or(".").to_owned(),
            primary_weight: primary_weight_text,
            full_length_evidence,
        }),
        line_number: line,
    })
}

fn optional_id(value: Option<&str>) -> Option<String> {
    value
        .filter(|value| !value.is_empty() && *value != ".")
        .map(str::to_owned)
}

fn validate_fractional_tu_id(
    path: &Path,
    line: usize,
    read_id: &str,
    tu_id: &str,
) -> Result<(), MembershipError> {
    if tu_id == "." {
        return Err(MembershipError::FractionalSentinelTuId {
            path: path.to_path_buf(),
            line,
            read_id: read_id.to_owned(),
        });
    }
    if tu_id.contains(',') {
        return Err(MembershipError::ReservedFractionalTuId {
            path: path.to_path_buf(),
            line,
            read_id: read_id.to_owned(),
            tu_id: tu_id.to_owned(),
        });
    }
    Ok(())
}

fn parse_primary_weight(
    path: &Path,
    line: usize,
    read_id: &str,
    value: &str,
) -> Result<f64, MembershipError> {
    if value.is_empty() || value == "." {
        return Ok(0.0);
    }
    let weight = value
        .parse::<f64>()
        .map_err(|_| MembershipError::InvalidPrimaryWeight {
            path: path.to_path_buf(),
            line,
            read_id: read_id.to_owned(),
            value: value.to_owned(),
        })?;
    if !weight.is_finite() || !(0.0..=1.0).contains(&weight) {
        return Err(MembershipError::InvalidPrimaryWeight {
            path: path.to_path_buf(),
            line,
            read_id: read_id.to_owned(),
            value: value.to_owned(),
        });
    }
    Ok(weight)
}

#[allow(clippy::too_many_arguments)]
fn validate_v2_invariants(
    path: &Path,
    line: usize,
    read_id: &str,
    status: MembershipStatus,
    hard_tu_id: Option<&str>,
    best_candidate_tu_id: Option<&str>,
    second_candidate_tu_id: Option<&str>,
    fractional_assignments: &[WeightedTu],
    primary_weight: f64,
) -> Result<(), MembershipError> {
    let inconsistent = |reason: String| MembershipError::InconsistentV2 {
        path: path.to_path_buf(),
        line,
        read_id: read_id.to_owned(),
        reason,
    };

    if best_candidate_tu_id.is_some() && best_candidate_tu_id == second_candidate_tu_id {
        return Err(inconsistent(
            "best and second candidate TU IDs must differ".to_owned(),
        ));
    }

    match status {
        MembershipStatus::Unique | MembershipStatus::Partial => {
            let hard_tu_id = hard_tu_id.expect("hard assignment checked before invariants");
            let Some(best_candidate_tu_id) = best_candidate_tu_id else {
                return Err(inconsistent(format!(
                    "{} assignment requires a best candidate TU ID",
                    status.as_str()
                )));
            };
            if hard_tu_id != best_candidate_tu_id {
                return Err(inconsistent(format!(
                    "hard TU ID {hard_tu_id:?} must equal best candidate {best_candidate_tu_id:?}"
                )));
            }
            if primary_weight != 1.0 {
                return Err(inconsistent(format!(
                    "{} primary_weight must be 1, got {primary_weight}",
                    status.as_str()
                )));
            }
            if fractional_assignments.len() != 1 || fractional_assignments[0].tu_id != hard_tu_id {
                return Err(inconsistent(format!(
                    "{} fractional assignments must contain exactly hard TU {hard_tu_id:?}",
                    status.as_str()
                )));
            }
            if fractional_assignments[0].weight != 1.0 {
                return Err(inconsistent(format!(
                    "hard-TU fractional weight must be 1, got {}",
                    fractional_assignments[0].weight
                )));
            }
        }
        MembershipStatus::Ambiguous => {
            let Some(best_candidate_tu_id) = best_candidate_tu_id else {
                return Err(inconsistent(
                    "ambiguous assignment requires a best candidate TU ID".to_owned(),
                ));
            };
            let Some(second_candidate_tu_id) = second_candidate_tu_id else {
                return Err(inconsistent(
                    "ambiguous assignment requires a second candidate TU ID".to_owned(),
                ));
            };
            if fractional_assignments.is_empty() {
                if primary_weight != 0.0 {
                    return Err(inconsistent(format!(
                        "primary_weight must be 0 without fractional assignments, got {primary_weight}"
                    )));
                }
            } else {
                let best_weight = fractional_assignments
                    .iter()
                    .find(|assignment| assignment.tu_id == best_candidate_tu_id)
                    .map(|assignment| assignment.weight)
                    .ok_or_else(|| {
                        inconsistent(format!(
                            "best candidate {best_candidate_tu_id:?} is absent from fractional assignments"
                        ))
                    })?;
                if !fractional_assignments
                    .iter()
                    .any(|assignment| assignment.tu_id == second_candidate_tu_id)
                {
                    return Err(inconsistent(format!(
                        "second candidate {second_candidate_tu_id:?} is absent from fractional assignments"
                    )));
                }
                if primary_weight != best_weight {
                    return Err(inconsistent(format!(
                        "primary_weight {primary_weight} must equal best-candidate fractional weight {best_weight}"
                    )));
                }
            }
        }
        MembershipStatus::Unassigned => {
            if best_candidate_tu_id.is_some() || second_candidate_tu_id.is_some() {
                return Err(inconsistent(
                    "unassigned status cannot have candidate TU IDs".to_owned(),
                ));
            }
            if !fractional_assignments.is_empty() {
                return Err(MembershipError::UnassignedFractions {
                    path: path.to_path_buf(),
                    line,
                    read_id: read_id.to_owned(),
                });
            }
            if primary_weight != 0.0 {
                return Err(inconsistent(format!(
                    "unassigned primary_weight must be 0, got {primary_weight}"
                )));
            }
        }
        MembershipStatus::Legacy => unreachable!("v2 rows cannot have legacy status"),
    }
    Ok(())
}

fn parse_fractional_assignments(
    path: &Path,
    value: &str,
    line: usize,
    read_id: &str,
) -> Result<Vec<WeightedTu>, MembershipError> {
    if value.is_empty() || value == "." {
        return Ok(Vec::new());
    }
    let mut assignments = Vec::new();
    let mut seen = HashSet::new();
    for item in value.split(',') {
        let Some((tu_id, weight)) = item.rsplit_once(':') else {
            return Err(MembershipError::InvalidFraction {
                path: path.to_path_buf(),
                line,
                read_id: read_id.to_owned(),
                value: item.to_owned(),
            });
        };
        validate_fractional_tu_id(path, line, read_id, tu_id)?;
        if tu_id.is_empty() || !seen.insert(tu_id.to_owned()) {
            return Err(MembershipError::DuplicateFractionalTu {
                path: path.to_path_buf(),
                line,
                read_id: read_id.to_owned(),
                tu_id: tu_id.to_owned(),
            });
        }
        let weight =
            weight
                .parse::<f64>()
                .map_err(|_| MembershipError::InvalidFractionalWeight {
                    path: path.to_path_buf(),
                    line,
                    read_id: read_id.to_owned(),
                    value: weight.to_owned(),
                })?;
        if !weight.is_finite() || weight <= 0.0 || weight > 1.0 {
            return Err(MembershipError::InvalidFractionalWeight {
                path: path.to_path_buf(),
                line,
                read_id: read_id.to_owned(),
                value: weight.to_string(),
            });
        }
        assignments.push(WeightedTu {
            tu_id: tu_id.to_owned(),
            weight,
        });
    }
    let sum = assignments.iter().map(|assignment| assignment.weight).sum();
    if sum != 1.0 {
        return Err(MembershipError::FractionalWeightSum {
            path: path.to_path_buf(),
            line,
            read_id: read_id.to_owned(),
            sum,
        });
    }
    Ok(assignments)
}

pub(crate) fn membership_metadata(
    ambiguity_margin: f64,
    fractional_assignment: bool,
) -> Vec<String> {
    vec![
        "#trackclustertu_membership_schema=v2".to_owned(),
        format!("#columns={}", V2_COLUMNS.join("\t")),
        "#legacy_score_columns=score1:span_jaccard,score2:overlap_over_longer".to_owned(),
        "#assignment_score=overlap_over_longer".to_owned(),
        "#ambiguity_rule=best_assignment_score-second_assignment_score<=ambiguity_margin"
            .to_owned(),
        format!("#ambiguity_margin={ambiguity_margin}"),
        format!("#fractional_assignment={fractional_assignment}"),
        "#lexical_tu_ids_are_ordering_only=true".to_owned(),
    ]
}

pub(crate) fn write_membership<W: Write>(
    writer: &mut DelimitedWriter<W>,
    record: &MembershipRecord,
) -> Result<(), MembershipError> {
    if record
        .fractional_assignments
        .iter()
        .any(|assignment| assignment.tu_id == ".")
    {
        return Err(MembershipError::FractionalSentinelTuIdWrite {
            read_id: record.read_id.clone(),
        });
    }
    if let Some(assignment) = record
        .fractional_assignments
        .iter()
        .find(|assignment| assignment.tu_id.contains(','))
    {
        return Err(MembershipError::ReservedFractionalTuIdWrite {
            read_id: record.read_id.clone(),
            tu_id: assignment.tu_id.clone(),
        });
    }
    let hard_tu_id = record.hard_tu_id.as_deref().unwrap_or(".");
    if let Some(v2) = &record.v2 {
        if !v2.complete {
            writer.write_record([
                record.read_id.as_str(),
                hard_tu_id,
                record.score1.as_str(),
                record.score2.as_str(),
                "v2",
                record.status.as_str(),
            ])?;
            return Ok(());
        }
        let fractional = if record.fractional_assignments.is_empty() {
            ".".to_owned()
        } else {
            record
                .fractional_assignments
                .iter()
                .map(|assignment| format!("{}:{}", assignment.tu_id, assignment.weight))
                .collect::<Vec<_>>()
                .join(",")
        };
        writer.write_record([
            record.read_id.as_str(),
            hard_tu_id,
            record.score1.as_str(),
            record.score2.as_str(),
            "v2",
            record.status.as_str(),
            v2.best_candidate_tu_id.as_deref().unwrap_or("."),
            v2.second_candidate_tu_id.as_deref().unwrap_or("."),
            v2.second_score1.as_str(),
            v2.second_score2.as_str(),
            v2.best_five_prime_delta_bp.as_str(),
            v2.best_three_prime_delta_bp.as_str(),
            v2.second_five_prime_delta_bp.as_str(),
            v2.second_three_prime_delta_bp.as_str(),
            v2.best_assignment_score.as_str(),
            v2.second_assignment_score.as_str(),
            v2.score_margin.as_str(),
            v2.primary_weight.as_str(),
            fractional.as_str(),
            if v2.full_length_evidence {
                "true"
            } else {
                "false"
            },
        ])?;
    } else {
        writer.write_record([
            record.read_id.as_str(),
            hard_tu_id,
            record.score1.as_str(),
            record.score2.as_str(),
        ])?;
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    fn parse_text(input: &str) -> Result<Vec<MembershipRecord>, MembershipError> {
        let reader = DelimitedReader::from_reader(
            Path::new("membership.tsv"),
            Delimiter::Tab,
            input.as_bytes(),
        );
        parse_memberships(reader, true)
    }

    #[test]
    fn metadata_scan_ignores_hash_lines_inside_quoted_multiline_fields() {
        let input = concat!(
            "#trackclustertu_membership_schema=v2\n",
            "\"read\n",
            "#not_metadata\"\tTU000001\t1\t1\n",
            "  #ambiguity_margin=0.02\r\n",
        );
        let mut reader = std::io::Cursor::new(input.as_bytes());
        let metadata = scan_membership_metadata(&mut reader).unwrap();
        assert_eq!(
            metadata,
            [
                "#trackclustertu_membership_schema=v2",
                "#ambiguity_margin=0.02",
            ]
        );
    }

    fn full_v2_row(
        read_id: &str,
        hard_tu_id: &str,
        status: &str,
        best: &str,
        second: &str,
        primary_weight: &str,
        fractions: &str,
    ) -> String {
        format!(
            "{read_id}\t{hard_tu_id}\t.\t.\tv2\t{status}\t{best}\t{second}\t.\t.\t.\t.\t.\t.\t.\t.\t.\t{primary_weight}\t{fractions}\tfalse\n"
        )
    }

    #[test]
    fn v1_and_v2_memberships_share_quoted_parser_and_no_hard_semantics() {
        let input = concat!(
            "\"read\tlegacy\"\t\"TU\"\"legacy\"\t1\t1\n",
            "\"read\nambiguous\"\t.\t.\t.\tv2\tambiguous\tTU_A\tTU_B\t.\t.\t.\t.\t.\t.\t.\t.\t0\t0.5\tTU_A:0.5,TU_B:0.5\tfalse\n",
            "unassigned\t.\t.\t.\tv2\tunassigned\t.\t.\t.\t.\t.\t.\t.\t.\t.\t.\t.\t0\t.\tfalse\n",
        );
        let reader = DelimitedReader::from_reader(
            Path::new("membership.tsv"),
            Delimiter::Tab,
            input.as_bytes(),
        );
        let records = parse_memberships(reader, true).unwrap();
        assert_eq!(records.len(), 3);
        assert_eq!(records[0].read_id, "read\tlegacy");
        assert_eq!(records[0].hard_tu_id.as_deref(), Some("TU\"legacy"));
        assert_eq!(records[1].read_id, "read\nambiguous");
        assert_eq!(records[1].hard_tu_id, None);
        assert_eq!(records[1].contributions().len(), 2);
        assert_eq!(records[2].status, MembershipStatus::Unassigned);
    }

    #[test]
    fn membership_writer_round_trips_quoted_v2_fields_and_fractions() {
        let path = Path::new("membership.tsv");
        let tu_id = " TU:colon\tquote\"line\nbreak ".to_owned();
        let record = MembershipRecord {
            read_id: " read,tab\tquote\"line\r\nbreak ".to_owned(),
            hard_tu_id: Some(tu_id.clone()),
            score1: "1.000000".to_owned(),
            score2: "1.000000".to_owned(),
            status: MembershipStatus::Unique,
            fractional_assignments: vec![WeightedTu {
                tu_id: tu_id.clone(),
                weight: 1.0,
            }],
            v2: Some(V2MembershipFields {
                complete: true,
                best_candidate_tu_id: Some(tu_id.clone()),
                second_candidate_tu_id: None,
                second_score1: ".".to_owned(),
                second_score2: ".".to_owned(),
                best_five_prime_delta_bp: "0".to_owned(),
                best_three_prime_delta_bp: "0".to_owned(),
                second_five_prime_delta_bp: ".".to_owned(),
                second_three_prime_delta_bp: ".".to_owned(),
                best_assignment_score: "1.000000".to_owned(),
                second_assignment_score: ".".to_owned(),
                score_margin: ".".to_owned(),
                primary_weight: "1".to_owned(),
                full_length_evidence: true,
            }),
            line_number: 0,
        };
        let mut bytes = Vec::new();
        {
            let mut writer = DelimitedWriter::from_writer(path, Delimiter::Tab, &mut bytes);
            write_membership(&mut writer, &record).unwrap();
            writer.flush().unwrap();
        }
        let reader = DelimitedReader::from_reader(path, Delimiter::Tab, bytes.as_slice());
        let reparsed = parse_memberships(reader, true).unwrap();
        assert_eq!(reparsed.len(), 1);
        assert_eq!(reparsed[0].read_id, record.read_id);
        assert_eq!(reparsed[0].hard_tu_id, record.hard_tu_id);
        assert_eq!(
            reparsed[0].fractional_assignments,
            record.fractional_assignments
        );
        assert!(reparsed[0]
            .v2
            .as_ref()
            .is_some_and(|v2| v2.full_length_evidence));
    }

    #[test]
    fn exact_identifier_bytes_and_supported_short_v2_rows_are_preserved() {
        let rows = parse_text(concat!(
            " read legacy \t TU legacy \t1\t1\n",
            "short\tTU_SHORT\t.\t.\tv2\tunique\n",
            "legacy-dot\tTU_DOT\t.\t.\t.\n",
        ))
        .unwrap();
        assert_eq!(rows[0].read_id, " read legacy ");
        assert_eq!(rows[0].hard_tu_id.as_deref(), Some(" TU legacy "));
        assert_eq!(rows[1].hard_tu_id.as_deref(), Some("TU_SHORT"));
        assert_eq!(rows[1].status, MembershipStatus::Unique);
        assert_eq!(rows[2].status, MembershipStatus::Legacy);
    }

    #[test]
    fn unsupported_schema_versions_are_not_silently_legacy() {
        let error = parse_text("read\tTU\t.\t.\tv3\tunique\n").unwrap_err();
        assert!(matches!(
            error,
            MembershipError::UnsupportedSchema { schema, .. } if schema == "v3"
        ));
    }

    #[test]
    fn fractional_comma_ids_are_rejected_on_read_and_write() {
        let input = full_v2_row("read", "TU,A", "unique", "TU,A", ".", "1", "TU,A:1");
        let error = parse_text(&input).unwrap_err();
        assert!(matches!(
            error,
            MembershipError::ReservedFractionalTuId { tu_id, .. } if tu_id == "TU,A"
        ));

        let path = Path::new("membership.tsv");
        let record = MembershipRecord {
            read_id: "read".to_owned(),
            hard_tu_id: Some("TU,A".to_owned()),
            score1: ".".to_owned(),
            score2: ".".to_owned(),
            status: MembershipStatus::Unique,
            fractional_assignments: vec![WeightedTu {
                tu_id: "TU,A".to_owned(),
                weight: 1.0,
            }],
            v2: Some(V2MembershipFields {
                complete: true,
                best_candidate_tu_id: Some("TU,A".to_owned()),
                second_candidate_tu_id: None,
                second_score1: ".".to_owned(),
                second_score2: ".".to_owned(),
                best_five_prime_delta_bp: ".".to_owned(),
                best_three_prime_delta_bp: ".".to_owned(),
                second_five_prime_delta_bp: ".".to_owned(),
                second_three_prime_delta_bp: ".".to_owned(),
                best_assignment_score: ".".to_owned(),
                second_assignment_score: ".".to_owned(),
                score_margin: ".".to_owned(),
                primary_weight: "1".to_owned(),
                full_length_evidence: false,
            }),
            line_number: 0,
        };
        let mut bytes = Vec::new();
        let mut writer = DelimitedWriter::from_writer(path, Delimiter::Tab, &mut bytes);
        let error = write_membership(&mut writer, &record).unwrap_err();
        assert!(matches!(
            error,
            MembershipError::ReservedFractionalTuIdWrite { tu_id, .. } if tu_id == "TU,A"
        ));
        drop(writer);
        assert!(bytes.is_empty());

        let dot_fraction = full_v2_row(
            "dot-fraction",
            ".",
            "ambiguous",
            "TU_A",
            "TU_B",
            "0.5",
            ".:0.5,TU_B:0.5",
        );
        assert!(matches!(
            parse_text(&dot_fraction).unwrap_err(),
            MembershipError::FractionalSentinelTuId { .. }
        ));
    }

    #[test]
    fn full_v2_rows_reject_contradictory_audit_fields() {
        let cases = [
            (
                full_v2_row("hard-best", "TU_A", "unique", "TU_B", ".", "1", "."),
                "must equal best candidate",
            ),
            (
                full_v2_row("unique-weight", "TU_A", "unique", "TU_A", ".", "0", "."),
                "primary_weight must be 1",
            ),
            (
                full_v2_row(
                    "unique-fractions",
                    "TU_A",
                    "unique",
                    "TU_A",
                    "TU_B",
                    "1",
                    "TU_A:0.5,TU_B:0.5",
                ),
                "must contain exactly hard TU",
            ),
            (
                full_v2_row(
                    "unique-missing-fraction",
                    "TU_A",
                    "unique",
                    "TU_A",
                    ".",
                    "1",
                    ".",
                ),
                "must contain exactly hard TU",
            ),
            (
                full_v2_row(
                    "ambiguous-weight",
                    ".",
                    "ambiguous",
                    "TU_A",
                    "TU_B",
                    "0.4",
                    "TU_A:0.5,TU_B:0.5",
                ),
                "must equal best-candidate fractional weight",
            ),
            (
                full_v2_row(
                    "ambiguous-second",
                    ".",
                    "ambiguous",
                    "TU_A",
                    "TU_B",
                    "0.5",
                    "TU_A:0.5,TU_C:0.5",
                ),
                "second candidate",
            ),
            (
                full_v2_row("ambiguous-candidates", ".", "ambiguous", ".", ".", "0", "."),
                "requires a best candidate",
            ),
            (
                full_v2_row(
                    "unassigned-candidate",
                    ".",
                    "unassigned",
                    "TU_A",
                    ".",
                    "0",
                    ".",
                ),
                "cannot have candidate",
            ),
            (
                full_v2_row("unassigned-weight", ".", "unassigned", ".", ".", "0.5", "."),
                "primary_weight must be 0",
            ),
        ];

        for (row, expected) in cases {
            let error = parse_text(&row).unwrap_err();
            let message = error.to_string();
            assert!(message.contains(expected), "{message}");
            assert!(message.contains("membership.tsv\":1"), "{message}");
        }
    }

    #[test]
    fn ambiguous_fractional_rows_allow_more_than_two_consistent_candidates() {
        let rows = [
            full_v2_row(
                "fractional",
                ".",
                "ambiguous",
                "TU_A",
                "TU_B",
                "0.5",
                "TU_A:0.5,TU_B:0.25,TU_C:0.25",
            ),
            full_v2_row(
                "ambiguous-no-fractions",
                ".",
                "ambiguous",
                "TU_A",
                "TU_B",
                "0",
                ".",
            ),
            full_v2_row("unique", "TU_A", "unique", "TU_A", ".", "1", "TU_A:1"),
            full_v2_row("unassigned", ".", "unassigned", ".", ".", "0", "."),
        ]
        .concat();
        let records = parse_text(&rows).unwrap();
        assert_eq!(records[0].fractional_assignments.len(), 3);
        assert_eq!(records[0].contributions().len(), 3);
        assert_eq!(records.len(), 4);
    }

    #[test]
    fn partially_present_v2_audit_columns_are_rejected() {
        let error = parse_text("read\tTU_A\t.\t.\tv2\tunique\tTU_A\n").unwrap_err();
        assert!(matches!(
            error,
            MembershipError::TooFewColumns {
                expected: 20,
                got: 7,
                ..
            }
        ));
    }
}
