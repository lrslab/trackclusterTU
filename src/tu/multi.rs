use std::collections::{HashMap, HashSet};
use std::path::{Path, PathBuf};

use thiserror::Error;

use crate::io::delimited::{DelimitedError, DelimitedReader, Delimiter};

const POOLED_ID_SEPARATOR: &str = "::";

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct SampleManifestRecord {
    pub sample: String,
    pub reads: PathBuf,
    pub group: Option<String>,
    pub evidence: Option<PathBuf>,
}

#[derive(Clone, Debug, PartialEq)]
pub struct MultiSampleCounts {
    tu_ids: Vec<String>,
    samples: Vec<String>,
    /// Non-ambiguous assignments retained as the legacy integer matrix.
    counts: Vec<Vec<u64>>,
    unique_counts: Vec<Vec<u64>>,
    full_length_evidence_counts: Vec<Vec<u64>>,
    fractional_counts: Vec<Vec<f64>>,
    groups: Vec<String>,
    /// Non-ambiguous assignments retained as the legacy integer group matrix.
    group_counts: Vec<Vec<u64>>,
    group_unique_counts: Vec<Vec<u64>>,
    group_full_length_evidence_counts: Vec<Vec<u64>>,
    group_fractional_counts: Vec<Vec<f64>>,
}

impl MultiSampleCounts {
    /// TU identifiers aligned one-to-one with every matrix row.
    pub fn tu_ids(&self) -> &[String] {
        &self.tu_ids
    }

    /// Sample identifiers aligned one-to-one with each sample-matrix column.
    pub fn samples(&self) -> &[String] {
        &self.samples
    }

    /// Group identifiers aligned one-to-one with each group-matrix column.
    pub fn groups(&self) -> &[String] {
        &self.groups
    }

    /// Total hard-assignment counts by TU row and sample column.
    pub fn counts(&self) -> &[Vec<u64>] {
        &self.counts
    }

    /// Unique hard-assignment counts by TU row and sample column.
    pub fn unique_counts(&self) -> &[Vec<u64>] {
        &self.unique_counts
    }

    /// Full-length-evidence counts by TU row and sample column.
    pub fn full_length_evidence_counts(&self) -> &[Vec<u64>] {
        &self.full_length_evidence_counts
    }

    /// Fractional assignment counts by TU row and sample column.
    pub fn fractional_counts(&self) -> &[Vec<f64>] {
        &self.fractional_counts
    }

    /// Total hard-assignment counts by TU row and group column.
    pub fn group_counts(&self) -> &[Vec<u64>] {
        &self.group_counts
    }

    /// Unique hard-assignment counts by TU row and group column.
    pub fn group_unique_counts(&self) -> &[Vec<u64>] {
        &self.group_unique_counts
    }

    /// Full-length-evidence counts by TU row and group column.
    pub fn group_full_length_evidence_counts(&self) -> &[Vec<u64>] {
        &self.group_full_length_evidence_counts
    }

    /// Fractional assignment counts by TU row and group column.
    pub fn group_fractional_counts(&self) -> &[Vec<f64>] {
        &self.group_fractional_counts
    }

    pub fn has_groups(&self) -> bool {
        !self.groups.is_empty()
    }

    /// Consume the result and retain only TU rows meeting `minimum` total hard assignments.
    ///
    /// Every sample/group metric matrix is filtered with the same row mask, preserving the
    /// alignment that callers could previously break by independently rebuilding public vectors.
    pub fn filter_min_total_count(mut self, minimum: u64) -> Self {
        let keep: Vec<bool> = self
            .counts
            .iter()
            .map(|row| row.iter().sum::<u64>() >= minimum)
            .collect();

        self.tu_ids = retain_rows(self.tu_ids, &keep);
        self.counts = retain_rows(self.counts, &keep);
        self.unique_counts = retain_rows(self.unique_counts, &keep);
        self.full_length_evidence_counts = retain_rows(self.full_length_evidence_counts, &keep);
        self.fractional_counts = retain_rows(self.fractional_counts, &keep);
        self.group_counts = retain_rows(self.group_counts, &keep);
        self.group_unique_counts = retain_rows(self.group_unique_counts, &keep);
        self.group_full_length_evidence_counts =
            retain_rows(self.group_full_length_evidence_counts, &keep);
        self.group_fractional_counts = retain_rows(self.group_fractional_counts, &keep);
        self.assert_aligned();
        self
    }

    fn assert_aligned(&self) {
        let row_count = self.tu_ids.len();
        for matrix_rows in [
            self.counts.len(),
            self.unique_counts.len(),
            self.full_length_evidence_counts.len(),
            self.fractional_counts.len(),
            self.group_counts.len(),
            self.group_unique_counts.len(),
            self.group_full_length_evidence_counts.len(),
            self.group_fractional_counts.len(),
        ] {
            assert_eq!(
                matrix_rows, row_count,
                "TU count matrices must stay row-aligned"
            );
        }
        for row in self
            .counts
            .iter()
            .chain(&self.unique_counts)
            .chain(&self.full_length_evidence_counts)
        {
            assert_eq!(
                row.len(),
                self.samples.len(),
                "sample count rows must match the sample columns"
            );
        }
        for row in &self.fractional_counts {
            assert_eq!(
                row.len(),
                self.samples.len(),
                "sample count rows must match the sample columns"
            );
        }
        for row in self
            .group_counts
            .iter()
            .chain(&self.group_unique_counts)
            .chain(&self.group_full_length_evidence_counts)
        {
            assert_eq!(
                row.len(),
                self.groups.len(),
                "group count rows must match the group columns"
            );
        }
        for row in &self.group_fractional_counts {
            assert_eq!(
                row.len(),
                self.groups.len(),
                "group count rows must match the group columns"
            );
        }
    }
}

fn retain_rows<T>(rows: Vec<T>, keep: &[bool]) -> Vec<T> {
    assert_eq!(rows.len(), keep.len(), "row mask must match matrix rows");
    rows.into_iter()
        .zip(keep.iter().copied())
        .filter_map(|(row, retain)| retain.then_some(row))
        .collect()
}

#[derive(Error, Debug)]
pub enum MultiSampleError {
    #[error(transparent)]
    Delimited(#[from] DelimitedError),

    #[error("manifest {path:?} is empty")]
    EmptyManifest { path: PathBuf },

    #[error("{path:?}:{line}: manifest is missing required column {column:?}")]
    MissingColumn {
        path: PathBuf,
        line: usize,
        column: &'static str,
    },

    #[error("{path:?}:{line}: duplicate manifest column {column:?}")]
    DuplicateColumn {
        path: PathBuf,
        line: usize,
        column: String,
    },

    #[error("{path:?}:{line}: expected at least {expected} columns, got {got}")]
    TooFewColumns {
        path: PathBuf,
        line: usize,
        expected: usize,
        got: usize,
    },

    #[error("{path:?}:{line}: sample name cannot be empty")]
    EmptySampleName { path: PathBuf, line: usize },

    #[error("{path:?}:{line}: reads path cannot be empty for sample {sample:?}")]
    EmptyReadsPath {
        path: PathBuf,
        line: usize,
        sample: String,
    },

    #[error("{path:?}:{line}: sample name {sample:?} contains reserved separator \"{separator}\"")]
    ReservedSeparatorInSample {
        path: PathBuf,
        line: usize,
        sample: String,
        separator: &'static str,
    },

    #[error("{path:?}:{line}: duplicate sample name {sample:?}")]
    DuplicateSample {
        path: PathBuf,
        line: usize,
        sample: String,
    },

    #[error("pooled read id {read_id:?} is missing the <sample>{separator}<read_id> tag")]
    InvalidPooledReadId {
        read_id: String,
        separator: &'static str,
    },

    #[error("pooled read id {read_id:?} references unknown sample {sample:?}")]
    UnknownSample { read_id: String, sample: String },

    #[error("pooled membership references unknown TU {tu_id:?}")]
    UnknownTuId { tu_id: String },

    #[error("fractional assignment weight must be finite and between 0 and 1, got {weight}")]
    InvalidAssignmentWeight { weight: f64 },
}

pub fn parse_manifest(path: &Path) -> Result<Vec<SampleManifestRecord>, MultiSampleError> {
    let path = path.to_path_buf();
    let mut reader = DelimitedReader::open(&path, Delimiter::Tab)?;
    let base_dir = path.parent().unwrap_or(Path::new("."));

    let mut header_line: Option<usize> = None;
    let mut sample_col: Option<usize> = None;
    let mut reads_col: Option<usize> = None;
    let mut group_col: Option<usize> = None;
    let mut evidence_col: Option<usize> = None;

    let mut seen_samples: HashSet<String> = HashSet::new();
    let mut records: Vec<SampleManifestRecord> = Vec::new();

    for table_record in reader.records() {
        let table_record = table_record?;
        let line_number = table_record.line_number() as usize;
        let fields = table_record.fields();
        if header_line.is_none() {
            let mut seen_columns: HashSet<String> = HashSet::new();
            for (idx, raw_column) in fields.iter().enumerate() {
                let column = raw_column.trim().to_ascii_lowercase();
                if !seen_columns.insert(column.clone()) {
                    return Err(MultiSampleError::DuplicateColumn {
                        path: path.clone(),
                        line: line_number,
                        column,
                    });
                }

                match column.as_str() {
                    "sample" => sample_col = Some(idx),
                    "reads" => reads_col = Some(idx),
                    "group" => group_col = Some(idx),
                    "evidence" | "evidence_tsv" => {
                        if evidence_col.replace(idx).is_some() {
                            return Err(MultiSampleError::DuplicateColumn {
                                path: path.clone(),
                                line: line_number,
                                column: "evidence/evidence_tsv".to_owned(),
                            });
                        }
                    }
                    _ => {}
                }
            }

            header_line = Some(line_number);
            if sample_col.is_none() {
                return Err(MultiSampleError::MissingColumn {
                    path: path.clone(),
                    line: line_number,
                    column: "sample",
                });
            }
            if reads_col.is_none() {
                return Err(MultiSampleError::MissingColumn {
                    path: path.clone(),
                    line: line_number,
                    column: "reads",
                });
            }
            continue;
        }

        let sample_idx = sample_col.expect("checked above");
        let reads_idx = reads_col.expect("checked above");
        let required_columns = sample_idx.max(reads_idx).max(group_col.unwrap_or(0)) + 1;
        if fields.len() < required_columns {
            return Err(MultiSampleError::TooFewColumns {
                path: path.clone(),
                line: line_number,
                expected: required_columns,
                got: fields.len(),
            });
        }

        let sample = fields.get(sample_idx).unwrap_or("").trim().to_owned();
        if sample.is_empty() {
            return Err(MultiSampleError::EmptySampleName {
                path: path.clone(),
                line: line_number,
            });
        }
        if sample.contains(POOLED_ID_SEPARATOR) {
            return Err(MultiSampleError::ReservedSeparatorInSample {
                path: path.clone(),
                line: line_number,
                sample,
                separator: POOLED_ID_SEPARATOR,
            });
        }
        if !seen_samples.insert(sample.clone()) {
            return Err(MultiSampleError::DuplicateSample {
                path: path.clone(),
                line: line_number,
                sample,
            });
        }

        let reads_value = fields.get(reads_idx).unwrap_or("").trim();
        if reads_value.is_empty() {
            return Err(MultiSampleError::EmptyReadsPath {
                path: path.clone(),
                line: line_number,
                sample,
            });
        }

        let reads = {
            let reads_path = PathBuf::from(reads_value);
            if reads_path.is_absolute() {
                reads_path
            } else {
                base_dir.join(reads_path)
            }
        };

        let group = group_col
            .and_then(|idx| fields.get(idx))
            .map(|value| value.trim())
            .filter(|value| !value.is_empty())
            .map(str::to_owned);
        let evidence = evidence_col
            .and_then(|idx| fields.get(idx))
            .map(|value| value.trim())
            .filter(|value| !value.is_empty() && *value != ".")
            .map(|value| {
                let evidence_path = PathBuf::from(value);
                if evidence_path.is_absolute() {
                    evidence_path
                } else {
                    base_dir.join(evidence_path)
                }
            });

        records.push(SampleManifestRecord {
            sample,
            reads,
            group,
            evidence,
        });
    }

    if header_line.is_none() {
        return Err(MultiSampleError::EmptyManifest { path });
    }

    Ok(records)
}

pub fn pooled_read_id(sample: &str, read_id: &str) -> String {
    format!("{sample}{POOLED_ID_SEPARATOR}{read_id}")
}

pub fn split_pooled_read_id(read_id: &str) -> Option<(&str, &str)> {
    read_id.split_once(POOLED_ID_SEPARATOR)
}

#[derive(Clone, Debug)]
pub struct MultiSampleCounter {
    tu_ids: Vec<String>,
    samples: Vec<String>,
    counts: Vec<Vec<u64>>,
    unique_counts: Vec<Vec<u64>>,
    full_length_evidence_counts: Vec<Vec<u64>>,
    fractional_counts: Vec<Vec<f64>>,
    groups: Vec<String>,
    group_counts: Vec<Vec<u64>>,
    group_unique_counts: Vec<Vec<u64>>,
    group_full_length_evidence_counts: Vec<Vec<u64>>,
    group_fractional_counts: Vec<Vec<f64>>,
    tu_to_idx: HashMap<String, usize>,
    sample_to_idx: HashMap<String, usize>,
    sample_to_group_idx: Vec<Option<usize>>,
}

impl MultiSampleCounter {
    pub fn new(manifest: &[SampleManifestRecord], tu_ids: Vec<String>) -> Self {
        let mut samples: Vec<String> = Vec::with_capacity(manifest.len());
        let mut sample_to_idx: HashMap<String, usize> = HashMap::with_capacity(manifest.len());
        for (idx, record) in manifest.iter().enumerate() {
            sample_to_idx.insert(record.sample.clone(), idx);
            samples.push(record.sample.clone());
        }

        let mut groups: Vec<String> = Vec::new();
        let mut group_to_idx: HashMap<String, usize> = HashMap::new();
        let mut sample_to_group_idx: Vec<Option<usize>> = Vec::with_capacity(manifest.len());
        for record in manifest {
            let maybe_group_idx = match record.group.as_ref() {
                Some(group) => {
                    let idx = match group_to_idx.get(group) {
                        Some(&idx) => idx,
                        None => {
                            let idx = groups.len();
                            groups.push(group.clone());
                            group_to_idx.insert(group.clone(), idx);
                            idx
                        }
                    };
                    Some(idx)
                }
                None => None,
            };
            sample_to_group_idx.push(maybe_group_idx);
        }

        let tu_to_idx = tu_ids
            .iter()
            .enumerate()
            .map(|(idx, tu_id)| (tu_id.clone(), idx))
            .collect();

        let counts = vec![vec![0; samples.len()]; tu_ids.len()];
        let unique_counts = vec![vec![0; samples.len()]; tu_ids.len()];
        let full_length_evidence_counts = vec![vec![0; samples.len()]; tu_ids.len()];
        let fractional_counts = vec![vec![0.0; samples.len()]; tu_ids.len()];
        let group_counts = vec![vec![0; groups.len()]; tu_ids.len()];
        let group_unique_counts = vec![vec![0; groups.len()]; tu_ids.len()];
        let group_full_length_evidence_counts = vec![vec![0; groups.len()]; tu_ids.len()];
        let group_fractional_counts = vec![vec![0.0; groups.len()]; tu_ids.len()];

        Self {
            tu_ids,
            samples,
            counts,
            unique_counts,
            full_length_evidence_counts,
            fractional_counts,
            groups,
            group_counts,
            group_unique_counts,
            group_full_length_evidence_counts,
            group_fractional_counts,
            tu_to_idx,
            sample_to_idx,
            sample_to_group_idx,
        }
    }

    #[cfg(test)]
    fn add_assignment(
        &mut self,
        pooled_read_id: &str,
        tu_id: &str,
    ) -> Result<(), MultiSampleError> {
        self.add_assignment_metrics(pooled_read_id, tu_id, false, false, true, 1.0)
    }

    /// Add one TU contribution with explicit assignment semantics.
    ///
    /// Legacy/recount inputs can set only `total_assignment`; clustering additionally records
    /// unique and evidence-backed subsets. `fractional_weight` is independent of the integer
    /// columns so ambiguous reads can contribute fractions without being counted as hard calls.
    pub fn add_assignment_metrics(
        &mut self,
        pooled_read_id: &str,
        tu_id: &str,
        unique_assignment: bool,
        full_length_evidence: bool,
        total_assignment: bool,
        fractional_weight: f64,
    ) -> Result<(), MultiSampleError> {
        if !fractional_weight.is_finite() || !(0.0..=1.0).contains(&fractional_weight) {
            return Err(MultiSampleError::InvalidAssignmentWeight {
                weight: fractional_weight,
            });
        }
        let (sample, _) = split_pooled_read_id(pooled_read_id).ok_or_else(|| {
            MultiSampleError::InvalidPooledReadId {
                read_id: pooled_read_id.to_owned(),
                separator: POOLED_ID_SEPARATOR,
            }
        })?;

        let sample_idx =
            *self
                .sample_to_idx
                .get(sample)
                .ok_or_else(|| MultiSampleError::UnknownSample {
                    read_id: pooled_read_id.to_owned(),
                    sample: sample.to_owned(),
                })?;
        let tu_idx = *self
            .tu_to_idx
            .get(tu_id)
            .ok_or_else(|| MultiSampleError::UnknownTuId {
                tu_id: tu_id.to_owned(),
            })?;

        if total_assignment {
            self.counts[tu_idx][sample_idx] += 1;
        }
        if unique_assignment {
            self.unique_counts[tu_idx][sample_idx] += 1;
        }
        if full_length_evidence {
            self.full_length_evidence_counts[tu_idx][sample_idx] += 1;
        }
        self.fractional_counts[tu_idx][sample_idx] += fractional_weight;
        if let Some(group_idx) = self.sample_to_group_idx[sample_idx] {
            if total_assignment {
                self.group_counts[tu_idx][group_idx] += 1;
            }
            if unique_assignment {
                self.group_unique_counts[tu_idx][group_idx] += 1;
            }
            if full_length_evidence {
                self.group_full_length_evidence_counts[tu_idx][group_idx] += 1;
            }
            self.group_fractional_counts[tu_idx][group_idx] += fractional_weight;
        }

        Ok(())
    }

    pub fn finish(self) -> MultiSampleCounts {
        let counts = MultiSampleCounts {
            tu_ids: self.tu_ids,
            samples: self.samples,
            counts: self.counts,
            unique_counts: self.unique_counts,
            full_length_evidence_counts: self.full_length_evidence_counts,
            fractional_counts: self.fractional_counts,
            groups: self.groups,
            group_counts: self.group_counts,
            group_unique_counts: self.group_unique_counts,
            group_full_length_evidence_counts: self.group_full_length_evidence_counts,
            group_fractional_counts: self.group_fractional_counts,
        };
        counts.assert_aligned();
        counts
    }
}

#[cfg(test)]
mod tests {
    use std::fs;
    use std::time::{SystemTime, UNIX_EPOCH};

    use super::*;
    use crate::io::delimited::{DelimitedWriter, Delimiter};

    fn unique_tmp_dir(prefix: &str) -> PathBuf {
        let nanos = SystemTime::now()
            .duration_since(UNIX_EPOCH)
            .expect("clock")
            .as_nanos();
        std::env::temp_dir().join(format!("{prefix}_{nanos}"))
    }

    #[test]
    fn parse_manifest_resolves_relative_paths_and_groups() {
        let tmp = unique_tmp_dir("trackclustertu_manifest_parse");
        fs::create_dir_all(tmp.join("reads")).unwrap();

        let manifest = tmp.join("samples.tsv");
        fs::write(
            &manifest,
            concat!(
                "sample\treads\tgroup\tevidence\n",
                "sampleA\treads/a.bed\tcase\tevidence/a.tsv\n",
                "sampleB\treads/b.bed\tcontrol\t.\n",
            ),
        )
        .unwrap();

        let records = parse_manifest(&manifest).unwrap();
        assert_eq!(
            records,
            vec![
                SampleManifestRecord {
                    sample: "sampleA".to_owned(),
                    reads: tmp.join("reads/a.bed"),
                    group: Some("case".to_owned()),
                    evidence: Some(tmp.join("evidence/a.tsv")),
                },
                SampleManifestRecord {
                    sample: "sampleB".to_owned(),
                    reads: tmp.join("reads/b.bed"),
                    group: Some("control".to_owned()),
                    evidence: None,
                },
            ]
        );

        let _ = fs::remove_dir_all(tmp);
    }

    #[test]
    fn parse_manifest_rejects_duplicate_samples() {
        let tmp = unique_tmp_dir("trackclustertu_manifest_dup");
        fs::create_dir_all(&tmp).unwrap();

        let manifest = tmp.join("samples.tsv");
        fs::write(&manifest, "sample\treads\nsampleA\ta.bed\nsampleA\tb.bed\n").unwrap();

        let error = parse_manifest(&manifest).unwrap_err();
        assert!(matches!(error, MultiSampleError::DuplicateSample { .. }));

        let _ = fs::remove_dir_all(tmp);
    }

    #[test]
    fn parse_manifest_preserves_empty_trailing_group_field() {
        let tmp = unique_tmp_dir("trackclustertu_manifest_empty_group");
        fs::create_dir_all(&tmp).unwrap();

        let manifest = tmp.join("samples.tsv");
        fs::write(
            &manifest,
            "sample\treads\tgroup\nsampleA\ta.bed\t\nsampleB\tb.bed\tcontrol\n",
        )
        .unwrap();

        let records = parse_manifest(&manifest).unwrap();
        assert_eq!(records.len(), 2);
        assert_eq!(records[0].group, None);
        assert_eq!(records[1].group.as_deref(), Some("control"));

        let _ = fs::remove_dir_all(tmp);
    }

    #[test]
    fn manifest_round_trips_quoted_delimiters_quotes_and_line_breaks() {
        let tmp = unique_tmp_dir("trackclustertu_manifest_quoted");
        fs::create_dir_all(&tmp).unwrap();
        let manifest = tmp.join("samples.tsv");
        let sample = "sample,tab\tquote\"line\nbreak";
        let group = "group,tab\tquote\"crlf\r\nvalue";
        let reads = "reads,tab\tquote\"line\nbreak.bed";
        let mut writer = DelimitedWriter::create(&manifest, Delimiter::Tab, &[]).unwrap();
        writer.write_record(["sample", "reads", "group"]).unwrap();
        writer.write_record([sample, reads, group]).unwrap();
        writer.flush().unwrap();

        let records = parse_manifest(&manifest).unwrap();
        assert_eq!(records.len(), 1);
        assert_eq!(records[0].sample, sample);
        assert_eq!(records[0].group.as_deref(), Some(group));
        assert_eq!(records[0].reads, tmp.join(reads));
        let _ = fs::remove_dir_all(tmp);
    }

    #[test]
    fn counter_aggregates_sample_and_group_counts() {
        let manifest = vec![
            SampleManifestRecord {
                sample: "sampleA".to_owned(),
                reads: PathBuf::from("a.bed"),
                group: Some("control".to_owned()),
                evidence: None,
            },
            SampleManifestRecord {
                sample: "sampleB".to_owned(),
                reads: PathBuf::from("b.bed"),
                group: Some("treated".to_owned()),
                evidence: None,
            },
        ];

        let mut counter = MultiSampleCounter::new(
            &manifest,
            vec![
                "TU000001".to_owned(),
                "TU000002".to_owned(),
                "TU000003".to_owned(),
            ],
        );

        counter.add_assignment("sampleA::r1", "TU000001").unwrap();
        counter.add_assignment("sampleA::r2", "TU000001").unwrap();
        counter.add_assignment("sampleB::r3", "TU000001").unwrap();
        counter.add_assignment("sampleB::r4", "TU000003").unwrap();

        let counts = counter.finish();
        assert_eq!(counts.samples(), &["sampleA", "sampleB"]);
        assert_eq!(counts.groups(), &["control", "treated"]);
        assert_eq!(counts.counts()[0], vec![2, 1]);
        assert_eq!(counts.counts()[1], vec![0, 0]);
        assert_eq!(counts.counts()[2], vec![0, 1]);
        assert_eq!(counts.group_counts()[0], vec![2, 1]);
        assert_eq!(counts.group_counts()[2], vec![0, 1]);
    }

    #[test]
    fn counter_keeps_hard_and_fractional_semantics_separate() {
        let manifest = vec![SampleManifestRecord {
            sample: "sampleA".to_owned(),
            reads: PathBuf::from("a.bed"),
            group: Some("case".to_owned()),
            evidence: None,
        }];
        let mut counter =
            MultiSampleCounter::new(&manifest, vec!["TU_A".to_owned(), "TU_B".to_owned()]);

        counter
            .add_assignment_metrics("sampleA::unique", "TU_A", true, false, true, 1.0)
            .unwrap();
        counter
            .add_assignment_metrics("sampleA::ambiguous", "TU_A", false, false, false, 0.5)
            .unwrap();
        counter
            .add_assignment_metrics("sampleA::ambiguous", "TU_B", false, false, false, 0.5)
            .unwrap();
        let counts = counter.finish();

        assert_eq!(counts.counts(), &[vec![1], vec![0]]);
        assert_eq!(counts.unique_counts(), &[vec![1], vec![0]]);
        assert_eq!(counts.full_length_evidence_counts(), &[vec![0], vec![0]]);
        assert_eq!(counts.fractional_counts(), &[vec![1.5], vec![0.5]]);
        assert_eq!(
            counts
                .fractional_counts()
                .iter()
                .map(|row| row[0])
                .sum::<f64>(),
            2.0,
            "each of the two reads contributes total fractional weight 1.0"
        );
        assert_eq!(counts.group_counts(), &[vec![1], vec![0]]);
        assert_eq!(counts.group_fractional_counts(), &[vec![1.5], vec![0.5]]);
    }

    #[test]
    fn filtering_retains_every_aligned_metric_row() {
        let manifest = vec![
            SampleManifestRecord {
                sample: "sampleA".to_owned(),
                reads: PathBuf::from("a.bed"),
                group: Some("case".to_owned()),
                evidence: None,
            },
            SampleManifestRecord {
                sample: "sampleB".to_owned(),
                reads: PathBuf::from("b.bed"),
                group: Some("control".to_owned()),
                evidence: None,
            },
        ];
        let mut counter =
            MultiSampleCounter::new(&manifest, vec!["TU_keep".to_owned(), "TU_drop".to_owned()]);
        counter
            .add_assignment_metrics("sampleA::r1", "TU_keep", true, true, true, 1.0)
            .unwrap();
        counter
            .add_assignment_metrics("sampleB::r2", "TU_keep", false, false, true, 1.0)
            .unwrap();
        counter
            .add_assignment_metrics("sampleA::r3", "TU_drop", true, false, true, 1.0)
            .unwrap();

        let counts = counter.finish().filter_min_total_count(2);
        assert_eq!(counts.tu_ids(), &["TU_keep"]);
        assert_eq!(counts.samples().len(), 2);
        assert_eq!(counts.groups().len(), 2);
        for rows in [
            counts.counts(),
            counts.unique_counts(),
            counts.full_length_evidence_counts(),
            counts.group_counts(),
            counts.group_unique_counts(),
            counts.group_full_length_evidence_counts(),
        ] {
            assert_eq!(rows.len(), 1);
            assert_eq!(rows[0].len(), 2);
        }
        for rows in [counts.fractional_counts(), counts.group_fractional_counts()] {
            assert_eq!(rows.len(), 1);
            assert_eq!(rows[0].len(), 2);
        }
    }
}
