use std::collections::{BTreeSet, HashMap, HashSet};
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::{Path, PathBuf};

use thiserror::Error;

const POOLED_ID_SEPARATOR: &str = "::";

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct SampleManifestRecord {
    pub sample: String,
    pub reads: PathBuf,
    pub group: Option<String>,
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct MultiSampleCounts {
    pub tu_ids: Vec<String>,
    pub samples: Vec<String>,
    pub counts: Vec<Vec<u64>>,
    pub groups: Vec<String>,
    pub group_counts: Vec<Vec<u64>>,
}

impl MultiSampleCounts {
    pub fn has_groups(&self) -> bool {
        !self.groups.is_empty()
    }
}

#[derive(Error, Debug)]
pub enum MultiSampleError {
    #[error("I/O error reading {path:?}: {source}")]
    IoRead {
        path: PathBuf,
        source: std::io::Error,
    },

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
}

fn split_fields(line: &str) -> Vec<&str> {
    let mut fields: Vec<&str> = line.split('\t').collect();
    if fields.len() == 1 {
        fields = line.split_whitespace().collect();
    }
    fields
}

pub fn parse_manifest(path: &Path) -> Result<Vec<SampleManifestRecord>, MultiSampleError> {
    let path = path.to_path_buf();
    let file = File::open(&path).map_err(|source| MultiSampleError::IoRead {
        path: path.clone(),
        source,
    })?;
    let reader = BufReader::new(file);
    let base_dir = path.parent().unwrap_or(Path::new("."));

    let mut header_line: Option<usize> = None;
    let mut sample_col: Option<usize> = None;
    let mut reads_col: Option<usize> = None;
    let mut group_col: Option<usize> = None;

    let mut seen_samples: HashSet<String> = HashSet::new();
    let mut records: Vec<SampleManifestRecord> = Vec::new();

    for (line_idx, line_result) in reader.lines().enumerate() {
        let line_number = line_idx + 1;
        let line = line_result.map_err(|source| MultiSampleError::IoRead {
            path: path.clone(),
            source,
        })?;
        let line = line.trim();
        if line.is_empty() || line.starts_with('#') {
            continue;
        }

        let fields = split_fields(line);
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

        let sample = fields[sample_idx].trim().to_owned();
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

        let reads_value = fields[reads_idx].trim();
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

        records.push(SampleManifestRecord {
            sample,
            reads,
            group,
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

pub fn sorted_tu_ids<I>(tu_ids: I) -> Vec<String>
where
    I: IntoIterator,
    I::Item: AsRef<str>,
{
    let mut ordered: BTreeSet<String> = BTreeSet::new();
    for tu_id in tu_ids {
        ordered.insert(tu_id.as_ref().to_owned());
    }
    ordered.into_iter().collect()
}

#[derive(Clone, Debug)]
pub struct MultiSampleCounter {
    tu_ids: Vec<String>,
    samples: Vec<String>,
    counts: Vec<Vec<u64>>,
    groups: Vec<String>,
    group_counts: Vec<Vec<u64>>,
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
        let group_counts = vec![vec![0; groups.len()]; tu_ids.len()];

        Self {
            tu_ids,
            samples,
            counts,
            groups,
            group_counts,
            tu_to_idx,
            sample_to_idx,
            sample_to_group_idx,
        }
    }

    pub fn add_assignment(
        &mut self,
        pooled_read_id: &str,
        tu_id: &str,
    ) -> Result<(), MultiSampleError> {
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

        self.counts[tu_idx][sample_idx] += 1;
        if let Some(group_idx) = self.sample_to_group_idx[sample_idx] {
            self.group_counts[tu_idx][group_idx] += 1;
        }

        Ok(())
    }

    pub fn finish(self) -> MultiSampleCounts {
        MultiSampleCounts {
            tu_ids: self.tu_ids,
            samples: self.samples,
            counts: self.counts,
            groups: self.groups,
            group_counts: self.group_counts,
        }
    }
}

#[cfg(test)]
mod tests {
    use std::fs;
    use std::time::{SystemTime, UNIX_EPOCH};

    use super::*;

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
            "sample\treads\tgroup\nsampleA\treads/a.bed\tcase\nsampleB\treads/b.bed\tcontrol\n",
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
                },
                SampleManifestRecord {
                    sample: "sampleB".to_owned(),
                    reads: tmp.join("reads/b.bed"),
                    group: Some("control".to_owned()),
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
    fn counter_aggregates_sample_and_group_counts() {
        let manifest = vec![
            SampleManifestRecord {
                sample: "sampleA".to_owned(),
                reads: PathBuf::from("a.bed"),
                group: Some("control".to_owned()),
            },
            SampleManifestRecord {
                sample: "sampleB".to_owned(),
                reads: PathBuf::from("b.bed"),
                group: Some("treated".to_owned()),
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
        assert_eq!(counts.samples, vec!["sampleA", "sampleB"]);
        assert_eq!(counts.groups, vec!["control", "treated"]);
        assert_eq!(counts.counts[0], vec![2, 1]);
        assert_eq!(counts.counts[1], vec![0, 0]);
        assert_eq!(counts.counts[2], vec![0, 1]);
        assert_eq!(counts.group_counts[0], vec![2, 1]);
        assert_eq!(counts.group_counts[2], vec![0, 1]);
    }
}
