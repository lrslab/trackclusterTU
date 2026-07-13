use std::collections::HashMap;
use std::ffi::OsString;
use std::fs;
use std::path::{Component, Path, PathBuf};
use std::sync::atomic::{AtomicU64, Ordering};

use anyhow::{Context, Result};

static NEXT_TEMP_ID: AtomicU64 = AtomicU64::new(0);

#[derive(Debug)]
struct OutputEntry {
    final_path: PathBuf,
    staged_path: PathBuf,
    backup_path: PathBuf,
}

/// Stages a related set of output files and publishes them only after the
/// complete command succeeds.
///
/// Output paths are compared using their canonical existing ancestor so that
/// aliases through symlinked directories are rejected before any final output
/// is truncated. Existing outputs are backed up during commit and restored if
/// publishing any staged file fails.
#[derive(Debug)]
pub(crate) struct OutputTransaction {
    entries: Vec<OutputEntry>,
    final_to_entry: HashMap<PathBuf, usize>,
    committed: bool,
}

impl OutputTransaction {
    pub(crate) fn new<I, O>(inputs: I, outputs: O) -> Result<Self>
    where
        I: IntoIterator,
        I::Item: AsRef<Path>,
        O: IntoIterator,
        O::Item: AsRef<Path>,
    {
        let input_paths: Vec<PathBuf> = inputs
            .into_iter()
            .map(|path| path.as_ref().to_path_buf())
            .collect();
        let output_paths: Vec<PathBuf> = outputs
            .into_iter()
            .map(|path| path.as_ref().to_path_buf())
            .collect();

        let mut input_identities: HashMap<PathBuf, PathBuf> = HashMap::new();
        for path in input_paths {
            let identity = path_identity(&path)
                .with_context(|| format!("failed to resolve input path {}", path.display()))?;
            input_identities.entry(identity).or_insert(path);
        }

        let mut output_identities: HashMap<PathBuf, PathBuf> = HashMap::new();
        let mut entries = Vec::with_capacity(output_paths.len());
        let mut final_to_entry = HashMap::with_capacity(output_paths.len());

        for final_path in output_paths {
            if final_path.is_dir() {
                anyhow::bail!(
                    "output path is an existing directory: {}",
                    final_path.display()
                );
            }
            let identity = path_identity(&final_path).with_context(|| {
                format!("failed to resolve output path {}", final_path.display())
            })?;

            if let Some(previous) = output_identities.insert(identity.clone(), final_path.clone()) {
                anyhow::bail!(
                    "output paths alias the same file: {} and {}",
                    previous.display(),
                    final_path.display()
                );
            }
            if let Some(input) = input_identities.get(&identity) {
                anyhow::bail!(
                    "output path {} aliases input {}",
                    final_path.display(),
                    input.display()
                );
            }

            let parent = final_path.parent().unwrap_or_else(|| Path::new("."));
            fs::create_dir_all(parent).with_context(|| {
                format!("failed to create output directory {}", parent.display())
            })?;

            let staged_path = unique_sibling_path(&final_path, "staged")?;
            let backup_path = unique_sibling_path(&final_path, "backup")?;
            let entry_idx = entries.len();
            final_to_entry.insert(final_path.clone(), entry_idx);
            entries.push(OutputEntry {
                final_path,
                staged_path,
                backup_path,
            });
        }

        Ok(Self {
            entries,
            final_to_entry,
            committed: false,
        })
    }

    pub(crate) fn staged_path(&self, final_path: &Path) -> Result<PathBuf> {
        let idx = self.final_to_entry.get(final_path).ok_or_else(|| {
            anyhow::anyhow!(
                "output path {} was not registered in the transaction",
                final_path.display()
            )
        })?;
        Ok(self.entries[*idx].staged_path.clone())
    }

    pub(crate) fn commit(mut self) -> Result<()> {
        for entry in &self.entries {
            if entry.staged_path.exists() {
                fs::OpenOptions::new()
                    .read(true)
                    .open(&entry.staged_path)
                    .and_then(|file| file.sync_all())
                    .with_context(|| {
                        format!(
                            "failed to flush staged output {}",
                            entry.staged_path.display()
                        )
                    })?;
            }
        }

        let mut backed_up: Vec<usize> = Vec::new();
        for (idx, entry) in self.entries.iter().enumerate() {
            if entry.final_path.exists() {
                if let Err(error) = fs::rename(&entry.final_path, &entry.backup_path) {
                    restore_backups(&self.entries, &backed_up);
                    return Err(error).with_context(|| {
                        format!(
                            "failed to stage existing output {} for replacement",
                            entry.final_path.display()
                        )
                    });
                }
                backed_up.push(idx);
            }
        }

        let mut published: Vec<usize> = Vec::new();
        for (idx, entry) in self.entries.iter().enumerate() {
            if !entry.staged_path.exists() {
                continue;
            }
            if let Err(error) = fs::rename(&entry.staged_path, &entry.final_path) {
                for published_idx in published.iter().rev().copied() {
                    let published_entry = &self.entries[published_idx];
                    let _ = fs::remove_file(&published_entry.final_path);
                }
                restore_backups(&self.entries, &backed_up);
                return Err(error).with_context(|| {
                    format!(
                        "failed to publish staged output {}",
                        entry.final_path.display()
                    )
                });
            }
            published.push(idx);
        }

        self.committed = true;
        for idx in backed_up {
            let backup = &self.entries[idx].backup_path;
            if let Err(error) = fs::remove_file(backup) {
                eprintln!(
                    "warning: published outputs successfully but could not remove backup {}: {error}",
                    backup.display()
                );
            }
        }

        Ok(())
    }
}

impl Drop for OutputTransaction {
    fn drop(&mut self) {
        if self.committed {
            return;
        }
        for entry in &self.entries {
            let _ = fs::remove_file(&entry.staged_path);
            if entry.backup_path.exists() && !entry.final_path.exists() {
                let _ = fs::rename(&entry.backup_path, &entry.final_path);
            } else {
                let _ = fs::remove_file(&entry.backup_path);
            }
        }
    }
}

fn restore_backups(entries: &[OutputEntry], backed_up: &[usize]) {
    for idx in backed_up.iter().rev().copied() {
        let entry = &entries[idx];
        if entry.final_path.exists() {
            let _ = fs::remove_file(&entry.final_path);
        }
        let _ = fs::rename(&entry.backup_path, &entry.final_path);
    }
}

fn unique_sibling_path(final_path: &Path, kind: &str) -> Result<PathBuf> {
    let parent = final_path.parent().unwrap_or_else(|| Path::new("."));
    let file_name = final_path
        .file_name()
        .ok_or_else(|| anyhow::anyhow!("output path has no file name: {}", final_path.display()))?;

    for _ in 0..1_000 {
        let id = NEXT_TEMP_ID.fetch_add(1, Ordering::Relaxed);
        let mut candidate_name = OsString::from(".");
        candidate_name.push(file_name);
        candidate_name.push(format!(
            ".trackclustertu.{kind}.{}.{}",
            std::process::id(),
            id
        ));
        let candidate = parent.join(candidate_name);
        if !candidate.exists() {
            return Ok(candidate);
        }
    }

    anyhow::bail!(
        "failed to allocate a temporary sibling for {}",
        final_path.display()
    )
}

fn path_identity(path: &Path) -> Result<PathBuf> {
    let absolute = if path.is_absolute() {
        path.to_path_buf()
    } else {
        std::env::current_dir()?.join(path)
    };

    if absolute.exists() {
        return fs::canonicalize(&absolute)
            .with_context(|| format!("failed to canonicalize {}", absolute.display()));
    }

    let mut ancestor = absolute.as_path();
    let mut missing_parts: Vec<OsString> = Vec::new();
    while !ancestor.exists() {
        let file_name = ancestor.file_name().ok_or_else(|| {
            anyhow::anyhow!(
                "could not find an existing ancestor for {}",
                absolute.display()
            )
        })?;
        missing_parts.push(file_name.to_os_string());
        ancestor = ancestor.parent().ok_or_else(|| {
            anyhow::anyhow!(
                "could not find an existing ancestor for {}",
                absolute.display()
            )
        })?;
    }

    let mut resolved = fs::canonicalize(ancestor)
        .with_context(|| format!("failed to canonicalize {}", ancestor.display()))?;
    for part in missing_parts.into_iter().rev() {
        resolved.push(part);
    }
    Ok(lexically_normalize(&resolved))
}

fn lexically_normalize(path: &Path) -> PathBuf {
    let mut normalized = PathBuf::new();
    for component in path.components() {
        match component {
            Component::CurDir => {}
            Component::ParentDir => {
                normalized.pop();
            }
            other => normalized.push(other.as_os_str()),
        }
    }
    normalized
}

#[cfg(test)]
mod tests {
    use std::fs;
    use std::time::{SystemTime, UNIX_EPOCH};

    use super::*;

    fn temp_dir(name: &str) -> PathBuf {
        let nonce = SystemTime::now()
            .duration_since(UNIX_EPOCH)
            .expect("clock")
            .as_nanos();
        std::env::temp_dir().join(format!("trackclustertu_{name}_{nonce}"))
    }

    #[test]
    fn rejects_duplicate_outputs_before_truncation() {
        let dir = temp_dir("duplicate_outputs");
        fs::create_dir_all(&dir).unwrap();
        let output = dir.join("result.tsv");
        fs::write(&output, "existing\n").unwrap();

        let error = OutputTransaction::new(Vec::<PathBuf>::new(), [&output, &output]).unwrap_err();
        assert!(error.to_string().contains("alias the same file"));
        assert_eq!(fs::read_to_string(&output).unwrap(), "existing\n");

        let _ = fs::remove_dir_all(dir);
    }

    #[test]
    fn rejects_output_that_aliases_input() {
        let dir = temp_dir("input_alias");
        fs::create_dir_all(&dir).unwrap();
        let input = dir.join("reads.bed");
        fs::write(&input, "input\n").unwrap();

        let error = OutputTransaction::new([&input], [&input]).unwrap_err();
        assert!(error.to_string().contains("aliases input"));
        assert_eq!(fs::read_to_string(&input).unwrap(), "input\n");

        let _ = fs::remove_dir_all(dir);
    }

    #[test]
    fn drop_preserves_existing_output_and_removes_staging() {
        let dir = temp_dir("drop_rollback");
        fs::create_dir_all(&dir).unwrap();
        let output = dir.join("result.tsv");
        fs::write(&output, "old\n").unwrap();

        {
            let transaction = OutputTransaction::new(Vec::<PathBuf>::new(), [&output]).unwrap();
            let staged = transaction.staged_path(&output).unwrap();
            fs::write(&staged, "new\n").unwrap();
        }

        assert_eq!(fs::read_to_string(&output).unwrap(), "old\n");
        assert_eq!(fs::read_dir(&dir).unwrap().count(), 1);

        let _ = fs::remove_dir_all(dir);
    }

    #[test]
    fn commit_replaces_outputs_as_a_set_and_removes_absent_optional_output() {
        let dir = temp_dir("commit");
        fs::create_dir_all(&dir).unwrap();
        let first = dir.join("first.tsv");
        let second = dir.join("second.tsv");
        fs::write(&first, "old-first\n").unwrap();
        fs::write(&second, "stale-optional\n").unwrap();

        let transaction = OutputTransaction::new(Vec::<PathBuf>::new(), [&first, &second]).unwrap();
        fs::write(transaction.staged_path(&first).unwrap(), "new-first\n").unwrap();
        transaction.commit().unwrap();

        assert_eq!(fs::read_to_string(&first).unwrap(), "new-first\n");
        assert!(!second.exists());
        assert_eq!(fs::read_dir(&dir).unwrap().count(), 1);

        let _ = fs::remove_dir_all(dir);
    }
}
