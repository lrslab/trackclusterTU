use std::collections::BTreeMap;
use std::ffi::OsStr;
use std::fs::{self, File};
use std::io::{BufReader, BufWriter, Read, Write};
use std::path::{Path, PathBuf};
use std::process::Command;
use std::time::{SystemTime, UNIX_EPOCH};

use anyhow::{Context, Result};
use serde::Serialize;
use serde_json::Value;
use sha2::{Digest, Sha256};

use crate::tools::output_transaction::OutputTransaction;

pub(crate) const RUN_MANIFEST_FILE_NAME: &str = "run_manifest.json";
const RUN_MANIFEST_SCHEMA: &str = "trackclustertu.run-manifest.v1";
const VERSION_OUTPUT_LIMIT: usize = 16 * 1024;

/// A JSON-safe, lossless representation of an operating-system string.
///
/// UTF-8 values remain directly readable. Non-UTF-8 values retain their raw Unix bytes or
/// Windows wide units as hexadecimal, while `display` remains a human-readable lossy rendering.
#[derive(Clone, Debug, Eq, PartialEq, Serialize)]
pub(crate) struct ManifestOsValue {
    pub display: String,
    pub encoding: &'static str,
    pub value: String,
}

impl ManifestOsValue {
    pub(crate) fn new(value: &OsStr) -> Self {
        if let Some(utf8) = value.to_str() {
            return Self {
                display: utf8.to_owned(),
                encoding: "utf-8",
                value: utf8.to_owned(),
            };
        }

        #[cfg(unix)]
        {
            use std::os::unix::ffi::OsStrExt;

            return Self {
                display: value.to_string_lossy().into_owned(),
                encoding: "unix-bytes-hex",
                value: hex_bytes(value.as_bytes()),
            };
        }

        #[cfg(windows)]
        {
            use std::os::windows::ffi::OsStrExt;

            let mut encoded = String::new();
            for unit in value.encode_wide() {
                use std::fmt::Write as _;
                let _ = write!(encoded, "{unit:04x}");
            }
            return Self {
                display: value.to_string_lossy().into_owned(),
                encoding: "windows-wide-hex",
                value: encoded,
            };
        }

        #[allow(unreachable_code)]
        Self {
            display: value.to_string_lossy().into_owned(),
            encoding: "platform-lossy",
            value: value.to_string_lossy().into_owned(),
        }
    }

    pub(crate) fn path(path: &Path) -> Self {
        Self::new(path.as_os_str())
    }
}

#[cfg(unix)]
fn hex_bytes(bytes: &[u8]) -> String {
    use std::fmt::Write as _;

    let mut encoded = String::with_capacity(bytes.len() * 2);
    for byte in bytes {
        let _ = write!(encoded, "{byte:02x}");
    }
    encoded
}

#[derive(Clone, Debug, Eq, PartialEq)]
pub(crate) struct Executable {
    requested: PathBuf,
    resolved: PathBuf,
}

impl Executable {
    pub(crate) fn resolve(requested: PathBuf, label: &str) -> Result<Self> {
        if requested.as_os_str().is_empty() {
            anyhow::bail!("--{label} executable path cannot be empty");
        }

        let resolved = resolve_executable_path(&requested).with_context(|| {
            format!(
                "failed to resolve --{label} executable {}",
                requested.display()
            )
        })?;
        let metadata = fs::metadata(&resolved).with_context(|| {
            format!(
                "failed to inspect --{label} executable {}",
                resolved.display()
            )
        })?;
        if !metadata.is_file() {
            anyhow::bail!(
                "--{label} executable is not a regular file: {}",
                resolved.display()
            );
        }

        Ok(Self {
            requested,
            resolved,
        })
    }

    pub(crate) fn requested(&self) -> &Path {
        &self.requested
    }

    pub(crate) fn resolved(&self) -> &Path {
        &self.resolved
    }
}

fn resolve_executable_path(requested: &Path) -> Result<PathBuf> {
    if requested.is_absolute() || requested.components().count() > 1 {
        return requested
            .canonicalize()
            .with_context(|| format!("executable not found: {}", requested.display()));
    }

    let path_value = std::env::var_os("PATH").context("PATH is not set")?;
    for directory in std::env::split_paths(&path_value) {
        let candidate = directory.join(requested);
        if candidate.is_file() {
            return candidate.canonicalize().with_context(|| {
                format!("failed to canonicalize executable {}", candidate.display())
            });
        }

        #[cfg(windows)]
        {
            if candidate.extension().is_none() {
                if let Some(path_ext) = std::env::var_os("PATHEXT") {
                    for extension in path_ext.to_string_lossy().split(';') {
                        let extension = extension.trim_start_matches('.');
                        let candidate = candidate.with_extension(extension);
                        if candidate.is_file() {
                            return candidate.canonicalize().with_context(|| {
                                format!("failed to canonicalize executable {}", candidate.display())
                            });
                        }
                    }
                }
            }
        }
    }

    anyhow::bail!(
        "executable {:?} was not found on PATH; pass an explicit path",
        requested
    )
}

#[derive(Clone, Debug, Serialize)]
pub(crate) struct InputFingerprint {
    role: String,
    #[serde(skip_serializing_if = "Option::is_none")]
    sample: Option<String>,
    canonical_path: ManifestOsValue,
    size_bytes: u64,
    sha256: String,
}

#[derive(Clone, Debug)]
pub(crate) struct InputDescriptor {
    role: String,
    sample: Option<String>,
    path: PathBuf,
}

impl InputDescriptor {
    pub(crate) fn new(role: impl Into<String>, path: impl Into<PathBuf>) -> Self {
        Self {
            role: role.into(),
            sample: None,
            path: path.into(),
        }
    }

    pub(crate) fn for_sample(
        role: impl Into<String>,
        sample: impl Into<String>,
        path: impl Into<PathBuf>,
    ) -> Self {
        Self {
            role: role.into(),
            sample: Some(sample.into()),
            path: path.into(),
        }
    }

    pub(crate) fn path(&self) -> &Path {
        &self.path
    }
}

#[derive(Clone, Debug, Serialize)]
pub(crate) struct ExternalToolRecord {
    name: String,
    requested_path: ManifestOsValue,
    resolved_path: ManifestOsValue,
    #[serde(skip_serializing_if = "Option::is_none")]
    version: Option<String>,
    #[serde(skip_serializing_if = "Option::is_none")]
    version_error: Option<String>,
}

#[derive(Clone, Debug, Serialize)]
struct PackageRecord {
    name: &'static str,
    version: &'static str,
}

#[derive(Clone, Debug, Default, Serialize)]
struct SourceRecord {
    #[serde(skip_serializing_if = "Option::is_none")]
    git_revision: Option<String>,
    #[serde(skip_serializing_if = "Option::is_none")]
    git_dirty: Option<bool>,
}

#[derive(Clone, Debug, Serialize)]
pub(crate) struct RunManifest {
    schema: &'static str,
    command: String,
    timestamp_unix_seconds: u64,
    package: PackageRecord,
    source: SourceRecord,
    effective_configuration: Value,
    library_profile: String,
    external_tools: Vec<ExternalToolRecord>,
    inputs: Vec<InputFingerprint>,
}

impl RunManifest {
    pub(crate) fn capture(
        command: impl Into<String>,
        effective_configuration: Value,
        library_profile: impl Into<String>,
        input_descriptors: &[InputDescriptor],
        tools: &[(&str, &Executable)],
        timestamp_unix_seconds: u64,
    ) -> Result<Self> {
        let inputs = fingerprint_inputs(input_descriptors)?;
        let external_tools = tools
            .iter()
            .map(|(name, executable)| capture_tool_version(name, executable))
            .collect();

        Ok(Self {
            schema: RUN_MANIFEST_SCHEMA,
            command: command.into(),
            timestamp_unix_seconds,
            package: PackageRecord {
                name: env!("CARGO_PKG_NAME"),
                version: env!("CARGO_PKG_VERSION"),
            },
            source: capture_source_revision(),
            effective_configuration,
            library_profile: library_profile.into(),
            external_tools,
            inputs,
        })
    }

    #[cfg(test)]
    fn to_pretty_json(&self) -> Result<Vec<u8>> {
        let mut bytes = serde_json::to_vec_pretty(self)?;
        bytes.push(b'\n');
        Ok(bytes)
    }
}

pub(crate) fn timestamp_now() -> Result<u64> {
    Ok(SystemTime::now()
        .duration_since(UNIX_EPOCH)
        .context("system clock is before the Unix epoch")?
        .as_secs())
}

pub(crate) struct RunManifestPublisher {
    final_path: PathBuf,
    transaction: OutputTransaction,
}

impl RunManifestPublisher {
    pub(crate) fn new(out_dir: &Path, inputs: &[InputDescriptor]) -> Result<Self> {
        let final_path = out_dir.join(RUN_MANIFEST_FILE_NAME);
        let transaction = OutputTransaction::new(
            inputs.iter().map(InputDescriptor::path),
            std::iter::once(&final_path),
        )?;
        Ok(Self {
            final_path,
            transaction,
        })
    }

    pub(crate) fn final_path(&self) -> &Path {
        &self.final_path
    }

    pub(crate) fn publish(self, manifest: &RunManifest) -> Result<()> {
        let staged_path = self.transaction.staged_path(&self.final_path)?;
        let file = File::create(&staged_path)
            .with_context(|| format!("failed to create staged run manifest {staged_path:?}"))?;
        let mut writer = BufWriter::new(file);
        serde_json::to_writer_pretty(&mut writer, manifest)
            .with_context(|| format!("failed to serialize run manifest {staged_path:?}"))?;
        writer.write_all(b"\n")?;
        writer.flush()?;
        drop(writer);
        self.transaction.commit()
    }
}

fn fingerprint_inputs(descriptors: &[InputDescriptor]) -> Result<Vec<InputFingerprint>> {
    let mut digest_cache: BTreeMap<PathBuf, (u64, String)> = BTreeMap::new();
    let mut fingerprints = Vec::with_capacity(descriptors.len());
    for descriptor in descriptors {
        let canonical_path = descriptor.path.canonicalize().with_context(|| {
            format!(
                "failed to canonicalize {} input {}",
                descriptor.role,
                descriptor.path.display()
            )
        })?;
        let (size_bytes, sha256) = if let Some(cached) = digest_cache.get(&canonical_path) {
            cached.clone()
        } else {
            let fingerprint = sha256_file(&canonical_path).with_context(|| {
                format!(
                    "failed to fingerprint {} input {}",
                    descriptor.role,
                    canonical_path.display()
                )
            })?;
            digest_cache.insert(canonical_path.clone(), fingerprint.clone());
            fingerprint
        };
        fingerprints.push(InputFingerprint {
            role: descriptor.role.clone(),
            sample: descriptor.sample.clone(),
            canonical_path: ManifestOsValue::path(&canonical_path),
            size_bytes,
            sha256,
        });
    }
    Ok(fingerprints)
}

fn sha256_file(path: &Path) -> Result<(u64, String)> {
    let file = File::open(path)?;
    let mut reader = BufReader::new(file);
    let mut hasher = Sha256::new();
    let mut buffer = [0_u8; 64 * 1024];
    let mut size_bytes = 0_u64;
    loop {
        let read = reader.read(&mut buffer)?;
        if read == 0 {
            break;
        }
        size_bytes += read as u64;
        hasher.update(&buffer[..read]);
    }
    Ok((size_bytes, format!("{:x}", hasher.finalize())))
}

fn capture_tool_version(name: &str, executable: &Executable) -> ExternalToolRecord {
    let output = Command::new(executable.resolved())
        .arg("--version")
        .output();
    let (version, version_error) = match output {
        Ok(output) => {
            let combined = if output.stdout.is_empty() {
                &output.stderr
            } else {
                &output.stdout
            };
            let version = bounded_text(combined);
            if output.status.success() {
                (nonempty(version), None)
            } else {
                let detail = nonempty(version)
                    .map(|text| format!("{}: {text}", output.status))
                    .unwrap_or_else(|| output.status.to_string());
                (None, Some(detail))
            }
        }
        Err(error) => (None, Some(error.to_string())),
    };

    ExternalToolRecord {
        name: name.to_owned(),
        requested_path: ManifestOsValue::path(executable.requested()),
        resolved_path: ManifestOsValue::path(executable.resolved()),
        version,
        version_error,
    }
}

fn bounded_text(bytes: &[u8]) -> String {
    let end = bytes.len().min(VERSION_OUTPUT_LIMIT);
    let mut text = String::from_utf8_lossy(&bytes[..end]).trim().to_owned();
    if bytes.len() > VERSION_OUTPUT_LIMIT {
        text.push_str(" [truncated]");
    }
    text
}

fn nonempty(value: String) -> Option<String> {
    (!value.is_empty()).then_some(value)
}

fn capture_source_revision() -> SourceRecord {
    let revision = option_env!("TRACKCLUSTERTU_GIT_REVISION")
        .filter(|value| !value.is_empty())
        .map(str::to_owned);
    let git_dirty = revision.as_ref().and_then(|_| {
        option_env!("TRACKCLUSTERTU_GIT_DIRTY").and_then(|value| match value {
            "true" => Some(true),
            "false" => Some(false),
            _ => None,
        })
    });
    SourceRecord {
        git_revision: revision,
        git_dirty,
    }
}

#[cfg(test)]
mod tests {
    use std::time::{SystemTime, UNIX_EPOCH};

    #[cfg(unix)]
    use std::os::unix::fs::PermissionsExt;

    use serde_json::json;

    use super::*;

    fn temp_dir(label: &str) -> PathBuf {
        let id = SystemTime::now()
            .duration_since(UNIX_EPOCH)
            .unwrap()
            .as_nanos();
        std::env::temp_dir().join(format!("trackclustertu_{label}_{id}"))
    }

    #[test]
    fn source_revision_uses_compile_time_metadata() {
        let source = capture_source_revision();
        assert_eq!(
            source.git_revision.as_deref(),
            option_env!("TRACKCLUSTERTU_GIT_REVISION")
        );
        assert_eq!(
            source.git_dirty,
            option_env!("TRACKCLUSTERTU_GIT_DIRTY").and_then(|value| match value {
                "true" => Some(true),
                "false" => Some(false),
                _ => None,
            })
        );
    }

    #[test]
    fn fingerprints_sha256_and_json_escapes_strings() {
        let directory = temp_dir("manifest_json");
        fs::create_dir_all(&directory).unwrap();
        let input = directory.join("input with spaces.txt");
        fs::write(&input, b"abc").unwrap();
        let descriptors = vec![InputDescriptor::for_sample(
            "fastq",
            "sample\t\"quoted\"\nline",
            &input,
        )];
        let manifest = RunManifest::capture(
            "map",
            json!({"argument": "tab\tquote\"newline\n"}),
            "direct-rna",
            &descriptors,
            &[],
            42,
        )
        .unwrap();

        let bytes = manifest.to_pretty_json().unwrap();
        let parsed: Value = serde_json::from_slice(&bytes).unwrap();
        assert_eq!(parsed["timestamp_unix_seconds"], 42);
        assert_eq!(
            parsed["inputs"][0]["sha256"],
            "ba7816bf8f01cfea414140de5dae2223b00361a396177a9cb410ff61f20015ad"
        );
        assert_eq!(parsed["inputs"][0]["sample"], "sample\t\"quoted\"\nline");
        assert_eq!(
            parsed["effective_configuration"]["argument"],
            "tab\tquote\"newline\n"
        );

        fs::remove_dir_all(directory).unwrap();
    }

    #[test]
    fn serialization_is_deterministic_except_for_the_supplied_timestamp() {
        let directory = temp_dir("manifest_determinism");
        fs::create_dir_all(&directory).unwrap();
        let input = directory.join("input");
        fs::write(&input, b"stable").unwrap();
        let descriptors = vec![InputDescriptor::new("manifest", &input)];
        let first = RunManifest::capture(
            "map",
            json!({"threads": 4}),
            "direct-rna",
            &descriptors,
            &[],
            1,
        )
        .unwrap();
        let second = RunManifest::capture(
            "map",
            json!({"threads": 4}),
            "direct-rna",
            &descriptors,
            &[],
            2,
        )
        .unwrap();
        let mut first: Value = serde_json::from_slice(&first.to_pretty_json().unwrap()).unwrap();
        let mut second: Value = serde_json::from_slice(&second.to_pretty_json().unwrap()).unwrap();
        first
            .as_object_mut()
            .unwrap()
            .remove("timestamp_unix_seconds");
        second
            .as_object_mut()
            .unwrap()
            .remove("timestamp_unix_seconds");
        assert_eq!(first, second);

        fs::remove_dir_all(directory).unwrap();
    }

    #[cfg(unix)]
    #[test]
    fn captures_version_from_an_executable_path_containing_spaces() {
        let directory = temp_dir("tool version spaces");
        fs::create_dir_all(&directory).unwrap();
        let tool = directory.join("fake mapper tool");
        fs::write(&tool, "#!/bin/sh\nprintf 'mapper 7.8.9\\n'\n").unwrap();
        let mut permissions = fs::metadata(&tool).unwrap().permissions();
        permissions.set_mode(0o755);
        fs::set_permissions(&tool, permissions).unwrap();

        let executable = Executable::resolve(tool.clone(), "minimap2").unwrap();
        let record = capture_tool_version("minimap2", &executable);
        let value = serde_json::to_value(record).unwrap();
        assert_eq!(value["version"], "mapper 7.8.9");
        assert_eq!(
            value["resolved_path"]["display"],
            tool.canonicalize().unwrap().display().to_string()
        );

        fs::remove_dir_all(directory).unwrap();
    }

    #[test]
    fn publisher_replaces_the_manifest_as_one_atomic_file() {
        let directory = temp_dir("manifest_atomic");
        fs::create_dir_all(&directory).unwrap();
        let input = directory.join("input");
        fs::write(&input, b"input").unwrap();
        let descriptors = vec![InputDescriptor::new("manifest", &input)];
        let publisher = RunManifestPublisher::new(&directory, &descriptors).unwrap();
        let manifest =
            RunManifest::capture("map", json!({}), "direct-rna", &descriptors, &[], 5).unwrap();
        publisher.publish(&manifest).unwrap();

        let parsed: Value =
            serde_json::from_slice(&fs::read(directory.join(RUN_MANIFEST_FILE_NAME)).unwrap())
                .unwrap();
        assert_eq!(parsed["timestamp_unix_seconds"], 5);
        assert!(!fs::read_dir(&directory)
            .unwrap()
            .filter_map(Result::ok)
            .any(|entry| entry.file_name().to_string_lossy().contains(".staged.")));

        fs::remove_dir_all(directory).unwrap();
    }
}
