use std::env;
use std::path::Path;
use std::process::Command;

const REVISION_ENV: &str = "TRACKCLUSTERTU_GIT_REVISION";
const DIRTY_ENV: &str = "TRACKCLUSTERTU_GIT_DIRTY";

fn explicit_revision() -> Option<String> {
    let value = env::var(REVISION_ENV).ok()?;
    let value = value.trim();
    if value.is_empty() {
        return None;
    }
    assert!(
        (7..=64).contains(&value.len()) && value.bytes().all(|byte| byte.is_ascii_hexdigit()),
        "{REVISION_ENV} must be a 7-64 character hexadecimal Git revision"
    );
    Some(value.to_ascii_lowercase())
}

fn explicit_dirty() -> Option<bool> {
    let value = env::var(DIRTY_ENV).ok()?;
    match value.trim().to_ascii_lowercase().as_str() {
        "true" | "1" => Some(true),
        "false" | "0" => Some(false),
        _ => panic!("{DIRTY_ENV} must be true, false, 1, or 0"),
    }
}

fn git_output(repository: &Path, args: &[&str]) -> Option<String> {
    let output = Command::new("git")
        .arg("-C")
        .arg(repository)
        .args(args)
        .env("GIT_OPTIONAL_LOCKS", "0")
        .output()
        .ok()?;
    if !output.status.success() {
        return None;
    }
    Some(String::from_utf8_lossy(&output.stdout).trim().to_owned())
}

fn main() {
    println!("cargo:rerun-if-env-changed={REVISION_ENV}");
    println!("cargo:rerun-if-env-changed={DIRTY_ENV}");

    let repository = env::var_os("CARGO_MANIFEST_DIR").map(std::path::PathBuf::from);
    let revision = explicit_revision().or_else(|| {
        repository
            .as_deref()
            .and_then(|path| git_output(path, &["rev-parse", "--verify", "HEAD"]))
            .filter(|value| {
                (7..=64).contains(&value.len())
                    && value.bytes().all(|byte| byte.is_ascii_hexdigit())
            })
    });
    let dirty = revision.as_ref().and_then(|_| {
        explicit_dirty().or_else(|| {
            repository.as_deref().and_then(|path| {
                git_output(
                    path,
                    &["status", "--porcelain=v1", "--untracked-files=normal"],
                )
                .map(|status| !status.is_empty())
            })
        })
    });

    if let Some(revision) = revision {
        println!("cargo:rustc-env={REVISION_ENV}={revision}");
    }
    if let Some(dirty) = dirty {
        println!("cargo:rustc-env={DIRTY_ENV}={dirty}");
    }
}
