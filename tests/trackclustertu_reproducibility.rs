#![cfg(unix)]

use std::fs;
use std::os::unix::fs::PermissionsExt;
use std::process::Command;
use std::time::{SystemTime, UNIX_EPOCH};

use noodles::{bam, sam};
use serde_json::Value;

fn unique_tmp_dir(prefix: &str) -> std::path::PathBuf {
    let nanos = SystemTime::now()
        .duration_since(UNIX_EPOCH)
        .expect("clock")
        .as_nanos();
    std::env::temp_dir().join(format!("{prefix}_{nanos}"))
}

fn make_executable(path: &std::path::Path) {
    let mut permissions = fs::metadata(path).unwrap().permissions();
    permissions.set_mode(0o755);
    fs::set_permissions(path, permissions).unwrap();
}

fn create_empty_bam(path: &std::path::Path) {
    let header = sam::Header::default();
    let mut writer = bam::io::Writer::new(fs::File::create(path).unwrap());
    writer.write_header(&header).unwrap();
    writer.try_finish().unwrap();
}

fn install_fake_tools(
    directory: &std::path::Path,
    fixture_bam: &std::path::Path,
    argument_log: &std::path::Path,
) -> (std::path::PathBuf, std::path::PathBuf) {
    fs::create_dir_all(directory).unwrap();
    let minimap2 = directory.join("fake minimap2 executable");
    fs::write(
        &minimap2,
        r#"#!/bin/sh
if [ "$1" = "--version" ]; then
  printf 'minimap2 fake 2.99\n'
  exit 0
fi
printf '%s\n' "$@" > "$MINIMAP2_ARGUMENT_LOG"
printf 'fake-sam\n'
"#,
    )
    .unwrap();
    make_executable(&minimap2);

    let samtools = directory.join("fake samtools executable");
    fs::write(
        &samtools,
        r#"#!/bin/sh
if [ "$1" = "--version" ]; then
  printf 'samtools fake 9.99\n'
  exit 0
fi
cmd="$1"
shift
case "$cmd" in
  view)
    cat
    ;;
  sort)
    out=""
    while [ $# -gt 0 ]; do
      case "$1" in
        -o)
          out="$2"
          shift 2
          ;;
        -@)
          shift 2
          ;;
        -)
          shift
          ;;
        *)
          shift
          ;;
      esac
    done
    cat >/dev/null
    cp "$FAKE_BAM_PATH" "$out"
    ;;
  index)
    : > "$1.bai"
    ;;
  *)
    printf 'unexpected samtools command: %s\n' "$cmd" >&2
    exit 1
    ;;
esac
"#,
    )
    .unwrap();
    make_executable(&samtools);

    assert!(!argument_log.exists());
    assert!(fixture_bam.exists());
    (minimap2, samtools)
}

struct CommandFixture<'a> {
    manifest: &'a std::path::Path,
    reference: &'a std::path::Path,
    minimap2: &'a std::path::Path,
    samtools: &'a std::path::Path,
    fixture_bam: &'a std::path::Path,
    argument_log: &'a std::path::Path,
}

fn base_command(command: &str, out_dir: &std::path::Path, fixture: &CommandFixture<'_>) -> Command {
    let mut process = Command::new(env!("CARGO_BIN_EXE_trackclustertu"));
    process
        .env("FAKE_BAM_PATH", fixture.fixture_bam)
        .env("MINIMAP2_ARGUMENT_LOG", fixture.argument_log)
        .arg(command)
        .arg("--manifest")
        .arg(fixture.manifest)
        .arg("--reference-fasta")
        .arg(fixture.reference)
        .arg("--out-dir")
        .arg(out_dir)
        .arg("--threads")
        .arg("3")
        .arg("--minimap2")
        .arg(fixture.minimap2)
        .arg("--samtools")
        .arg(fixture.samtools)
        .arg("--minimap2-arg")
        .arg("-ax")
        .arg("--minimap2-arg")
        .arg("one minimap argument with spaces");
    process
}

fn assert_success(output: &std::process::Output) {
    assert!(
        output.status.success(),
        "stdout:\n{}\nstderr:\n{}",
        String::from_utf8_lossy(&output.stdout),
        String::from_utf8_lossy(&output.stderr)
    );
}

#[test]
fn map_and_run_write_reproducible_atomic_manifests() {
    let tmp = unique_tmp_dir("trackclustertu reproducibility with spaces");
    let tools_dir = tmp.join("tool directory with spaces");
    let fixture_bam = tmp.join("fixture bam.bam");
    let argument_log = tmp.join("minimap argument log.txt");
    fs::create_dir_all(&tmp).unwrap();
    create_empty_bam(&fixture_bam);
    let (minimap2, samtools) = install_fake_tools(&tools_dir, &fixture_bam, &argument_log);

    let reference = tmp.join("reference genome.fa");
    fs::write(&reference, ">chr1\nACGTACGT\n").unwrap();
    let fastq = tmp.join("reads with spaces.fastq");
    fs::write(&fastq, "@read1\nACGT\n+\nIIII\n").unwrap();
    let input_manifest = tmp.join("input manifest.tsv");
    fs::write(
        &input_manifest,
        format!(
            "sample\tgroup\treads\nsample one\tcontrol\t{}\n",
            fastq.display()
        ),
    )
    .unwrap();

    let map_out = tmp.join("map output with spaces");
    let fixture = CommandFixture {
        manifest: &input_manifest,
        reference: &reference,
        minimap2: &minimap2,
        samtools: &samtools,
        fixture_bam: &fixture_bam,
        argument_log: &argument_log,
    };
    let map_output = base_command("map", &map_out, &fixture).output().unwrap();
    assert_success(&map_output);

    assert_eq!(
        fs::read_to_string(&argument_log).unwrap(),
        format!(
            "-ax\nmap-ont\n-ax\none minimap argument with spaces\n-t\n2\n{}\n{}\n",
            reference.canonicalize().unwrap().display(),
            fastq.canonicalize().unwrap().display(),
        )
    );
    let map_manifest_path = map_out.join("run_manifest.json");
    let map_manifest_bytes = fs::read(&map_manifest_path).unwrap();
    let map_manifest: Value = serde_json::from_slice(&map_manifest_bytes).unwrap();
    assert_eq!(map_manifest["schema"], "trackclustertu.run-manifest.v1");
    assert_eq!(map_manifest["command"], "map");
    assert_eq!(map_manifest["library_profile"], "direct-rna");
    if let Some(revision) = option_env!("TRACKCLUSTERTU_GIT_REVISION") {
        assert_eq!(map_manifest["source"]["git_revision"], revision);
    }
    if let Some(dirty) = option_env!("TRACKCLUSTERTU_GIT_DIRTY") {
        assert_eq!(map_manifest["source"]["git_dirty"], dirty == "true");
    }
    assert_eq!(
        map_manifest["external_tools"][0]["version"],
        "minimap2 fake 2.99"
    );
    assert_eq!(
        map_manifest["external_tools"][1]["version"],
        "samtools fake 9.99"
    );
    assert_eq!(
        map_manifest["effective_configuration"]["minimap2_args"][3]["value"],
        "one minimap argument with spaces"
    );
    assert!(map_manifest["inputs"]
        .as_array()
        .unwrap()
        .iter()
        .all(|input| {
            input["canonical_path"]["display"]
                .as_str()
                .is_some_and(|path| std::path::Path::new(path).is_absolute())
                && input["sha256"]
                    .as_str()
                    .is_some_and(|hash| hash.len() == 64)
        }));

    // A failed rerun must not replace the last successful manifest.
    fs::write(
        &minimap2,
        "#!/bin/sh\nif [ \"$1\" = \"--version\" ]; then echo 'broken 1'; exit 0; fi\nexit 17\n",
    )
    .unwrap();
    make_executable(&minimap2);
    let failed = base_command("map", &map_out, &fixture).output().unwrap();
    assert!(!failed.status.success());
    assert_eq!(fs::read(&map_manifest_path).unwrap(), map_manifest_bytes);

    // Restore the mapper and verify the full command publishes a run-level manifest only after
    // clustering succeeds.
    fs::write(
        &minimap2,
        r#"#!/bin/sh
if [ "$1" = "--version" ]; then printf 'minimap2 fake 2.99\n'; exit 0; fi
printf '%s\n' "$@" > "$MINIMAP2_ARGUMENT_LOG"
printf 'fake-sam\n'
"#,
    )
    .unwrap();
    make_executable(&minimap2);
    let run_out = tmp.join("run output with spaces");
    let run_output = base_command("run", &run_out, &fixture).output().unwrap();
    assert_success(&run_output);
    let run_manifest: Value =
        serde_json::from_slice(&fs::read(run_out.join("run_manifest.json")).unwrap()).unwrap();
    assert_eq!(run_manifest["command"], "run");
    assert_eq!(
        run_manifest["effective_configuration"]["clustering"]["threads"],
        3
    );
    assert!(run_manifest["effective_configuration"]["mapping"].is_object());

    fs::remove_dir_all(tmp).unwrap();
}
