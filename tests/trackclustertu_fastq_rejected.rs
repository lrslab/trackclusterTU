use std::fs;
use std::io::Write;
use std::process::Command;
use std::time::{SystemTime, UNIX_EPOCH};

fn unique_tmp_dir(prefix: &str) -> std::path::PathBuf {
    let nanos = SystemTime::now()
        .duration_since(UNIX_EPOCH)
        .expect("clock")
        .as_nanos();
    std::env::temp_dir().join(format!("{prefix}_{nanos}"))
}

#[test]
fn trackclustertu_rejects_fastq_manifest_inputs() {
    let tmp = unique_tmp_dir("trackclustertu_fastq_rejected_test");
    fs::create_dir_all(&tmp).unwrap();

    let fastq_path = tmp.join("sampleA.fq");
    let manifest_path = tmp.join("samples.tsv");
    let out_dir = tmp.join("out");

    let mut fastq = fs::File::create(&fastq_path).unwrap();
    writeln!(fastq, "@r1").unwrap();
    writeln!(fastq, "ACGTACGT").unwrap();
    writeln!(fastq, "+").unwrap();
    writeln!(fastq, "IIIIIIII").unwrap();

    fs::write(
        &manifest_path,
        format!(
            "sample\tgroup\treads\nsampleA\tcontrol\t{}\n",
            fastq_path.display()
        ),
    )
    .unwrap();

    let exe = env!("CARGO_BIN_EXE_trackclustertu");
    let output = Command::new(exe)
        .args([
            "cluster",
            "--manifest",
            manifest_path.to_str().unwrap(),
            "--out-dir",
            out_dir.to_str().unwrap(),
        ])
        .output()
        .unwrap();

    assert!(
        !output.status.success(),
        "stdout:\n{}\nstderr:\n{}",
        String::from_utf8_lossy(&output.stdout),
        String::from_utf8_lossy(&output.stderr)
    );

    let stderr = String::from_utf8_lossy(&output.stderr);
    assert!(
        stderr.contains("cannot be FASTQ here")
            || stderr.contains("FASTQ inputs are not supported directly"),
        "stderr did not explain FASTQ preprocessing requirement:\n{stderr}",
    );

    let _ = fs::remove_dir_all(&tmp);
}
