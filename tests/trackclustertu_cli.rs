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
fn trackclustertu_clusters_and_writes_outputs() {
    let tmp = unique_tmp_dir("trackclustertu_test");
    fs::create_dir_all(&tmp).unwrap();

    let input_path = tmp.join("reads.bed");
    let out_tu = tmp.join("tu.bed");
    let out_membership = tmp.join("membership.tsv");
    let out_tu_count = tmp.join("tu_count.csv");

    let mut input = fs::File::create(&input_path).unwrap();
    writeln!(input, "chr1\t100\t200\tr1\t0\t+").unwrap();
    writeln!(input, "chr1\t101\t201\tr2\t0\t+").unwrap();
    writeln!(input, "chr1\t120\t180\tr3\t0\t+").unwrap();
    writeln!(input, "chr1\t300\t400\tr4\t0\t+").unwrap();
    writeln!(input, "chr1\t301\t401\tr5\t0\t+").unwrap();
    writeln!(input, "chr1\t320\t360\tr6\t0\t+").unwrap();
    writeln!(input, "chr1\t100\t200\tr7\t0\t-").unwrap();

    let exe = env!("CARGO_BIN_EXE_trackclustertu");
    let output = Command::new(exe)
        .args([
            "cluster",
            "--in",
            input_path.to_str().unwrap(),
            "--format",
            "bed6",
            "--score1-threshold",
            "0.95",
            "--score2-threshold",
            "0.99",
            "--out-tu",
            out_tu.to_str().unwrap(),
            "--out-membership",
            out_membership.to_str().unwrap(),
            "--out-tu-count",
            out_tu_count.to_str().unwrap(),
        ])
        .output()
        .unwrap();

    assert!(
        output.status.success(),
        "stdout:\n{}\nstderr:\n{}",
        String::from_utf8_lossy(&output.stdout),
        String::from_utf8_lossy(&output.stderr)
    );

    let tu_text = fs::read_to_string(&out_tu).unwrap();
    assert_eq!(
        tu_text,
        concat!(
            "chr1\t100\t200\tTU000001\t0\t+\n",
            "chr1\t120\t180\tTU000002\t0\t+\n",
            "chr1\t300\t400\tTU000003\t0\t+\n",
            "chr1\t320\t360\tTU000004\t0\t+\n",
            "chr1\t100\t200\tTU000005\t0\t-\n",
        )
    );

    let membership_text = fs::read_to_string(&out_membership).unwrap();
    assert_eq!(
        membership_text,
        concat!(
            "r1\tTU000001\t1.000000\t1.000000\n",
            "r2\tTU000001\t0.980198\t0.990000\n",
            "r3\tTU000002\t1.000000\t1.000000\n",
            "r4\tTU000003\t1.000000\t1.000000\n",
            "r5\tTU000003\t0.980198\t0.990000\n",
            "r6\tTU000004\t1.000000\t1.000000\n",
            "r7\tTU000005\t1.000000\t1.000000\n",
        )
    );

    let count_text = fs::read_to_string(&out_tu_count).unwrap();
    assert_eq!(
        count_text,
        concat!(
            "tu_id,count\n",
            "TU000001,2\n",
            "TU000002,1\n",
            "TU000003,2\n",
            "TU000004,1\n",
            "TU000005,1\n",
        )
    );

    let _ = fs::remove_dir_all(&tmp);
}

#[test]
fn trackclustertu_can_skip_score2_attachment() {
    let tmp = unique_tmp_dir("trackclustertu_test_skip_score2");
    fs::create_dir_all(&tmp).unwrap();

    let input_path = tmp.join("reads.bed");
    let out_tu = tmp.join("tu.bed");
    let out_membership = tmp.join("membership.tsv");
    let out_tu_count = tmp.join("tu_count.csv");

    let mut input = fs::File::create(&input_path).unwrap();
    writeln!(input, "chr1\t100\t200\tr1\t0\t+").unwrap();
    writeln!(input, "chr1\t101\t201\tr2\t0\t+").unwrap();
    writeln!(input, "chr1\t108\t200\tr3\t0\t+").unwrap();

    let exe = env!("CARGO_BIN_EXE_trackclustertu");
    let output = Command::new(exe)
        .args([
            "cluster",
            "--in",
            input_path.to_str().unwrap(),
            "--format",
            "bed6",
            "--score1-threshold",
            "0.95",
            "--score2-threshold",
            "0.90",
            "--skip-score2-attachment",
            "--out-tu",
            out_tu.to_str().unwrap(),
            "--out-membership",
            out_membership.to_str().unwrap(),
            "--out-tu-count",
            out_tu_count.to_str().unwrap(),
        ])
        .output()
        .unwrap();

    assert!(
        output.status.success(),
        "stdout:\n{}\nstderr:\n{}",
        String::from_utf8_lossy(&output.stdout),
        String::from_utf8_lossy(&output.stderr)
    );

    let tu_text = fs::read_to_string(&out_tu).unwrap();
    assert_eq!(
        tu_text,
        concat!(
            "chr1\t100\t200\tTU000001\t0\t+\n",
            "chr1\t108\t200\tTU000002\t0\t+\n",
        )
    );

    let membership_text = fs::read_to_string(&out_membership).unwrap();
    assert_eq!(
        membership_text,
        concat!(
            "r1\tTU000001\t1.000000\t1.000000\n",
            "r2\tTU000001\t0.980198\t0.990000\n",
            "r3\tTU000002\t1.000000\t1.000000\n",
        )
    );

    let count_text = fs::read_to_string(&out_tu_count).unwrap();
    assert_eq!(
        count_text,
        concat!("tu_id,count\n", "TU000001,2\n", "TU000002,1\n",)
    );

    let _ = fs::remove_dir_all(&tmp);
}

#[test]
fn trackclustertu_cluster_can_override_three_prime_tolerance() {
    let tmp = unique_tmp_dir("trackclustertu_test_three_prime_tolerance");
    fs::create_dir_all(&tmp).unwrap();

    let input_path = tmp.join("reads.bed");
    let out_tu = tmp.join("tu.bed");

    let mut input = fs::File::create(&input_path).unwrap();
    writeln!(input, "chr1\t100\t205\tparent\t0\t+").unwrap();
    writeln!(input, "chr1\t112\t210\tchild\t0\t+").unwrap();

    let exe = env!("CARGO_BIN_EXE_trackclustertu");
    let output = Command::new(exe)
        .args([
            "cluster",
            "--in",
            input_path.to_str().unwrap(),
            "--format",
            "bed6",
            "--score1-threshold",
            "0.95",
            "--score2-threshold",
            "0.60",
            "--three-prime-tolerance-bp",
            "0",
            "--out-tu",
            out_tu.to_str().unwrap(),
        ])
        .output()
        .unwrap();

    assert!(
        output.status.success(),
        "stdout:\n{}\nstderr:\n{}",
        String::from_utf8_lossy(&output.stdout),
        String::from_utf8_lossy(&output.stderr)
    );

    let tu_text = fs::read_to_string(&out_tu).unwrap();
    assert_eq!(
        tu_text,
        concat!(
            "chr1\t100\t205\tTU000001\t0\t+\n",
            "chr1\t112\t210\tTU000002\t0\t+\n",
        )
    );

    let _ = fs::remove_dir_all(&tmp);
}

#[test]
fn trackclustertu_cluster_can_override_five_prime_delta() {
    let tmp = unique_tmp_dir("trackclustertu_test_five_prime_delta");
    fs::create_dir_all(&tmp).unwrap();

    let input_path = tmp.join("reads.bed");
    let out_tu = tmp.join("tu.bed");

    let mut input = fs::File::create(&input_path).unwrap();
    writeln!(input, "chr1\t100\t205\tparent\t0\t+").unwrap();
    writeln!(input, "chr1\t150\t210\tchild\t0\t+").unwrap();

    let exe = env!("CARGO_BIN_EXE_trackclustertu");
    let output = Command::new(exe)
        .args([
            "cluster",
            "--in",
            input_path.to_str().unwrap(),
            "--format",
            "bed6",
            "--score1-threshold",
            "0.95",
            "--score2-threshold",
            "0.60",
            "--max-5p-delta",
            "50",
            "--out-tu",
            out_tu.to_str().unwrap(),
        ])
        .output()
        .unwrap();

    assert!(
        output.status.success(),
        "stdout:\n{}\nstderr:\n{}",
        String::from_utf8_lossy(&output.stdout),
        String::from_utf8_lossy(&output.stderr)
    );

    let tu_text = fs::read_to_string(&out_tu).unwrap();
    assert_eq!(tu_text, "chr1\t100\t205\tTU000001\t0\t+\n");

    let _ = fs::remove_dir_all(&tmp);
}

#[test]
fn trackclustertu_cluster_help_hides_recount_only_and_fastq_options() {
    let exe = env!("CARGO_BIN_EXE_trackclustertu");
    let output = Command::new(exe)
        .args(["cluster", "--help"])
        .output()
        .unwrap();

    assert!(
        output.status.success(),
        "stdout:\n{}\nstderr:\n{}",
        String::from_utf8_lossy(&output.stdout),
        String::from_utf8_lossy(&output.stderr)
    );

    let help = String::from_utf8_lossy(&output.stdout);
    assert!(help.contains("Cluster bacterial directRNA reads"));
    assert!(help.contains("--three-prime-tolerance-bp"));
    assert!(help.contains("--max-5p-delta"));
    assert!(!help.contains("--pooled-membership"));
    assert!(!help.contains("fastq"));
}

#[test]
fn trackclustertu_recount_help_is_count_only() {
    let exe = env!("CARGO_BIN_EXE_trackclustertu");
    let output = Command::new(exe)
        .args(["recount", "--help"])
        .output()
        .unwrap();

    assert!(
        output.status.success(),
        "stdout:\n{}\nstderr:\n{}",
        String::from_utf8_lossy(&output.stdout),
        String::from_utf8_lossy(&output.stderr)
    );

    let help = String::from_utf8_lossy(&output.stdout);
    assert!(help.contains("Recompute TU count tables"));
    assert!(help.contains("--pooled-membership"));
    assert!(help.contains("--out-tu-count"));
    assert!(!help.contains("--in <"));
    assert!(!help.contains("--annotation-bed"));
    assert!(!help.contains("--score1-threshold"));
}
