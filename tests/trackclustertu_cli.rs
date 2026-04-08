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
fn trackclustertu_diagnose_missed_tus_reports_same_3p_missing_5p_mode() {
    let tmp = unique_tmp_dir("trackclustertu_diagnose_missed_tus");
    fs::create_dir_all(&tmp).unwrap();

    let reads_path = tmp.join("reads.bed");
    let tus_path = tmp.join("tus.bed");
    let genes_path = tmp.join("genes.bed");
    let out_tsv = tmp.join("missed.tsv");
    let out_bed = tmp.join("missed.bed");

    let mut reads = fs::File::create(&reads_path).unwrap();
    writeln!(reads, "chr1\t100\t210\tr1\t0\t+").unwrap();
    writeln!(reads, "chr1\t102\t210\tr2\t0\t+").unwrap();
    writeln!(reads, "chr1\t104\t210\tr3\t0\t+").unwrap();
    writeln!(reads, "chr1\t150\t210\tr4\t0\t+").unwrap();
    writeln!(reads, "chr1\t151\t210\tr5\t0\t+").unwrap();

    let mut tus = fs::File::create(&tus_path).unwrap();
    writeln!(tus, "chr1\t150\t210\tTU000001\t0\t+").unwrap();

    let mut genes = fs::File::create(&genes_path).unwrap();
    writeln!(genes, "chr1\t105\t205\tgeneA\t0\t+").unwrap();

    let exe = env!("CARGO_BIN_EXE_trackclustertu");
    let output = Command::new(exe)
        .args([
            "diagnose-missed-tus",
            "--in",
            reads_path.to_str().unwrap(),
            "--existing-tu",
            tus_path.to_str().unwrap(),
            "--annotation-bed",
            genes_path.to_str().unwrap(),
            "--out-tsv",
            out_tsv.to_str().unwrap(),
            "--out-bed",
            out_bed.to_str().unwrap(),
            "--three-prime-window-bp",
            "12",
            "--five-prime-window-bp",
            "10",
            "--min-family-support",
            "2",
            "--min-mode-support",
            "2",
            "--min-mode-fraction",
            "0.20",
        ])
        .output()
        .unwrap();

    assert!(
        output.status.success(),
        "stdout:\n{}\nstderr:\n{}",
        String::from_utf8_lossy(&output.stdout),
        String::from_utf8_lossy(&output.stderr)
    );

    let report = fs::read_to_string(&out_tsv).unwrap();
    assert!(report.contains("same_3p_missing_5p_mode"), "{report}");
    assert!(
        report.contains("\tchr1\t+\t100\t210\t210\t100\t5\t3\t0.600000\t"),
        "{report}"
    );
    assert!(
        report.contains("\tTU000001\t150\t210\t0\t50\t0.545455\t0.545455\tgeneA"),
        "{report}"
    );

    let bed = fs::read_to_string(&out_bed).unwrap();
    assert!(bed.contains("chr1\t100\t210\tMISS"), "{bed}");

    let _ = fs::remove_dir_all(&tmp);
}

#[test]
fn trackclustertu_rescue_missed_tus_promotes_boundary_mode() {
    let tmp = unique_tmp_dir("trackclustertu_rescue_missed_tus");
    fs::create_dir_all(&tmp).unwrap();

    let reads_path = tmp.join("reads.bed");
    let tus_path = tmp.join("tus.bed");
    let membership_path = tmp.join("membership.tsv");
    let genes_path = tmp.join("genes.bed");
    let out_tu = tmp.join("rescued.tus.bed");
    let out_membership = tmp.join("rescued.membership.tsv");
    let out_tu_count = tmp.join("rescued.tu_count.csv");
    let out_candidates_tsv = tmp.join("rescued.candidates.tsv");

    let mut reads = fs::File::create(&reads_path).unwrap();
    writeln!(reads, "chr1\t100\t210\tr1\t0\t+").unwrap();
    writeln!(reads, "chr1\t102\t210\tr2\t0\t+").unwrap();
    writeln!(reads, "chr1\t104\t210\tr3\t0\t+").unwrap();
    writeln!(reads, "chr1\t150\t210\tr4\t0\t+").unwrap();
    writeln!(reads, "chr1\t151\t210\tr5\t0\t+").unwrap();

    let mut tus = fs::File::create(&tus_path).unwrap();
    writeln!(tus, "chr1\t150\t210\tTU000001\t0\t+").unwrap();

    let mut membership = fs::File::create(&membership_path).unwrap();
    writeln!(membership, "r1\tTU000001\t0.545455\t0.545455").unwrap();
    writeln!(membership, "r2\tTU000001\t0.522523\t0.522523").unwrap();
    writeln!(membership, "r3\tTU000001\t0.500000\t0.500000").unwrap();
    writeln!(membership, "r4\tTU000001\t1.000000\t1.000000").unwrap();
    writeln!(membership, "r5\tTU000001\t0.983333\t0.983333").unwrap();

    let mut genes = fs::File::create(&genes_path).unwrap();
    writeln!(genes, "chr1\t105\t205\tgeneA\t0\t+").unwrap();

    let exe = env!("CARGO_BIN_EXE_trackclustertu");
    let output = Command::new(exe)
        .args([
            "rescue-missed-tus",
            "--in",
            reads_path.to_str().unwrap(),
            "--existing-tu",
            tus_path.to_str().unwrap(),
            "--existing-membership",
            membership_path.to_str().unwrap(),
            "--annotation-bed",
            genes_path.to_str().unwrap(),
            "--out-tu",
            out_tu.to_str().unwrap(),
            "--out-membership",
            out_membership.to_str().unwrap(),
            "--out-tu-count",
            out_tu_count.to_str().unwrap(),
            "--out-candidates-tsv",
            out_candidates_tsv.to_str().unwrap(),
            "--three-prime-window-bp",
            "12",
            "--five-prime-window-bp",
            "10",
            "--min-family-support",
            "2",
            "--min-mode-support",
            "2",
            "--min-mode-fraction",
            "0.20",
        ])
        .output()
        .unwrap();

    assert!(
        output.status.success(),
        "stdout:\n{}\nstderr:\n{}",
        String::from_utf8_lossy(&output.stdout),
        String::from_utf8_lossy(&output.stderr)
    );

    let rescued_tus = fs::read_to_string(&out_tu).unwrap();
    assert_eq!(
        rescued_tus,
        concat!(
            "chr1\t100\t210\tRESC0001\t0\t+\n",
            "chr1\t150\t210\tTU000001\t0\t+\n",
        )
    );

    let rescued_membership = fs::read_to_string(&out_membership).unwrap();
    assert_eq!(
        rescued_membership,
        concat!(
            "r1\tRESC0001\t1.000000\t1.000000\n",
            "r2\tRESC0001\t0.981818\t0.981818\n",
            "r3\tRESC0001\t0.963636\t0.963636\n",
            "r4\tTU000001\t1.000000\t1.000000\n",
            "r5\tTU000001\t0.983333\t0.983333\n",
        )
    );

    let rescued_counts = fs::read_to_string(&out_tu_count).unwrap();
    assert_eq!(
        rescued_counts,
        concat!("tu_id,count\n", "RESC0001,3\n", "TU000001,2\n",)
    );

    let candidate_report = fs::read_to_string(&out_candidates_tsv).unwrap();
    assert!(candidate_report.contains("MISS0001"), "{candidate_report}");
    assert!(
        candidate_report.contains("same_3p_missing_5p_mode"),
        "{candidate_report}"
    );
    assert!(candidate_report.contains("geneA"), "{candidate_report}");

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

#[test]
fn trackclustertu_top_level_help_lists_diagnose_missed_tus() {
    let exe = env!("CARGO_BIN_EXE_trackclustertu");
    let output = Command::new(exe).args(["--help"]).output().unwrap();

    assert!(
        output.status.success(),
        "stdout:\n{}\nstderr:\n{}",
        String::from_utf8_lossy(&output.stdout),
        String::from_utf8_lossy(&output.stderr)
    );

    let help = String::from_utf8_lossy(&output.stdout);
    assert!(help.contains("diagnose-missed-tus"));
    assert!(help.contains("rescue-missed-tus"));
}
