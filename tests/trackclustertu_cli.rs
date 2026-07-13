use std::fs;
use std::io::Write;
use std::process::Command;
use std::time::{SystemTime, UNIX_EPOCH};

mod support;
use support::{count_v1_projection, membership_v1_projection};

fn unique_tmp_dir(prefix: &str) -> std::path::PathBuf {
    let nanos = SystemTime::now()
        .duration_since(UNIX_EPOCH)
        .expect("clock")
        .as_nanos();
    std::env::temp_dir().join(format!("{prefix}_{nanos}"))
}

#[test]
fn cluster_rejects_duplicate_output_paths_without_truncating_existing_file() {
    let tmp = unique_tmp_dir("trackclustertu_duplicate_outputs");
    fs::create_dir_all(&tmp).unwrap();
    let input_path = tmp.join("reads.bed");
    let shared_output = tmp.join("shared.tsv");
    fs::write(&input_path, "chr1\t10\t20\tr1\t0\t+\n").unwrap();
    fs::write(&shared_output, "previous-success\n").unwrap();

    let output = Command::new(env!("CARGO_BIN_EXE_trackclustertu"))
        .args([
            "cluster",
            "--in",
            input_path.to_str().unwrap(),
            "--format",
            "bed6",
            "--out-tu",
            shared_output.to_str().unwrap(),
            "--out-membership",
            shared_output.to_str().unwrap(),
        ])
        .output()
        .unwrap();

    assert!(!output.status.success());
    assert!(
        String::from_utf8_lossy(&output.stderr).contains("alias the same file"),
        "stderr:\n{}",
        String::from_utf8_lossy(&output.stderr)
    );
    assert_eq!(
        fs::read_to_string(&shared_output).unwrap(),
        "previous-success\n"
    );
    let _ = fs::remove_dir_all(tmp);
}

#[test]
fn cluster_rejects_output_aliasing_input_without_truncating_input() {
    let tmp = unique_tmp_dir("trackclustertu_input_output_alias");
    fs::create_dir_all(&tmp).unwrap();
    let input_path = tmp.join("reads.bed");
    let membership = tmp.join("membership.tsv");
    let original = "chr1\t10\t20\tr1\t0\t+\n";
    fs::write(&input_path, original).unwrap();

    let output = Command::new(env!("CARGO_BIN_EXE_trackclustertu"))
        .args([
            "cluster",
            "--in",
            input_path.to_str().unwrap(),
            "--format",
            "bed6",
            "--out-tu",
            input_path.to_str().unwrap(),
            "--out-membership",
            membership.to_str().unwrap(),
        ])
        .output()
        .unwrap();

    assert!(!output.status.success());
    assert!(
        String::from_utf8_lossy(&output.stderr).contains("aliases input"),
        "stderr:\n{}",
        String::from_utf8_lossy(&output.stderr)
    );
    assert_eq!(fs::read_to_string(&input_path).unwrap(), original);
    assert!(!membership.exists());
    let _ = fs::remove_dir_all(tmp);
}

#[test]
fn cluster_skips_unknown_strand_records_it_and_continues() {
    let tmp = unique_tmp_dir("trackclustertu_unknown_strand");
    fs::create_dir_all(&tmp).unwrap();
    let input_path = tmp.join("reads.bed");
    let out_dir = tmp.join("results");
    fs::write(
        &input_path,
        concat!(
            "chr1\t0\t10\tgood_before\t0\t+\n",
            "chr1\t10\t20\tunknown_read\t0\t.\n",
            "chr1\t20\t30\tgood_after\t0\t+\n",
        ),
    )
    .unwrap();

    let output = Command::new(env!("CARGO_BIN_EXE_trackclustertu"))
        .args([
            "cluster",
            "--in",
            input_path.to_str().unwrap(),
            "--format",
            "bed6",
            "--out-dir",
            out_dir.to_str().unwrap(),
        ])
        .output()
        .unwrap();

    assert!(
        output.status.success(),
        "stderr:\n{}",
        String::from_utf8_lossy(&output.stderr)
    );
    let stderr = String::from_utf8_lossy(&output.stderr);
    assert!(
        stderr.contains("read_input_counts\ttotal=3\tretained=2\trejected=1"),
        "stderr:\n{stderr}"
    );
    assert!(stderr.contains("unknown_strand=1"), "stderr:\n{stderr}");

    let membership = fs::read_to_string(out_dir.join("membership.tsv")).unwrap();
    assert!(membership.contains("good_before"), "{membership}");
    assert!(membership.contains("good_after"), "{membership}");
    assert!(!membership.contains("unknown_read"), "{membership}");

    let rejections = fs::read_to_string(out_dir.join("read_rejections.tsv")).unwrap();
    assert!(
        rejections.contains(input_path.to_str().unwrap()),
        "{rejections}"
    );
    assert!(
        rejections.contains("\t2\tunknown_read\tvalidation\tunknown_strand\t"),
        "{rejections}"
    );
    let _ = fs::remove_dir_all(tmp);
}

#[test]
fn strict_read_errors_restores_transactional_failure() {
    let tmp = unique_tmp_dir("trackclustertu_strict_read_errors");
    fs::create_dir_all(&tmp).unwrap();
    let input_path = tmp.join("reads.bed");
    let out_dir = tmp.join("results");
    fs::write(
        &input_path,
        "chr1\t0\t10\tgood\t0\t+\nchr1\t10\t20\tbad\t0\t.\n",
    )
    .unwrap();

    let output = Command::new(env!("CARGO_BIN_EXE_trackclustertu"))
        .args([
            "cluster",
            "--in",
            input_path.to_str().unwrap(),
            "--format",
            "bed6",
            "--out-dir",
            out_dir.to_str().unwrap(),
            "--strict-read-errors",
        ])
        .output()
        .unwrap();

    assert!(!output.status.success());
    let stderr = String::from_utf8_lossy(&output.stderr);
    assert!(stderr.contains("strict read"), "stderr:\n{stderr}");
    assert!(stderr.contains("unknown_strand"), "stderr:\n{stderr}");
    assert!(!out_dir.join("tus.bed").exists());
    assert!(!out_dir.join("read_rejections.tsv").exists());
    let _ = fs::remove_dir_all(tmp);
}

#[test]
fn all_invalid_reads_publish_empty_outputs_and_complete_rejection_report() {
    let tmp = unique_tmp_dir("trackclustertu_all_invalid_reads");
    fs::create_dir_all(&tmp).unwrap();
    let input_path = tmp.join("reads.bed");
    let out_dir = tmp.join("results");
    fs::write(&input_path, "malformed\nchr1\t10\t20\tunknown\t0\t.\n").unwrap();

    let output = Command::new(env!("CARGO_BIN_EXE_trackclustertu"))
        .args([
            "cluster",
            "--in",
            input_path.to_str().unwrap(),
            "--format",
            "bed6",
            "--out-dir",
            out_dir.to_str().unwrap(),
        ])
        .output()
        .unwrap();

    assert!(
        output.status.success(),
        "stderr:\n{}",
        String::from_utf8_lossy(&output.stderr)
    );
    assert_eq!(fs::read_to_string(out_dir.join("tus.bed")).unwrap(), "");
    let rejections = fs::read_to_string(out_dir.join("read_rejections.tsv")).unwrap();
    assert_eq!(rejections.matches("\ttoo_few_columns\t").count(), 1);
    assert_eq!(rejections.matches("\tunknown_strand\t").count(), 1);
    let stderr = String::from_utf8_lossy(&output.stderr);
    assert!(
        stderr.contains("total=2\tretained=0\trejected=2"),
        "stderr:\n{stderr}"
    );
    assert!(stderr.contains("every input read was rejected"), "{stderr}");
    let _ = fs::remove_dir_all(tmp);
}

#[test]
fn late_annotation_error_preserves_prior_outputs_and_publishes_no_partial_set() {
    let tmp = unique_tmp_dir("trackclustertu_atomic_annotation_failure");
    fs::create_dir_all(&tmp).unwrap();
    let input_path = tmp.join("reads.bed");
    let annotation_path = tmp.join("genes.bed");
    let out_dir = tmp.join("nested/results");
    fs::write(&input_path, "chr1\t10\t20\tr1\t0\t+\n").unwrap();
    fs::write(&annotation_path, "malformed\n").unwrap();
    fs::create_dir_all(&out_dir).unwrap();
    fs::write(out_dir.join("tus.bed"), "previous-tu-output\n").unwrap();
    fs::write(
        out_dir.join("membership.tsv"),
        "previous-membership-output\n",
    )
    .unwrap();

    let output = Command::new(env!("CARGO_BIN_EXE_trackclustertu"))
        .args([
            "cluster",
            "--in",
            input_path.to_str().unwrap(),
            "--format",
            "bed6",
            "--annotation-bed",
            annotation_path.to_str().unwrap(),
            "--out-dir",
            out_dir.to_str().unwrap(),
        ])
        .output()
        .unwrap();

    assert!(!output.status.success());
    assert_eq!(
        fs::read_to_string(out_dir.join("tus.bed")).unwrap(),
        "previous-tu-output\n"
    );
    assert_eq!(
        fs::read_to_string(out_dir.join("membership.tsv")).unwrap(),
        "previous-membership-output\n"
    );
    for name in [
        "tu_endpoint_stats.tsv",
        "tu_id_map.tsv",
        "tu_count.csv",
        "read_rejections.tsv",
        "tu_gene.tsv",
        "tu_semantics.tsv",
        "tus.gff3",
        "tus.anchored.bed12",
        "gene_count.csv",
    ] {
        assert!(!out_dir.join(name).exists(), "unexpected output {name}");
    }
    assert!(fs::read_dir(&out_dir).unwrap().all(|entry| !entry
        .unwrap()
        .file_name()
        .to_string_lossy()
        .contains("staged")));
    let _ = fs::remove_dir_all(tmp);
}

#[test]
fn cluster_creates_nested_output_directories() {
    let tmp = unique_tmp_dir("trackclustertu_nested_outputs");
    fs::create_dir_all(&tmp).unwrap();
    let input_path = tmp.join("reads.bed");
    let out_tu = tmp.join("a/b/c/tus.bed");
    let out_membership = tmp.join("d/e/membership.tsv");
    fs::write(&input_path, "chr1\t10\t20\tr1\t0\t+\n").unwrap();

    let output = Command::new(env!("CARGO_BIN_EXE_trackclustertu"))
        .args([
            "cluster",
            "--tu-id-style",
            "sequential",
            "--in",
            input_path.to_str().unwrap(),
            "--format",
            "bed6",
            "--out-tu",
            out_tu.to_str().unwrap(),
            "--out-membership",
            out_membership.to_str().unwrap(),
        ])
        .output()
        .unwrap();

    assert!(
        output.status.success(),
        "stdout:\n{}\nstderr:\n{}",
        String::from_utf8_lossy(&output.stdout),
        String::from_utf8_lossy(&output.stderr)
    );
    assert!(out_tu.exists());
    assert!(out_membership.exists());
    let _ = fs::remove_dir_all(tmp);
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
    // Deliberately exercise the deprecated aliases once for script compatibility. Routine tests
    // below use the canonical metric names.
    let output = Command::new(exe)
        .args([
            "cluster",
            "--tu-id-style",
            "sequential",
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
        membership_v1_projection(&membership_text),
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
        count_v1_projection(&count_text),
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
            "--tu-id-style",
            "sequential",
            "--in",
            input_path.to_str().unwrap(),
            "--format",
            "bed6",
            "--span-jaccard-threshold",
            "0.95",
            "--overlap-over-longer-threshold",
            "0.90",
            "--skip-overlap-over-longer-attachment",
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
        membership_v1_projection(&membership_text),
        concat!(
            "r1\tTU000001\t1.000000\t1.000000\n",
            "r2\tTU000001\t0.980198\t0.990000\n",
            "r3\tTU000002\t1.000000\t1.000000\n",
        )
    );

    let count_text = fs::read_to_string(&out_tu_count).unwrap();
    assert_eq!(
        count_v1_projection(&count_text),
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
            "--tu-id-style",
            "sequential",
            "--in",
            input_path.to_str().unwrap(),
            "--format",
            "bed6",
            "--span-jaccard-threshold",
            "0.95",
            "--overlap-over-longer-threshold",
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
            "--tu-id-style",
            "sequential",
            "--in",
            input_path.to_str().unwrap(),
            "--format",
            "bed6",
            "--span-jaccard-threshold",
            "0.95",
            "--overlap-over-longer-threshold",
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
    let out_tsv = tmp.join("diagnose/reports/missed.tsv");
    let out_bed = tmp.join("diagnose/tracks/missed.bed");

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
    let rejection_report = out_tsv
        .parent()
        .unwrap()
        .join("diagnose_read_rejections.tsv");
    let rejection_text = fs::read_to_string(rejection_report).unwrap();
    assert!(rejection_text.contains("#trackclustertu_read_rejections_schema="));
    assert_eq!(rejection_text.lines().count(), 2);

    let _ = fs::remove_dir_all(&tmp);
}

#[test]
fn diagnose_skips_bad_bed_reads_and_writes_overridden_rejection_report() {
    let tmp = unique_tmp_dir("trackclustertu_diagnose_tolerant_reads");
    fs::create_dir_all(&tmp).unwrap();

    let reads_path = tmp.join("reads.bed");
    let tus_path = tmp.join("tus.bed");
    let out_tsv = tmp.join("diagnose/missed.tsv");
    let rejection_path = tmp.join("audit/custom_rejected_reads.tsv");
    fs::write(
        &reads_path,
        concat!(
            "chr1\t100\t210\tr1\t0\t+\n",
            "chr1\tbad\t210\tr1\t0\t+\n",
            "chr1\t102\t210\tr2\t0\t+\n",
            "chr1\t103\t210\tr1\t0\t.\n",
            "chr1\t104\t210\tr3\t0\t+\n",
        ),
    )
    .unwrap();
    fs::write(&tus_path, "chr1\t150\t210\tTU000001\t0\t+\n").unwrap();

    let exe = env!("CARGO_BIN_EXE_trackclustertu");
    let output = Command::new(exe)
        .args([
            "diagnose-missed-tus",
            "--in",
            reads_path.to_str().unwrap(),
            "--existing-tu",
            tus_path.to_str().unwrap(),
            "--out-tsv",
            out_tsv.to_str().unwrap(),
            "--out-read-rejections",
            rejection_path.to_str().unwrap(),
            "--min-family-support",
            "2",
            "--min-mode-support",
            "2",
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
    let rejections = fs::read_to_string(&rejection_path).unwrap();
    assert!(rejections.contains("\t2\t.\tparse\tinvalid_integer\t"));
    assert!(rejections.contains("\t4\tr1\tvalidate\tunknown_strand\t"));
    assert!(!rejections.contains("\tduplicate_read_id\t"));
    let stderr = String::from_utf8_lossy(&output.stderr);
    assert!(
        stderr.contains("total=5\tretained=3\trejected=2"),
        "stderr:\n{stderr}"
    );
    let stdout = String::from_utf8_lossy(&output.stdout);
    assert!(
        stdout.contains(rejection_path.to_string_lossy().as_ref()),
        "stdout:\n{stdout}"
    );

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
    let out_tu = tmp.join("rescue/tu/rescued.tus.bed");
    let out_membership = tmp.join("rescue/membership/rescued.membership.tsv");
    let out_tu_count = tmp.join("rescue/counts/rescued.tu_count.csv");
    let out_candidates_tsv = tmp.join("rescue/reports/rescued.candidates.tsv");
    let out_candidates_bed = tmp.join("rescue/reports/rescued.candidates.bed");

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
            "--out-candidates-bed",
            out_candidates_bed.to_str().unwrap(),
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
    assert!(
        candidate_report.starts_with("candidate_id\trescued_tu_id\trescued_hard_count\t"),
        "{candidate_report}"
    );
    assert!(candidate_report.contains("MISS0001"), "{candidate_report}");
    assert!(
        candidate_report.contains("MISS0001\tRESC0001\t3\t"),
        "{candidate_report}"
    );
    assert!(
        candidate_report.contains("same_3p_missing_5p_mode"),
        "{candidate_report}"
    );
    assert!(candidate_report.contains("geneA"), "{candidate_report}");
    let candidate_bed = fs::read_to_string(&out_candidates_bed).unwrap();
    assert!(
        candidate_bed.contains("\tMISS0001|RESC0001|same_3p_missing_5p_mode|geneA\t"),
        "{candidate_bed}"
    );
    let rejection_report = out_tu.parent().unwrap().join("rescue_read_rejections.tsv");
    let rejection_text = fs::read_to_string(rejection_report).unwrap();
    assert!(rejection_text.contains("#trackclustertu_read_rejections_schema="));
    assert_eq!(rejection_text.lines().count(), 2);

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
    assert!(help.contains("--span-jaccard-threshold"));
    assert!(help.contains("--score1-threshold"));
    assert!(help.contains("--overlap-over-longer-threshold"));
    assert!(help.contains("--score2-threshold"));
    assert!(help.contains("--skip-overlap-over-longer-attachment"));
    assert!(help.contains("--skip-score2-attachment"));
    assert!(help.contains("--ambiguity-margin"));
    assert!(help.contains("--fractional-assignment"));
    assert!(help.contains("--out-tu-endpoint-stats"));
    assert!(help.contains("--tu-id-style"));
    assert!(help.contains("--out-tu-id-map"));
    assert!(help.contains("--gene-min-overlap-bp"));
    assert!(help.contains("--gene-min-tu-fraction"));
    assert!(help.contains("--gene-min-gene-fraction"));
    assert!(help.contains("--out-tu-semantics"));
    assert!(help.contains("--out-tu-gff3"));
    assert!(!help.contains("--pooled-membership"));
    assert!(!help.contains("fastq"));
}

#[test]
fn trackclustertu_run_help_shows_canonical_score_flags_and_compatibility_aliases() {
    let output = Command::new(env!("CARGO_BIN_EXE_trackclustertu"))
        .args(["run", "--help"])
        .output()
        .unwrap();
    assert!(
        output.status.success(),
        "stderr:\n{}",
        String::from_utf8_lossy(&output.stderr)
    );
    let help = String::from_utf8_lossy(&output.stdout);
    assert!(help.contains("--span-jaccard-threshold"));
    assert!(help.contains("--overlap-over-longer-threshold"));
    assert!(help.contains("--score1-threshold"));
    assert!(help.contains("--score2-threshold"));
    assert!(help.contains("--skip-overlap-over-longer-attachment"));
    assert!(help.contains("--skip-score2-attachment"));
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
    assert!(!help.contains("--span-jaccard-threshold"));
    assert!(!help.contains("--score1-threshold"));
}

#[test]
fn trackclustertu_rescue_help_exposes_documented_defaults() {
    let output = Command::new(env!("CARGO_BIN_EXE_trackclustertu"))
        .args(["rescue-missed-tus", "--help"])
        .output()
        .unwrap();

    assert!(
        output.status.success(),
        "stdout:\n{}\nstderr:\n{}",
        String::from_utf8_lossy(&output.stdout),
        String::from_utf8_lossy(&output.stderr)
    );

    let help = String::from_utf8_lossy(&output.stdout);
    for option in [
        "--rescue-prefix",
        "--three-prime-window-bp",
        "--max-three-prime-family-diameter-bp",
        "--five-prime-window-bp",
        "--min-family-support",
        "--min-mode-support",
        "--min-mode-fraction",
        "--max-candidates-per-family",
    ] {
        assert!(help.contains(option), "missing {option} in help:\n{help}");
    }
    assert!(help.contains("[default: RESC]"), "help:\n{help}");
    assert!(help.contains("[default: 12]"), "help:\n{help}");
    assert!(help.contains("[default: 10]"), "help:\n{help}");
    assert!(help.contains("[default: 20]"), "help:\n{help}");
    assert!(help.contains("[default: 0.02]"), "help:\n{help}");
    assert!(help.contains("[default: 3]"), "help:\n{help}");
    assert!(
        help.contains("Defaults to `--three-prime-window-bp`"),
        "help:\n{help}"
    );
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

#[test]
fn trackclustertu_top_level_version_includes_binary_name() {
    let output = Command::new(env!("CARGO_BIN_EXE_trackclustertu"))
        .arg("--version")
        .output()
        .unwrap();

    assert!(
        output.status.success(),
        "stderr:\n{}",
        String::from_utf8_lossy(&output.stderr)
    );
    assert_eq!(
        String::from_utf8_lossy(&output.stdout),
        format!("trackclustertu {}\n", env!("CARGO_PKG_VERSION"))
    );
}
