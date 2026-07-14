use std::collections::{HashMap, HashSet};
use std::fs;
use std::path::{Path, PathBuf};
use std::process::{Command, Output};
use std::sync::atomic::{AtomicU64, Ordering};
use std::time::{SystemTime, UNIX_EPOCH};

static NEXT_TEMP_ID: AtomicU64 = AtomicU64::new(0);

struct Fixture {
    root: PathBuf,
    reads: PathBuf,
    tus: PathBuf,
    membership: PathBuf,
}

fn fixture(name: &str, reads: &str, tus: &str, membership: &str) -> Fixture {
    let nanos = SystemTime::now()
        .duration_since(UNIX_EPOCH)
        .expect("clock")
        .as_nanos();
    let sequence = NEXT_TEMP_ID.fetch_add(1, Ordering::Relaxed);
    let root = std::env::temp_dir().join(format!(
        "trackclustertu_{name}_{}_{nanos}_{sequence}",
        std::process::id()
    ));
    fs::create_dir_all(&root).unwrap();

    let reads_path = root.join("reads.bed");
    let tus_path = root.join("tus.bed");
    let membership_path = root.join("membership.tsv");
    fs::write(&reads_path, reads).unwrap();
    fs::write(&tus_path, tus).unwrap();
    fs::write(&membership_path, membership).unwrap();

    Fixture {
        root,
        reads: reads_path,
        tus: tus_path,
        membership: membership_path,
    }
}

fn standard_reads() -> &'static str {
    concat!(
        "chr1\t100\t210\tr1\t0\t+\n",
        "chr1\t102\t210\tr2\t0\t+\n",
        "chr1\t104\t210\tr3\t0\t+\n",
        "chr1\t150\t210\tr4\t0\t+\n",
        "chr1\t151\t210\tr5\t0\t+\n",
    )
}

fn standard_membership(tu_id: &str) -> String {
    (1..=5)
        .map(|idx| format!("r{idx}\t{tu_id}\t1.0\t1.0\n"))
        .collect()
}

fn full_v2_row(
    read_id: &str,
    hard_tu_id: &str,
    status: &str,
    best: &str,
    second: &str,
    primary_weight: &str,
    fractions: &str,
) -> String {
    full_v2_row_with_evidence(
        read_id,
        hard_tu_id,
        status,
        best,
        second,
        primary_weight,
        fractions,
        false,
    )
}

#[allow(clippy::too_many_arguments)]
fn full_v2_row_with_evidence(
    read_id: &str,
    hard_tu_id: &str,
    status: &str,
    best: &str,
    second: &str,
    primary_weight: &str,
    fractions: &str,
    full_length_evidence: bool,
) -> String {
    format!(
        "{read_id}\t{hard_tu_id}\t.\t.\tv2\t{status}\t{best}\t{second}\t.\t.\t.\t.\t.\t.\t.\t.\t.\t{primary_weight}\t{fractions}\t{full_length_evidence}\n"
    )
}

fn run_rescue(
    fixture: &Fixture,
    out_tu: &Path,
    out_membership: &Path,
    out_counts: Option<&Path>,
    rescue_prefix: Option<&str>,
) -> Output {
    let exe = env!("CARGO_BIN_EXE_trackclustertu");
    let mut command = Command::new(exe);
    command
        .arg("rescue-missed-tus")
        .arg("--in")
        .arg(&fixture.reads)
        .arg("--existing-tu")
        .arg(&fixture.tus)
        .arg("--existing-membership")
        .arg(&fixture.membership)
        .arg("--out-tu")
        .arg(out_tu)
        .arg("--out-membership")
        .arg(out_membership)
        .args([
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
        ]);
    if let Some(out_counts) = out_counts {
        command.arg("--out-tu-count").arg(out_counts);
    }
    if let Some(rescue_prefix) = rescue_prefix {
        command.arg("--rescue-prefix").arg(rescue_prefix);
    }
    command.output().unwrap()
}

fn assert_failed_with_line(output: &Output, path: &Path, line: usize, message: &str) {
    assert!(!output.status.success(), "command unexpectedly succeeded");
    let stderr = String::from_utf8_lossy(&output.stderr);
    assert!(stderr.contains(message), "stderr:\n{stderr}");
    assert!(
        stderr.contains(&format!("{path:?}:{line}")),
        "missing path and line context in stderr:\n{stderr}"
    );
}

#[test]
fn rescue_ids_skip_existing_ids_and_outputs_remain_consistent() {
    let fixture = fixture(
        "rescue_collision",
        standard_reads(),
        "chr1\t150\t210\tRESC0001\t0\t+\n",
        &standard_membership("RESC0001"),
    );
    let out_tu = fixture.root.join("rescued.tus.bed");
    let out_membership = fixture.root.join("rescued.membership.tsv");
    let out_counts = fixture.root.join("rescued.counts.csv");

    let output = run_rescue(&fixture, &out_tu, &out_membership, Some(&out_counts), None);
    assert!(
        output.status.success(),
        "stdout:\n{}\nstderr:\n{}",
        String::from_utf8_lossy(&output.stdout),
        String::from_utf8_lossy(&output.stderr)
    );

    let tu_text = fs::read_to_string(&out_tu).unwrap();
    let tu_ids: Vec<&str> = tu_text
        .lines()
        .map(|line| line.split('\t').nth(3).unwrap())
        .collect();
    assert_eq!(tu_ids, ["RESC0002", "RESC0001"]);
    assert_eq!(tu_ids.iter().copied().collect::<HashSet<_>>().len(), 2);

    let membership_text = fs::read_to_string(&out_membership).unwrap();
    let mut membership_counts: HashMap<&str, u64> = HashMap::new();
    let mut membership_rows = 0u64;
    for line in membership_text.lines() {
        let tu_id = line.split('\t').nth(1).unwrap();
        assert!(tu_ids.contains(&tu_id), "unknown membership TU ID {tu_id}");
        *membership_counts.entry(tu_id).or_default() += 1;
        membership_rows += 1;
    }

    let counts_text = fs::read_to_string(&out_counts).unwrap();
    let emitted_counts: HashMap<&str, u64> = counts_text
        .lines()
        .skip(1)
        .map(|line| {
            let (tu_id, count) = line.split_once(',').unwrap();
            (tu_id, count.parse().unwrap())
        })
        .collect();
    assert_eq!(emitted_counts.values().sum::<u64>(), membership_rows);
    assert_eq!(emitted_counts, membership_counts);

    fs::remove_dir_all(&fixture.root).unwrap();
}

#[test]
fn v2_ambiguous_and_unassigned_rows_have_no_hard_existing_tu_foreign_key() {
    let membership = concat!(
        "#trackclustertu_membership_schema=v2\n",
        "#columns=read_id\ttu_id\tscore1\tscore2\tschema_version\tassignment_status\n",
        "r1\t.\t.\t.\tv2\tambiguous\n",
        "r2\t.\t.\t.\tv2\tunassigned\n",
        "r3\t.\t.\t.\tv2\tambiguous\n",
        "r4\tTU000001\t1.0\t1.0\tv2\tunique\n",
        "r5\tTU000001\t0.98\t0.98\tv2\tpartial\n",
    );
    let fixture = fixture(
        "v2_no_hard_assignment",
        standard_reads(),
        "chr1\t150\t210\tTU000001\t0\t+\n",
        membership,
    );
    let out_tu = fixture.root.join("rescued.tus.bed");
    let out_membership = fixture.root.join("rescued.membership.tsv");

    let output = run_rescue(&fixture, &out_tu, &out_membership, None, None);
    assert!(
        output.status.success(),
        "stdout:\n{}\nstderr:\n{}",
        String::from_utf8_lossy(&output.stdout),
        String::from_utf8_lossy(&output.stderr)
    );
    let rescued_membership = fs::read_to_string(&out_membership).unwrap();
    assert!(rescued_membership.contains("#trackclustertu_membership_schema=v2\n"));
    assert!(rescued_membership.contains("#row_widths=6,20\n"));
    assert!(rescued_membership.contains(
        "#compatible_v2_columns=read_id\ttu_id\tscore1\tscore2\tschema_version\tassignment_status\n"
    ));
    assert!(rescued_membership.contains(
        "#canonical_v2_columns=read_id\ttu_id\tscore1\tscore2\tschema_version\tassignment_status\t"
    ));
    assert!(rescued_membership.contains("#contains_compatible_v2_rows=true\n"));
    assert!(!rescued_membership.contains("#columns=read_id\ttu_id\tscore1\tscore2\tschema_version\tassignment_status\tbest_candidate_tu_id"));
    assert!(rescued_membership.contains("r1\tRESC0001\t"));
    assert!(rescued_membership.contains("r2\tRESC0001\t"));
    assert!(rescued_membership.contains("r3\tRESC0001\t"));
    assert!(rescued_membership.contains("r4\tTU000001\t"));
    assert!(rescued_membership.contains("r5\tTU000001\t"));
    let row_widths: Vec<usize> = rescued_membership
        .lines()
        .filter(|line| !line.starts_with('#'))
        .map(|line| line.split('\t').count())
        .collect();
    assert_eq!(row_widths, [20, 20, 20, 6, 6]);

    let _ = fs::remove_dir_all(&fixture.root);
}

#[test]
fn declared_empty_v2_membership_promotes_missing_rows_as_canonical_v2() {
    let membership = concat!(
        "#trackclustertu_membership_schema=v2\n",
        "#columns=read_id\ttu_id\tscore1\tscore2\tschema_version\tassignment_status\n",
        "#ambiguity_margin=0.03\n",
    );
    let fixture = fixture(
        "declared_empty_v2",
        standard_reads(),
        "chr1\t150\t210\tTU000001\t0\t+\n",
        membership,
    );
    let out_tu = fixture.root.join("rescued.tus.bed");
    let out_membership = fixture.root.join("rescued.membership.tsv");

    let output = run_rescue(&fixture, &out_tu, &out_membership, None, None);
    assert!(
        output.status.success(),
        "stdout:\n{}\nstderr:\n{}",
        String::from_utf8_lossy(&output.stdout),
        String::from_utf8_lossy(&output.stderr)
    );

    let rescued_membership = fs::read_to_string(&out_membership).unwrap();
    assert!(rescued_membership.contains("#trackclustertu_membership_schema=v2\n"));
    assert!(rescued_membership.contains("#columns=read_id\ttu_id\tscore1\tscore2\tschema_version\tassignment_status\tbest_candidate_tu_id\t"));
    assert!(rescued_membership.contains("#ambiguity_margin=0.03\n"));
    let rows: Vec<Vec<&str>> = rescued_membership
        .lines()
        .filter(|line| !line.starts_with('#'))
        .map(|line| line.split('\t').collect())
        .collect();
    assert_eq!(rows.len(), 3);
    for row in rows {
        assert_eq!(row.len(), 20);
        assert_eq!(row[1], "RESC0001");
        assert_eq!(row[4], "v2");
        assert_eq!(row[5], "unique");
        assert_eq!(row[18], "RESC0001:1");
    }

    fs::remove_dir_all(&fixture.root).unwrap();
}

#[test]
fn mixed_marker_reflects_retained_output_rows_not_quarantined_input_rows() {
    let reads = concat!(
        "chr1\t100\t210\tr1\t0\t+\n",
        "chr1\tbad\t210\tbad_read\t0\t+\n",
    );
    let mut membership = full_v2_row(
        "r1",
        "TU000001",
        "unique",
        "TU000001",
        ".",
        "1",
        "TU000001:1",
    );
    membership.push_str("bad_read\tTU000001\t1\t1\n");
    let fixture = fixture(
        "actual_output_schema",
        reads,
        "chr1\t100\t210\tTU000001\t0\t+\n",
        &membership,
    );
    let out_tu = fixture.root.join("output/rescued.tus.bed");
    let out_membership = fixture.root.join("output/rescued.membership.tsv");

    let output = run_rescue(&fixture, &out_tu, &out_membership, None, None);
    assert!(
        output.status.success(),
        "stdout:\n{}\nstderr:\n{}",
        String::from_utf8_lossy(&output.stdout),
        String::from_utf8_lossy(&output.stderr)
    );

    let rescued_membership = fs::read_to_string(&out_membership).unwrap();
    assert!(rescued_membership.contains("#trackclustertu_membership_schema=v2\n"));
    assert!(rescued_membership.contains("#columns=read_id\ttu_id\tscore1\tscore2\tschema_version\tassignment_status\tbest_candidate_tu_id\t"));
    assert!(!rescued_membership.contains("#contains_legacy_rows=true"));
    assert!(!rescued_membership.contains("#row_widths="));
    let rows: Vec<&str> = rescued_membership
        .lines()
        .filter(|line| !line.starts_with('#'))
        .collect();
    assert_eq!(rows.len(), 1);
    assert_eq!(rows[0].split('\t').count(), 20);
    assert!(rows[0].starts_with("r1\tTU000001\t"));

    fs::remove_dir_all(&fixture.root).unwrap();
}

#[test]
fn no_rescue_preserves_canonical_v2_rows_and_zero_hard_count_tus() {
    let reads = concat!(
        "chr1\t100\t210\tr1\t0\t+\n",
        "chr1\t100\t210\tr2\t0\t+\n",
        "chr1\t100\t210\tr3\t0\t+\n",
    );
    let tus = concat!(
        "chr1\t100\t210\tTU000001\t0\t+\n",
        "chr1\t101\t210\tTU000002\t0\t+\n",
        "chr1\t300\t400\tTU000003\t0\t+\n",
    );
    let mut membership = concat!(
        "#trackclustertu_membership_schema=v2\n",
        "#ambiguity_margin=0.03\n",
        "#fractional_assignment=true\n",
    )
    .to_owned();
    let r1_row = "r1\tTU000001\t0.910001\t0.920002\tv2\tunique\tTU000001\t.\t.\t.\t4\t5\t.\t.\t0.920002\t.\t.\t1\tTU000001:1\ttrue\n";
    let r2_row = "r2\t.\t0.770003\t0.880004\tv2\tambiguous\tTU000001\tTU000002\t0.660005\t0.550006\t7\t8\t9\t10\t0.880004\t0.879004\t0.001\t0.5\tTU000001:0.5,TU000002:0.25,TU000003:0.25\tfalse\n";
    membership.push_str(r1_row);
    membership.push_str(r2_row);
    membership.push_str(&full_v2_row("r3", ".", "unassigned", ".", ".", "0", "."));

    let fixture = fixture("v2_no_rescue_preserves_rows", reads, tus, &membership);
    let out_tu = fixture.root.join("rescued.tus.bed");
    let out_membership = fixture.root.join("rescued.membership.tsv");
    let out_counts = fixture.root.join("rescued.counts.csv");
    let output = run_rescue(&fixture, &out_tu, &out_membership, Some(&out_counts), None);
    assert!(
        output.status.success(),
        "stdout:\n{}\nstderr:\n{}",
        String::from_utf8_lossy(&output.stdout),
        String::from_utf8_lossy(&output.stderr)
    );
    assert!(
        String::from_utf8_lossy(&output.stdout).contains("rescued_candidate_count=0"),
        "stdout:\n{}",
        String::from_utf8_lossy(&output.stdout)
    );
    assert_eq!(fs::read_to_string(&out_tu).unwrap(), tus);
    assert_eq!(
        fs::read_to_string(&out_counts).unwrap(),
        concat!(
            "tu_id,count\n",
            "TU000001,1\n",
            "TU000002,0\n",
            "TU000003,0\n",
        )
    );

    let rescued_membership = fs::read_to_string(&out_membership).unwrap();
    assert!(rescued_membership.contains("#trackclustertu_membership_schema=v2\n"));
    assert!(rescued_membership.contains("#columns=read_id\ttu_id\tscore1\tscore2\t"));
    assert!(rescued_membership.contains("#ambiguity_margin=0.03\n"));
    assert!(rescued_membership.contains("#fractional_assignment=true\n"));
    assert!(rescued_membership.contains("#rescue_unmodified_rows=preserved_from_input\n"));

    let rows_by_id: HashMap<&str, &str> = rescued_membership
        .lines()
        .filter(|line| !line.starts_with('#'))
        .map(|line| (line.split('\t').next().unwrap(), line))
        .collect();
    assert_eq!(rows_by_id["r1"], r1_row.trim_end());
    assert_eq!(rows_by_id["r2"], r2_row.trim_end());

    let rows: HashMap<&str, Vec<&str>> = rescued_membership
        .lines()
        .filter(|line| !line.starts_with('#'))
        .map(|line| {
            let fields: Vec<&str> = line.split('\t').collect();
            assert_eq!(fields.len(), 20, "row: {line}");
            (fields[0], fields)
        })
        .collect();
    assert_eq!(rows.len(), 3);
    assert_eq!(rows["r1"][1], "TU000001");
    assert_eq!(rows["r1"][5], "unique");
    assert_eq!(rows["r1"][18], "TU000001:1");
    assert_eq!(rows["r1"][19], "true");
    assert_eq!(rows["r2"][1], ".");
    assert_eq!(rows["r2"][5], "ambiguous");
    assert_eq!(rows["r2"][6], "TU000001");
    assert_eq!(rows["r2"][7], "TU000002");
    assert_eq!(rows["r2"][18], "TU000001:0.5,TU000002:0.25,TU000003:0.25");
    assert_eq!(rows["r3"][1], ".");
    assert_eq!(rows["r3"][5], "unassigned");

    fs::remove_dir_all(&fixture.root).unwrap();
}

#[test]
fn rescued_v2_rows_are_canonical_and_preserve_full_length_evidence() {
    let mut membership = concat!(
        "#trackclustertu_membership_schema=v2\n",
        "#ambiguity_margin=0.02\n",
        "#fractional_assignment=true\n",
    )
    .to_owned();
    membership.push_str(&full_v2_row_with_evidence(
        "r1",
        ".",
        "ambiguous",
        "TU000001",
        "TU000002",
        "0.5",
        "TU000001:0.5,TU000002:0.5",
        true,
    ));
    membership.push_str(&full_v2_row("r2", ".", "unassigned", ".", ".", "0", "."));
    membership.push_str(&full_v2_row(
        "r3",
        ".",
        "ambiguous",
        "TU000001",
        "TU000002",
        "0",
        ".",
    ));
    membership.push_str(&full_v2_row_with_evidence(
        "r4",
        "TU000001",
        "unique",
        "TU000001",
        ".",
        "1",
        "TU000001:1",
        true,
    ));
    membership.push_str(&full_v2_row(
        "r5",
        "TU000001",
        "partial",
        "TU000001",
        ".",
        "1",
        "TU000001:1",
    ));
    let fixture = fixture(
        "canonical_v2_rescue",
        standard_reads(),
        concat!(
            "chr1\t150\t210\tTU000001\t0\t+\n",
            "chr1\t300\t400\tTU000002\t0\t+\n",
        ),
        &membership,
    );
    let out_tu = fixture.root.join("rescued.tus.bed");
    let out_membership = fixture.root.join("rescued.membership.tsv");
    let out_counts = fixture.root.join("rescued.counts.csv");
    let output = run_rescue(&fixture, &out_tu, &out_membership, Some(&out_counts), None);
    assert!(
        output.status.success(),
        "stdout:\n{}\nstderr:\n{}",
        String::from_utf8_lossy(&output.stdout),
        String::from_utf8_lossy(&output.stderr)
    );
    assert_eq!(
        fs::read_to_string(&out_tu).unwrap(),
        concat!(
            "chr1\t100\t210\tRESC0001\t0\t+\n",
            "chr1\t150\t210\tTU000001\t0\t+\n",
            "chr1\t300\t400\tTU000002\t0\t+\n",
        )
    );
    assert_eq!(
        fs::read_to_string(&out_counts).unwrap(),
        concat!(
            "tu_id,count\n",
            "RESC0001,3\n",
            "TU000001,2\n",
            "TU000002,0\n",
        )
    );

    let rescued_membership = fs::read_to_string(&out_membership).unwrap();
    let rows: HashMap<&str, Vec<&str>> = rescued_membership
        .lines()
        .filter(|line| !line.starts_with('#'))
        .map(|line| {
            let fields: Vec<&str> = line.split('\t').collect();
            assert_eq!(fields.len(), 20, "row: {line}");
            (fields[0], fields)
        })
        .collect();
    assert_eq!(rows.len(), 5);
    for read_id in ["r1", "r2", "r3"] {
        assert_eq!(rows[read_id][1], "RESC0001");
        assert_eq!(rows[read_id][5], "unique");
        assert_eq!(rows[read_id][6], "RESC0001");
        assert_eq!(rows[read_id][7], ".");
        assert_eq!(rows[read_id][17], "1");
        assert_eq!(rows[read_id][18], "RESC0001:1");
    }
    assert_eq!(rows["r1"][19], "true");
    assert_eq!(rows["r2"][19], "false");
    assert_eq!(rows["r3"][19], "false");
    assert_eq!(rows["r4"][1], "TU000001");
    assert_eq!(rows["r4"][5], "unique");
    assert_eq!(rows["r4"][19], "true");
    assert_eq!(rows["r5"][1], "TU000001");
    assert_eq!(rows["r5"][5], "partial");

    fs::remove_dir_all(&fixture.root).unwrap();
}

#[test]
fn duplicate_read_id_copies_are_quarantined_and_recorded() {
    let fixture = fixture(
        "duplicate_read",
        concat!("chr1\t100\t210\tr1\t0\t+\n", "chr1\t102\t210\tr1\t0\t+\n",),
        "chr1\t150\t210\tTU000001\t0\t+\n",
        "r1\tTU000001\n",
    );
    let output_dir = fixture.root.join("must-not-be-created");
    let out_tu = output_dir.join("rescued.tus.bed");
    let out_membership = output_dir.join("rescued.membership.tsv");

    let output = run_rescue(&fixture, &out_tu, &out_membership, None, None);
    assert!(
        output.status.success(),
        "stdout:\n{}\nstderr:\n{}",
        String::from_utf8_lossy(&output.stdout),
        String::from_utf8_lossy(&output.stderr)
    );
    assert_eq!(
        fs::read_to_string(&out_tu).unwrap(),
        "chr1\t150\t210\tTU000001\t0\t+\n"
    );
    assert_eq!(fs::read_to_string(&out_membership).unwrap(), "");

    let rejection_path = output_dir.join("rescue_read_rejections.tsv");
    let rejection_text = fs::read_to_string(&rejection_path).unwrap();
    assert_eq!(rejection_text.matches("\tduplicate_read_id\t").count(), 2);
    assert!(rejection_text.contains("\t1\tr1\tdeduplicate\t"));
    assert!(rejection_text.contains("\t2\tr1\tdeduplicate\t"));
    let stderr = String::from_utf8_lossy(&output.stderr);
    assert!(
        stderr.contains("total=2\tretained=0\trejected=2\tduplicate_read_id=2"),
        "stderr:\n{stderr}"
    );

    fs::remove_dir_all(&fixture.root).unwrap();
}

#[test]
fn malformed_bed_read_is_skipped_without_turning_its_membership_into_a_global_error() {
    let fixture = fixture(
        "malformed_bed_read_membership",
        concat!(
            "chr1\t100\t210\tr1\t0\t+\n",
            "chr1\tbad\t210\tbad_read\t0\t+\n",
            "chr1\t102\t210\tr2\t0\t+\n",
            "chr1\t104\t210\tr3\t0\t+\n",
        ),
        "chr1\t150\t210\tTU000001\t0\t+\n",
        concat!(
            "r1\tTU000001\n",
            "bad_read\tTU000001\n",
            "r2\tTU000001\n",
            "r3\tTU000001\n",
        ),
    );
    let out_tu = fixture.root.join("output/rescued.tus.bed");
    let out_membership = fixture.root.join("output/rescued.membership.tsv");

    let output = run_rescue(&fixture, &out_tu, &out_membership, None, None);
    assert!(
        output.status.success(),
        "stdout:\n{}\nstderr:\n{}",
        String::from_utf8_lossy(&output.stdout),
        String::from_utf8_lossy(&output.stderr)
    );
    let rescued_membership = fs::read_to_string(&out_membership).unwrap();
    assert!(!rescued_membership.contains("bad_read"));
    assert_eq!(rescued_membership.lines().count(), 3);
    let rejection_text =
        fs::read_to_string(fixture.root.join("output/rescue_read_rejections.tsv")).unwrap();
    assert!(rejection_text.contains("\t2\t.\tparse\tinvalid_integer\t"));
    assert!(String::from_utf8_lossy(&output.stderr)
        .contains("total=4\tretained=3\trejected=1\tinvalid_integer=1"));

    fs::remove_dir_all(&fixture.root).unwrap();
}

#[test]
fn invalid_same_id_rows_do_not_remove_valid_rescue_read_or_break_tu_membership() {
    let fixture = fixture(
        "invalid_same_id_rescue",
        concat!(
            "chr1\t100\t210\tr1\t0\t+\n",
            "chr1\t100\t210\tr1\t0\t.\n",
            "chr1\t220\t210\tr1\t0\t+\n",
            "chr1\t102\t210\tr2\t0\t+\n",
            "chr1\t104\t210\tr3\t0\t+\n",
        ),
        "chr1\t150\t210\tTU000001\t0\t+\n",
        concat!("r1\tTU000001\n", "r2\tTU000001\n", "r3\tTU000001\n",),
    );
    let out_tu = fixture.root.join("output/rescued.tus.bed");
    let out_membership = fixture.root.join("output/rescued.membership.tsv");

    let output = run_rescue(&fixture, &out_tu, &out_membership, None, None);
    assert!(
        output.status.success(),
        "stdout:\n{}\nstderr:\n{}",
        String::from_utf8_lossy(&output.stdout),
        String::from_utf8_lossy(&output.stderr)
    );
    let rescued_membership = fs::read_to_string(&out_membership).unwrap();
    assert_eq!(rescued_membership.lines().count(), 3);
    assert_eq!(
        rescued_membership
            .lines()
            .filter(|line| line.starts_with("r1\t"))
            .count(),
        1
    );
    let rescued_tus = fs::read_to_string(&out_tu).unwrap();
    let emitted_tu_ids: HashSet<&str> = rescued_tus
        .lines()
        .map(|line| line.split('\t').nth(3).unwrap())
        .collect();
    for line in rescued_membership.lines() {
        let tu_id = line.split('\t').nth(1).unwrap();
        assert!(emitted_tu_ids.contains(tu_id), "missing rescued TU {tu_id}");
    }

    let rejection_text =
        fs::read_to_string(fixture.root.join("output/rescue_read_rejections.tsv")).unwrap();
    assert!(rejection_text.contains("\t2\tr1\tvalidate\tunknown_strand\t"));
    assert!(rejection_text.contains("\t3\t.\tparse\tinvalid_interval\t"));
    assert!(!rejection_text.contains("\tduplicate_read_id\t"));
    assert!(String::from_utf8_lossy(&output.stderr).contains("total=5\tretained=3\trejected=2"));

    fs::remove_dir_all(&fixture.root).unwrap();
}

#[test]
fn strict_read_errors_rejects_duplicate_records_atomically() {
    let fixture = fixture(
        "strict_duplicate_read",
        concat!("chr1\t100\t210\tr1\t0\t+\n", "chr1\t102\t210\tr1\t0\t+\n",),
        "chr1\t150\t210\tTU000001\t0\t+\n",
        "r1\tTU000001\n",
    );
    let output_dir = fixture.root.join("must-not-be-created");
    let out_tu = output_dir.join("rescued.tus.bed");
    let out_membership = output_dir.join("rescued.membership.tsv");

    let exe = env!("CARGO_BIN_EXE_trackclustertu");
    let output = Command::new(exe)
        .arg("rescue-missed-tus")
        .arg("--in")
        .arg(&fixture.reads)
        .arg("--existing-tu")
        .arg(&fixture.tus)
        .arg("--existing-membership")
        .arg(&fixture.membership)
        .arg("--out-tu")
        .arg(&out_tu)
        .arg("--out-membership")
        .arg(&out_membership)
        .arg("--strict-read-errors")
        .output()
        .unwrap();

    assert!(!output.status.success());
    let stderr = String::from_utf8_lossy(&output.stderr);
    assert!(stderr.contains("--strict-read-errors rejected 2 read record(s)"));
    assert!(stderr.contains("duplicate read ID"));
    assert!(!output_dir.exists());

    fs::remove_dir_all(&fixture.root).unwrap();
}

#[test]
fn duplicate_existing_tu_id_fails_before_outputs_are_created() {
    let fixture = fixture(
        "duplicate_tu",
        standard_reads(),
        concat!(
            "chr1\t150\t210\tTU000001\t0\t+\n",
            "chr1\t300\t400\tTU000001\t0\t+\n",
        ),
        &standard_membership("TU000001"),
    );
    let out_tu = fixture.root.join("rescued.tus.bed");
    let out_membership = fixture.root.join("rescued.membership.tsv");

    let output = run_rescue(&fixture, &out_tu, &out_membership, None, None);
    assert_failed_with_line(&output, &fixture.tus, 2, "duplicate existing TU ID");
    assert!(!out_tu.exists());
    assert!(!out_membership.exists());

    fs::remove_dir_all(&fixture.root).unwrap();
}

#[test]
fn duplicate_membership_read_id_fails_instead_of_overwriting_a_row() {
    let fixture = fixture(
        "duplicate_membership",
        standard_reads(),
        "chr1\t150\t210\tTU000001\t0\t+\n",
        concat!("r1\tTU000001\t1.0\t1.0\n", "r1\tTU000001\t1.0\t1.0\n",),
    );
    let out_tu = fixture.root.join("rescued.tus.bed");
    let out_membership = fixture.root.join("rescued.membership.tsv");

    let output = run_rescue(&fixture, &out_tu, &out_membership, None, None);
    assert_failed_with_line(
        &output,
        &fixture.membership,
        2,
        "duplicate membership for read ID",
    );
    assert!(!out_tu.exists());
    assert!(!out_membership.exists());

    fs::remove_dir_all(&fixture.root).unwrap();
}

#[test]
fn membership_foreign_keys_fail_with_line_context() {
    let fixture = fixture(
        "membership_foreign_keys",
        standard_reads(),
        "chr1\t150\t210\tTU000001\t0\t+\n",
        "missing_read\tTU000001\n",
    );
    let out_tu = fixture.root.join("rescued.tus.bed");
    let out_membership = fixture.root.join("rescued.membership.tsv");

    let unknown_read = run_rescue(&fixture, &out_tu, &out_membership, None, None);
    assert_failed_with_line(
        &unknown_read,
        &fixture.membership,
        1,
        "does not exist in the input BED",
    );
    assert!(!out_tu.exists());
    assert!(!out_membership.exists());

    fs::write(&fixture.membership, "r1\tmissing_tu\n").unwrap();
    let unknown_tu = run_rescue(&fixture, &out_tu, &out_membership, None, None);
    assert_failed_with_line(
        &unknown_tu,
        &fixture.membership,
        1,
        "does not exist in the existing TU BED",
    );
    assert!(!out_tu.exists());
    assert!(!out_membership.exists());

    fs::remove_dir_all(&fixture.root).unwrap();
}

#[test]
fn all_v2_membership_tu_references_are_foreign_keys() {
    let fixture = fixture(
        "all_membership_foreign_keys",
        standard_reads(),
        concat!(
            "chr1\t150\t210\tTU000001\t0\t+\n",
            "chr1\t300\t400\tTU000002\t0\t+\n",
        ),
        "placeholder\tTU000001\n",
    );
    let cases = [
        (
            "best",
            full_v2_row(
                "r1",
                ".",
                "ambiguous",
                "MISSING_BEST",
                "TU000001",
                "0.5",
                "MISSING_BEST:0.5,TU000001:0.5",
            ),
            "MISSING_BEST",
        ),
        (
            "second",
            full_v2_row(
                "r1",
                ".",
                "ambiguous",
                "TU000001",
                "MISSING_SECOND",
                "0.5",
                "TU000001:0.5,MISSING_SECOND:0.5",
            ),
            "MISSING_SECOND",
        ),
        (
            "third-fraction",
            full_v2_row(
                "r1",
                ".",
                "ambiguous",
                "TU000001",
                "TU000002",
                "0.5",
                "TU000001:0.5,TU000002:0.25,MISSING_THIRD:0.25",
            ),
            "MISSING_THIRD",
        ),
    ];

    for (label, membership, missing_tu) in cases {
        fs::write(&fixture.membership, membership).unwrap();
        let out_tu = fixture.root.join(format!("{label}.tus.bed"));
        let out_membership = fixture.root.join(format!("{label}.membership.tsv"));
        let output = run_rescue(&fixture, &out_tu, &out_membership, None, None);
        assert_failed_with_line(
            &output,
            &fixture.membership,
            1,
            "does not exist in the existing TU BED",
        );
        let stderr = String::from_utf8_lossy(&output.stderr);
        assert!(stderr.contains(missing_tu), "stderr:\n{stderr}");
        assert!(stderr.contains("read \"r1\""), "stderr:\n{stderr}");
        assert!(!out_tu.exists());
        assert!(!out_membership.exists());
    }

    fs::remove_dir_all(&fixture.root).unwrap();
}

#[test]
fn membership_foreign_keys_preserve_surrounding_identifier_spaces() {
    let fixture = fixture(
        "membership_spaced_foreign_keys",
        concat!(
            "chr1\t100\t210\t r1 \t0\t+\n",
            "chr1\t102\t210\t r2 \t0\t+\n",
            "chr1\t104\t210\t r3 \t0\t+\n",
            "chr1\t150\t210\t r4 \t0\t+\n",
            "chr1\t151\t210\t r5 \t0\t+\n",
        ),
        "chr1\t150\t210\t TU000001 \t0\t+\n",
        concat!(
            " r1 \t TU000001 \t1.0\t1.0\n",
            " r2 \t TU000001 \t1.0\t1.0\n",
            " r3 \t TU000001 \t1.0\t1.0\n",
            " r4 \t TU000001 \t1.0\t1.0\n",
            " r5 \t TU000001 \t1.0\t1.0\n",
        ),
    );
    let out_tu = fixture.root.join("rescued.tus.bed");
    let out_membership = fixture.root.join("rescued.membership.tsv");

    let output = run_rescue(&fixture, &out_tu, &out_membership, None, None);
    assert!(
        output.status.success(),
        "stdout:\n{}\nstderr:\n{}",
        String::from_utf8_lossy(&output.stdout),
        String::from_utf8_lossy(&output.stderr)
    );
    let membership = fs::read_to_string(&out_membership).unwrap();
    assert!(membership.contains(" r4 \t TU000001 \t"), "{membership}");
    assert!(membership.contains(" r5 \t TU000001 \t"), "{membership}");

    fs::remove_dir_all(&fixture.root).unwrap();
}

#[test]
fn invalid_rescue_prefixes_fail_before_outputs_are_created() {
    let fixture = fixture(
        "invalid_rescue_prefix",
        standard_reads(),
        "chr1\t150\t210\tTU000001\t0\t+\n",
        &standard_membership("TU000001"),
    );

    for (idx, prefix) in ["", "   ", "BAD\tID", "BAD\nID", "BAD\rID", "BAD,ID"]
        .into_iter()
        .enumerate()
    {
        let output_dir = fixture.root.join(format!("must-not-be-created-{idx}"));
        let out_tu = output_dir.join("rescued.tus.bed");
        let out_membership = output_dir.join("rescued.membership.tsv");
        let output = run_rescue(&fixture, &out_tu, &out_membership, None, Some(prefix));
        assert!(!output.status.success(), "prefix {prefix:?} was accepted");
        let stderr = String::from_utf8_lossy(&output.stderr);
        assert!(stderr.contains("--rescue-prefix"), "stderr:\n{stderr}");
        assert!(!out_tu.exists());
        assert!(!out_membership.exists());
        assert!(!output_dir.exists());
    }

    fs::remove_dir_all(&fixture.root).unwrap();
}
