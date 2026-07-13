use std::fs;
use std::io::Write;
use std::process::Command;
use std::time::{SystemTime, UNIX_EPOCH};

mod support;
use support::{
    count_v1_projection, membership_v1_projection, metric_long_total_projection,
    metric_matrix_total_projection,
};

fn unique_tmp_dir(prefix: &str) -> std::path::PathBuf {
    let nanos = SystemTime::now()
        .duration_since(UNIX_EPOCH)
        .expect("clock")
        .as_nanos();
    std::env::temp_dir().join(format!("{prefix}_{nanos}"))
}

fn write_sample_bed(path: &std::path::Path, lines: &[&str]) {
    let mut file = fs::File::create(path).unwrap();
    for line in lines {
        writeln!(file, "{line}").unwrap();
    }
}

#[allow(clippy::too_many_arguments)]
fn canonical_v2_membership_row(
    read_id: &str,
    hard_tu_id: &str,
    status: &str,
    best_candidate: &str,
    second_candidate: &str,
    primary_weight: &str,
    fractional_assignments: &str,
    full_length_evidence: bool,
) -> String {
    format!(
        "{read_id}\t{hard_tu_id}\t.\t.\tv2\t{status}\t{best_candidate}\t{second_candidate}\t.\t.\t.\t.\t.\t.\t.\t.\t.\t{primary_weight}\t{fractional_assignments}\t{full_length_evidence}\n"
    )
}

#[test]
fn trackclustertu_clusters_manifest_and_writes_sample_tables() {
    let tmp = unique_tmp_dir("trackclustertu_multisample_test");
    let reads_dir = tmp.join("reads");
    fs::create_dir_all(&reads_dir).unwrap();

    write_sample_bed(
        &reads_dir.join("sampleA.bed"),
        &[
            "chr1\t100\t200\tr1\t0\t+",
            "chr1\t101\t201\tr2\t0\t+",
            "chr1\t300\t400\tr3\t0\t+",
        ],
    );
    write_sample_bed(
        &reads_dir.join("sampleB.bed"),
        &[
            "chr1\t120\t180\tr1\t0\t+",
            "chr1\t300\t400\tr2\t0\t+",
            "chr1\t500\t600\tr3\t0\t+",
        ],
    );

    let manifest = tmp.join("samples.tsv");
    fs::write(
        &manifest,
        concat!(
            "sample\treads\tgroup\n",
            "sampleA\treads/sampleA.bed\tcontrol\n",
            "sampleB\treads/sampleB.bed\ttreated\n",
        ),
    )
    .unwrap();

    let out_tu = tmp.join("tus.bed");
    let out_membership = tmp.join("membership.tsv");
    let out_pooled_reads = tmp.join("pooled.bed");
    let out_tu_count = tmp.join("tu_count.csv");
    let out_sample_long = tmp.join("sample_long.tsv");
    let out_sample_matrix = tmp.join("sample_matrix.tsv");
    let out_group_matrix = tmp.join("group_matrix.tsv");
    let genes_path = tmp.join("genes.bed");
    let out_gene_count = tmp.join("gene_count.csv");
    let out_gene_sample_matrix = tmp.join("gene_sample_matrix.tsv");
    let out_gene_group_matrix = tmp.join("gene_group_matrix.tsv");

    write_sample_bed(
        &genes_path,
        &[
            "chr1\t90\t150\tgeneA\t0\t+",
            "chr1\t160\t220\tgeneB\t0\t+",
            "chr1\t290\t310\tgeneC\t0\t+",
            "chr1\t490\t520\tgeneD\t0\t+",
        ],
    );

    let exe = env!("CARGO_BIN_EXE_trackclustertu");
    let output = Command::new(exe)
        .args([
            "cluster",
            "--tu-id-style",
            "sequential",
            "--manifest",
            manifest.to_str().unwrap(),
            "--format",
            "bed6",
            "--out-tu",
            out_tu.to_str().unwrap(),
            "--out-membership",
            out_membership.to_str().unwrap(),
            "--out-pooled-reads",
            out_pooled_reads.to_str().unwrap(),
            "--out-tu-count",
            out_tu_count.to_str().unwrap(),
            "--out-tu-sample-count-long",
            out_sample_long.to_str().unwrap(),
            "--out-tu-sample-count-matrix",
            out_sample_matrix.to_str().unwrap(),
            "--out-tu-group-count-matrix",
            out_group_matrix.to_str().unwrap(),
            "--annotation-bed",
            genes_path.to_str().unwrap(),
            "--out-gene-count",
            out_gene_count.to_str().unwrap(),
            "--out-gene-sample-count-matrix",
            out_gene_sample_matrix.to_str().unwrap(),
            "--out-gene-group-count-matrix",
            out_gene_group_matrix.to_str().unwrap(),
        ])
        .output()
        .unwrap();

    assert!(
        output.status.success(),
        "stdout:\n{}\nstderr:\n{}",
        String::from_utf8_lossy(&output.stdout),
        String::from_utf8_lossy(&output.stderr)
    );

    let sample_matrix_header = fs::read_to_string(&out_sample_matrix)
        .unwrap()
        .lines()
        .next()
        .unwrap()
        .to_owned();
    assert_eq!(
        sample_matrix_header,
        concat!(
            "tu_id",
            "\tsampleA.unique_count\tsampleA.full_length_evidence_count\tsampleA.total_count\tsampleA.fractional_count",
            "\tsampleB.unique_count\tsampleB.full_length_evidence_count\tsampleB.total_count\tsampleB.fractional_count",
        )
    );
    let group_matrix_header = fs::read_to_string(&out_group_matrix)
        .unwrap()
        .lines()
        .next()
        .unwrap()
        .to_owned();
    assert_eq!(
        group_matrix_header,
        concat!(
            "tu_id",
            "\tcontrol.unique_count\tcontrol.full_length_evidence_count\tcontrol.total_count\tcontrol.fractional_count",
            "\ttreated.unique_count\ttreated.full_length_evidence_count\ttreated.total_count\ttreated.fractional_count",
        )
    );
    for path in [
        &out_gene_count,
        &out_gene_sample_matrix,
        &out_gene_group_matrix,
    ] {
        assert!(
            fs::read_to_string(path)
                .unwrap()
                .starts_with("#count_semantics=nonexclusive_same_strand_tu_overlap;"),
            "missing nonexclusive multi-gene count metadata in {}",
            path.display()
        );
    }

    assert_eq!(
        fs::read_to_string(&out_tu).unwrap(),
        concat!(
            "chr1\t100\t200\tTU000001\t0\t+\n",
            "chr1\t120\t180\tTU000002\t0\t+\n",
            "chr1\t300\t400\tTU000003\t0\t+\n",
            "chr1\t500\t600\tTU000004\t0\t+\n",
        )
    );

    assert_eq!(
        membership_v1_projection(&fs::read_to_string(&out_membership).unwrap()),
        concat!(
            "sampleA::r1\tTU000001\t1.000000\t1.000000\n",
            "sampleA::r2\tTU000001\t0.980198\t0.990000\n",
            "sampleB::r1\tTU000002\t1.000000\t1.000000\n",
            "sampleA::r3\tTU000003\t1.000000\t1.000000\n",
            "sampleB::r2\tTU000003\t1.000000\t1.000000\n",
            "sampleB::r3\tTU000004\t1.000000\t1.000000\n",
        )
    );

    assert_eq!(
        fs::read_to_string(&out_pooled_reads).unwrap(),
        concat!(
            "chr1\t100\t200\tsampleA::r1\t0\t+\n",
            "chr1\t101\t201\tsampleA::r2\t0\t+\n",
            "chr1\t300\t400\tsampleA::r3\t0\t+\n",
            "chr1\t120\t180\tsampleB::r1\t0\t+\n",
            "chr1\t300\t400\tsampleB::r2\t0\t+\n",
            "chr1\t500\t600\tsampleB::r3\t0\t+\n",
        )
    );

    assert_eq!(
        count_v1_projection(&fs::read_to_string(&out_tu_count).unwrap()),
        concat!(
            "tu_id,count\n",
            "TU000001,2\n",
            "TU000002,1\n",
            "TU000003,2\n",
            "TU000004,1\n",
        )
    );

    assert_eq!(
        metric_long_total_projection(&fs::read_to_string(&out_sample_long).unwrap()),
        concat!(
            "tu_id\tsample\tcount\n",
            "TU000001\tsampleA\t2\n",
            "TU000002\tsampleB\t1\n",
            "TU000003\tsampleA\t1\n",
            "TU000003\tsampleB\t1\n",
            "TU000004\tsampleB\t1\n",
        )
    );

    assert_eq!(
        metric_matrix_total_projection(&fs::read_to_string(&out_sample_matrix).unwrap()),
        concat!(
            "tu_id\tsampleA\tsampleB\n",
            "TU000001\t2\t0\n",
            "TU000002\t0\t1\n",
            "TU000003\t1\t1\n",
            "TU000004\t0\t1\n",
        )
    );

    assert_eq!(
        metric_matrix_total_projection(&fs::read_to_string(&out_group_matrix).unwrap()),
        concat!(
            "tu_id\tcontrol\ttreated\n",
            "TU000001\t2\t0\n",
            "TU000002\t0\t1\n",
            "TU000003\t1\t1\n",
            "TU000004\t0\t1\n",
        )
    );

    assert_eq!(
        count_v1_projection(&fs::read_to_string(&out_gene_count).unwrap()),
        concat!(
            "gene_id,count\n",
            "geneA,3\n",
            "geneB,3\n",
            "geneC,2\n",
            "geneD,1\n",
        )
    );

    assert_eq!(
        metric_matrix_total_projection(&fs::read_to_string(&out_gene_sample_matrix).unwrap()),
        concat!(
            "gene_id\tsampleA\tsampleB\n",
            "geneA\t2\t1\n",
            "geneB\t2\t1\n",
            "geneC\t1\t1\n",
            "geneD\t0\t1\n",
        )
    );

    assert_eq!(
        metric_matrix_total_projection(&fs::read_to_string(&out_gene_group_matrix).unwrap()),
        concat!(
            "gene_id\tcontrol\ttreated\n",
            "geneA\t2\t1\n",
            "geneB\t2\t1\n",
            "geneC\t1\t1\n",
            "geneD\t0\t1\n",
        )
    );

    let _ = fs::remove_dir_all(&tmp);
}

#[test]
fn trackclustertu_recounts_multi_sample_tables_from_pooled_membership() {
    let tmp = unique_tmp_dir("trackclustertu_multisample_recount_test");
    let reads_dir = tmp.join("reads");
    fs::create_dir_all(&reads_dir).unwrap();

    write_sample_bed(
        &reads_dir.join("sampleA.bed"),
        &[
            "chr1\t100\t200\tr1\t0\t+",
            "chr1\t101\t201\tr2\t0\t+",
            "chr1\t300\t400\tr3\t0\t+",
        ],
    );
    write_sample_bed(
        &reads_dir.join("sampleB.bed"),
        &[
            "chr1\t120\t180\tr1\t0\t+",
            "chr1\t300\t400\tr2\t0\t+",
            "chr1\t500\t600\tr3\t0\t+",
        ],
    );

    let manifest = tmp.join("samples.tsv");
    fs::write(
        &manifest,
        concat!(
            "sample\treads\tgroup\n",
            "sampleA\treads/sampleA.bed\tcontrol\n",
            "sampleB\treads/sampleB.bed\ttreated\n",
        ),
    )
    .unwrap();

    let membership = tmp.join("membership.tsv");
    let exe = env!("CARGO_BIN_EXE_trackclustertu");
    let initial = Command::new(exe)
        .args([
            "cluster",
            "--tu-id-style",
            "sequential",
            "--manifest",
            manifest.to_str().unwrap(),
            "--format",
            "bed6",
            "--out-tu",
            tmp.join("tus.bed").to_str().unwrap(),
            "--out-membership",
            membership.to_str().unwrap(),
        ])
        .output()
        .unwrap();
    assert!(
        initial.status.success(),
        "stdout:\n{}\nstderr:\n{}",
        String::from_utf8_lossy(&initial.stdout),
        String::from_utf8_lossy(&initial.stderr)
    );

    let out_tu_count = tmp.join("recount_tu_count.csv");
    let out_sample_matrix = tmp.join("recount_sample_matrix.tsv");
    let out_group_matrix = tmp.join("recount_group_matrix.tsv");
    let output = Command::new(exe)
        .args([
            "recount",
            "--manifest",
            manifest.to_str().unwrap(),
            "--pooled-membership",
            membership.to_str().unwrap(),
            "--out-tu-count",
            out_tu_count.to_str().unwrap(),
            "--out-tu-sample-count-matrix",
            out_sample_matrix.to_str().unwrap(),
            "--out-tu-group-count-matrix",
            out_group_matrix.to_str().unwrap(),
        ])
        .output()
        .unwrap();

    assert!(
        output.status.success(),
        "stdout:\n{}\nstderr:\n{}",
        String::from_utf8_lossy(&output.stdout),
        String::from_utf8_lossy(&output.stderr)
    );

    assert_eq!(
        count_v1_projection(&fs::read_to_string(&out_tu_count).unwrap()),
        concat!(
            "tu_id,count\n",
            "TU000001,2\n",
            "TU000002,1\n",
            "TU000003,2\n",
            "TU000004,1\n",
        )
    );
    assert_eq!(
        metric_matrix_total_projection(&fs::read_to_string(&out_sample_matrix).unwrap()),
        concat!(
            "tu_id\tsampleA\tsampleB\n",
            "TU000001\t2\t0\n",
            "TU000002\t0\t1\n",
            "TU000003\t1\t1\n",
            "TU000004\t0\t1\n",
        )
    );
    assert_eq!(
        metric_matrix_total_projection(&fs::read_to_string(&out_group_matrix).unwrap()),
        concat!(
            "tu_id\tcontrol\ttreated\n",
            "TU000001\t2\t0\n",
            "TU000002\t0\t1\n",
            "TU000003\t1\t1\n",
            "TU000004\t0\t1\n",
        )
    );

    let _ = fs::remove_dir_all(&tmp);
}

#[test]
fn rescued_v2_membership_recounts_promotions_and_untouched_fractions() {
    let tmp = unique_tmp_dir("trackclustertu_rescue_recount_test");
    let reads_dir = tmp.join("reads");
    fs::create_dir_all(&reads_dir).unwrap();

    write_sample_bed(
        &reads_dir.join("sampleA.bed"),
        &[
            "chr1\t100\t210\tr1\t0\t+",
            "chr1\t102\t210\tr2\t0\t+",
            "chr1\t150\t210\tr4\t0\t+",
            "chr1\t500\t600\tamb\t0\t+",
        ],
    );
    write_sample_bed(
        &reads_dir.join("sampleB.bed"),
        &["chr1\t104\t210\tr3\t0\t+", "chr1\t151\t210\tr5\t0\t+"],
    );

    let manifest = tmp.join("samples.tsv");
    fs::write(
        &manifest,
        concat!(
            "sample\treads\tgroup\n",
            "sampleA\treads/sampleA.bed\tcontrol\n",
            "sampleB\treads/sampleB.bed\ttreated\n",
        ),
    )
    .unwrap();

    let pooled_reads = tmp.join("pooled.bed");
    write_sample_bed(
        &pooled_reads,
        &[
            "chr1\t100\t210\tsampleA::r1\t0\t+",
            "chr1\t102\t210\tsampleA::r2\t0\t+",
            "chr1\t104\t210\tsampleB::r3\t0\t+",
            "chr1\t150\t210\tsampleA::r4\t0\t+",
            "chr1\t151\t210\tsampleB::r5\t0\t+",
            "chr1\t500\t600\tsampleA::amb\t0\t+",
        ],
    );

    let existing_tus = tmp.join("tus.bed");
    write_sample_bed(
        &existing_tus,
        &[
            "chr1\t150\t210\tTU_MAIN\t0\t+",
            "chr1\t500\t600\tTU_A\t0\t+",
            "chr1\t501\t600\tTU_B\t0\t+",
            "chr1\t800\t900\tTU_ZERO\t0\t+",
        ],
    );

    let mut membership = concat!(
        "#trackclustertu_membership_schema=v2\n",
        "#fractional_assignment=true\n",
    )
    .to_owned();
    for read_id in ["sampleA::r1", "sampleA::r2", "sampleB::r3"] {
        membership.push_str(&canonical_v2_membership_row(
            read_id,
            ".",
            "unassigned",
            ".",
            ".",
            "0",
            ".",
            false,
        ));
    }
    membership.push_str(&canonical_v2_membership_row(
        "sampleA::r4",
        "TU_MAIN",
        "unique",
        "TU_MAIN",
        ".",
        "1",
        "TU_MAIN:1",
        true,
    ));
    membership.push_str(&canonical_v2_membership_row(
        "sampleB::r5",
        "TU_MAIN",
        "unique",
        "TU_MAIN",
        ".",
        "1",
        "TU_MAIN:1",
        false,
    ));
    membership.push_str(&canonical_v2_membership_row(
        "sampleA::amb",
        ".",
        "ambiguous",
        "TU_A",
        "TU_B",
        "0.5",
        "TU_A:0.5,TU_B:0.5",
        false,
    ));
    let membership_path = tmp.join("membership.tsv");
    fs::write(&membership_path, membership).unwrap();

    let rescued_tus = tmp.join("rescue/rescued.tus.bed");
    let rescued_membership = tmp.join("rescue/rescued.membership.tsv");
    let rescued_counts = tmp.join("rescue/rescued.tu_count.csv");
    let exe = env!("CARGO_BIN_EXE_trackclustertu");
    let rescue = Command::new(exe)
        .args([
            "rescue-missed-tus",
            "--in",
            pooled_reads.to_str().unwrap(),
            "--existing-tu",
            existing_tus.to_str().unwrap(),
            "--existing-membership",
            membership_path.to_str().unwrap(),
            "--out-tu",
            rescued_tus.to_str().unwrap(),
            "--out-membership",
            rescued_membership.to_str().unwrap(),
            "--out-tu-count",
            rescued_counts.to_str().unwrap(),
            "--min-family-support",
            "2",
            "--min-mode-support",
            "2",
            "--min-mode-fraction",
            "0.2",
        ])
        .output()
        .unwrap();
    assert!(
        rescue.status.success(),
        "stdout:\n{}\nstderr:\n{}",
        String::from_utf8_lossy(&rescue.stdout),
        String::from_utf8_lossy(&rescue.stderr)
    );
    assert!(String::from_utf8_lossy(&rescue.stdout).contains("rescued_candidate_count=1"));

    let rescued_membership_text = fs::read_to_string(&rescued_membership).unwrap();
    for read_id in ["sampleA::r1", "sampleA::r2", "sampleB::r3"] {
        let row = rescued_membership_text
            .lines()
            .find(|line| line.starts_with(&format!("{read_id}\t")))
            .unwrap();
        let fields: Vec<&str> = row.split('\t').collect();
        assert_eq!(fields[1], "RESC0001");
        assert_eq!(fields[5], "unique");
        assert_eq!(fields[18], "RESC0001:1");
    }
    let ambiguous_row = rescued_membership_text
        .lines()
        .find(|line| line.starts_with("sampleA::amb\t"))
        .unwrap();
    assert_eq!(
        ambiguous_row.split('\t').collect::<Vec<_>>()[18],
        "TU_A:0.5,TU_B:0.5"
    );

    let rescued_count_text = fs::read_to_string(&rescued_counts).unwrap();
    assert!(rescued_count_text.contains("RESC0001,3\n"));
    assert!(rescued_count_text.contains("TU_MAIN,2\n"));
    assert!(rescued_count_text.contains("TU_ZERO,0\n"));

    let recount_count = tmp.join("recount/tu_count.csv");
    let recount_sample_matrix = tmp.join("recount/tu_sample_matrix.tsv");
    let recount = Command::new(exe)
        .args([
            "recount",
            "--manifest",
            manifest.to_str().unwrap(),
            "--pooled-membership",
            rescued_membership.to_str().unwrap(),
            "--out-tu-count",
            recount_count.to_str().unwrap(),
            "--out-tu-sample-count-matrix",
            recount_sample_matrix.to_str().unwrap(),
        ])
        .output()
        .unwrap();
    assert!(
        recount.status.success(),
        "stdout:\n{}\nstderr:\n{}",
        String::from_utf8_lossy(&recount.stdout),
        String::from_utf8_lossy(&recount.stderr)
    );

    let recount_count_text = fs::read_to_string(&recount_count).unwrap();
    assert!(recount_count_text.contains("RESC0001,3,3,0,3,3\n"));
    assert!(recount_count_text.contains("TU_MAIN,2,2,1,2,2\n"));
    assert!(recount_count_text.contains("TU_A,0,0,0,0,0.5\n"));
    assert!(recount_count_text.contains("TU_B,0,0,0,0,0.5\n"));
    assert!(!recount_count_text.contains("TU_ZERO"));

    let sample_matrix_text = fs::read_to_string(&recount_sample_matrix).unwrap();
    assert!(sample_matrix_text.contains("RESC0001\t2\t0\t2\t2\t1\t0\t1\t1\n"));
    assert!(sample_matrix_text.contains("TU_A\t0\t0\t0\t0.5\t0\t0\t0\t0\n"));

    let _ = fs::remove_dir_all(&tmp);
}

#[test]
fn trackclustertu_skips_group_matrix_when_manifest_has_no_group_column() {
    let tmp = unique_tmp_dir("trackclustertu_multisample_no_group_test");
    let reads_dir = tmp.join("reads");
    fs::create_dir_all(&reads_dir).unwrap();

    write_sample_bed(
        &reads_dir.join("sampleA.bed"),
        &["chr1\t100\t200\tr1\t0\t+"],
    );
    write_sample_bed(
        &reads_dir.join("sampleB.bed"),
        &["chr1\t100\t200\tr1\t0\t+"],
    );

    let manifest = tmp.join("samples.tsv");
    fs::write(
        &manifest,
        concat!(
            "sample\treads\n",
            "sampleA\treads/sampleA.bed\n",
            "sampleB\treads/sampleB.bed\n",
        ),
    )
    .unwrap();

    let out_group_matrix = tmp.join("group_matrix.tsv");
    let out_sample_matrix = tmp.join("sample_matrix.tsv");

    let exe = env!("CARGO_BIN_EXE_trackclustertu");
    let output = Command::new(exe)
        .args([
            "cluster",
            "--tu-id-style",
            "sequential",
            "--manifest",
            manifest.to_str().unwrap(),
            "--format",
            "bed6",
            "--out-tu",
            tmp.join("tus.bed").to_str().unwrap(),
            "--out-membership",
            tmp.join("membership.tsv").to_str().unwrap(),
            "--out-tu-sample-count-matrix",
            out_sample_matrix.to_str().unwrap(),
            "--out-tu-group-count-matrix",
            out_group_matrix.to_str().unwrap(),
        ])
        .output()
        .unwrap();

    assert!(
        output.status.success(),
        "stdout:\n{}\nstderr:\n{}",
        String::from_utf8_lossy(&output.stdout),
        String::from_utf8_lossy(&output.stderr)
    );

    assert_eq!(
        metric_matrix_total_projection(&fs::read_to_string(&out_sample_matrix).unwrap()),
        concat!("tu_id\tsampleA\tsampleB\n", "TU000001\t1\t1\n",)
    );
    assert!(!out_group_matrix.exists());

    let _ = fs::remove_dir_all(&tmp);
}
