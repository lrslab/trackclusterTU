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

fn write_sample_bed(path: &std::path::Path, lines: &[&str]) {
    let mut file = fs::File::create(path).unwrap();
    for line in lines {
        writeln!(file, "{line}").unwrap();
    }
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

    assert_eq!(
        fs::read_to_string(&out_tu).unwrap(),
        concat!(
            "chr1\t100\t200\tTU000001\t0\t+\n",
            "chr1\t300\t400\tTU000002\t0\t+\n",
            "chr1\t500\t600\tTU000003\t0\t+\n",
        )
    );

    assert_eq!(
        fs::read_to_string(&out_membership).unwrap(),
        concat!(
            "sampleA::r1\tTU000001\t1.000000\t1.000000\n",
            "sampleA::r2\tTU000001\t0.980198\t0.990000\n",
            "sampleB::r1\tTU000001\t0.600000\t0.600000\n",
            "sampleA::r3\tTU000002\t1.000000\t1.000000\n",
            "sampleB::r2\tTU000002\t1.000000\t1.000000\n",
            "sampleB::r3\tTU000003\t1.000000\t1.000000\n",
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
        fs::read_to_string(&out_tu_count).unwrap(),
        concat!(
            "tu_id,count\n",
            "TU000001,3\n",
            "TU000002,2\n",
            "TU000003,1\n",
        )
    );

    assert_eq!(
        fs::read_to_string(&out_sample_long).unwrap(),
        concat!(
            "tu_id\tsample\tcount\n",
            "TU000001\tsampleA\t2\n",
            "TU000001\tsampleB\t1\n",
            "TU000002\tsampleA\t1\n",
            "TU000002\tsampleB\t1\n",
            "TU000003\tsampleB\t1\n",
        )
    );

    assert_eq!(
        fs::read_to_string(&out_sample_matrix).unwrap(),
        concat!(
            "tu_id\tsampleA\tsampleB\n",
            "TU000001\t2\t1\n",
            "TU000002\t1\t1\n",
            "TU000003\t0\t1\n",
        )
    );

    assert_eq!(
        fs::read_to_string(&out_group_matrix).unwrap(),
        concat!(
            "tu_id\tcontrol\ttreated\n",
            "TU000001\t2\t1\n",
            "TU000002\t1\t1\n",
            "TU000003\t0\t1\n",
        )
    );

    assert_eq!(
        fs::read_to_string(&out_gene_count).unwrap(),
        concat!(
            "gene_id,count\n",
            "geneA,3\n",
            "geneB,3\n",
            "geneC,2\n",
            "geneD,1\n",
        )
    );

    assert_eq!(
        fs::read_to_string(&out_gene_sample_matrix).unwrap(),
        concat!(
            "gene_id\tsampleA\tsampleB\n",
            "geneA\t2\t1\n",
            "geneB\t2\t1\n",
            "geneC\t1\t1\n",
            "geneD\t0\t1\n",
        )
    );

    assert_eq!(
        fs::read_to_string(&out_gene_group_matrix).unwrap(),
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
        fs::read_to_string(&out_tu_count).unwrap(),
        concat!(
            "tu_id,count\n",
            "TU000001,3\n",
            "TU000002,2\n",
            "TU000003,1\n",
        )
    );
    assert_eq!(
        fs::read_to_string(&out_sample_matrix).unwrap(),
        concat!(
            "tu_id\tsampleA\tsampleB\n",
            "TU000001\t2\t1\n",
            "TU000002\t1\t1\n",
            "TU000003\t0\t1\n",
        )
    );
    assert_eq!(
        fs::read_to_string(&out_group_matrix).unwrap(),
        concat!(
            "tu_id\tcontrol\ttreated\n",
            "TU000001\t2\t1\n",
            "TU000002\t1\t1\n",
            "TU000003\t0\t1\n",
        )
    );

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
        fs::read_to_string(&out_sample_matrix).unwrap(),
        concat!("tu_id\tsampleA\tsampleB\n", "TU000001\t1\t1\n",)
    );
    assert!(!out_group_matrix.exists());

    let _ = fs::remove_dir_all(&tmp);
}
