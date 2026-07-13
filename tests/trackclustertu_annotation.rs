use std::fs;
use std::io::Write;
use std::process::Command;
use std::time::{SystemTime, UNIX_EPOCH};

mod support;
use support::count_v1_projection;

fn unique_tmp_dir(prefix: &str) -> std::path::PathBuf {
    let nanos = SystemTime::now()
        .duration_since(UNIX_EPOCH)
        .expect("clock")
        .as_nanos();
    std::env::temp_dir().join(format!("{prefix}_{nanos}"))
}

#[test]
fn trackclustertu_writes_gene_anchored_outputs() {
    let tmp = unique_tmp_dir("trackclustertu_annot_test");
    fs::create_dir_all(&tmp).unwrap();

    let reads_path = tmp.join("reads.bed");
    let genes_path = tmp.join("genes.bed");
    let out_tu = tmp.join("tu.bed");
    let out_membership = tmp.join("membership.tsv");
    let out_tu_gene = tmp.join("tu_gene.tsv");
    let out_tu_bed12 = tmp.join("tu_anchored.bed12");
    let out_gene_count = tmp.join("gene_count.csv");

    let mut reads = fs::File::create(&reads_path).unwrap();
    writeln!(reads, "chr1\t100\t200\tr1\t0\t+").unwrap();
    writeln!(reads, "chr1\t101\t201\tr2\t0\t+").unwrap();
    writeln!(reads, "chr1\t300\t400\tr3\t0\t+").unwrap();

    let mut genes = fs::File::create(&genes_path).unwrap();
    writeln!(genes, "chr1\t90\t150\tgeneA\t0\t+").unwrap();
    writeln!(genes, "chr1\t160\t220\tgeneB\t0\t+").unwrap();
    writeln!(genes, "chr1\t290\t310\tgeneC\t0\t+").unwrap();

    let exe = env!("CARGO_BIN_EXE_trackclustertu");
    let output = Command::new(exe)
        .args([
            "cluster",
            "--tu-id-style",
            "sequential",
            "--in",
            reads_path.to_str().unwrap(),
            "--format",
            "bed6",
            "--span-jaccard-threshold",
            "0.95",
            "--overlap-over-longer-threshold",
            "0.99",
            "--out-tu",
            out_tu.to_str().unwrap(),
            "--out-membership",
            out_membership.to_str().unwrap(),
            "--annotation-bed",
            genes_path.to_str().unwrap(),
            "--out-tu-gene",
            out_tu_gene.to_str().unwrap(),
            "--out-tu-bed12",
            out_tu_bed12.to_str().unwrap(),
            "--out-gene-count",
            out_gene_count.to_str().unwrap(),
        ])
        .output()
        .unwrap();

    assert!(
        output.status.success(),
        "stdout:\n{}\nstderr:\n{}",
        String::from_utf8_lossy(&output.stdout),
        String::from_utf8_lossy(&output.stderr)
    );

    let tu_gene_text = fs::read_to_string(&out_tu_gene).unwrap();
    assert_eq!(
        tu_gene_text,
        concat!(
            "#trackclustertu_tu_gene_schema=v2\n",
            "#gene_context_order=transcription_direction\n",
            "#gene_min_overlap_bp=1\n",
            "#gene_min_tu_fraction=0\n",
            "#gene_min_gene_fraction=0\n",
            "#contig\tstrand\ttu_id\ttu_start\ttu_end\trelation\tcontext_order\tgene_id\tgene_start\tgene_end\toverlap_bp\ttu_fraction\tgene_fraction\n",
            "chr1\t+\tTU000001\t100\t200\tsame_strand\t1\tgeneA\t90\t150\t50\t0.500000\t0.833333\n",
            "chr1\t+\tTU000001\t100\t200\tsame_strand\t2\tgeneB\t160\t220\t40\t0.400000\t0.666667\n",
            "chr1\t+\tTU000002\t300\t400\tsame_strand\t1\tgeneC\t290\t310\t10\t0.100000\t0.500000\n",
        )
    );

    let bed12_text = fs::read_to_string(&out_tu_bed12).unwrap();
    assert_eq!(
        bed12_text,
        concat!(
            "chr1\t100\t200\tTU000001\t0\t+\t100\t200\t0\t2\t50,40,\t0,60,\tr1,r2,|2\tgeneA,geneB\n",
            "chr1\t300\t400\tTU000002\t0\t+\t300\t400\t0\t1\t10,\t0,\tr3,|1\tgeneC\n",
        )
    );

    let gene_count_text = fs::read_to_string(&out_gene_count).unwrap();
    assert_eq!(
        count_v1_projection(&gene_count_text),
        concat!("gene_id,count\n", "geneA,2\n", "geneB,2\n", "geneC,1\n",)
    );

    let _ = fs::remove_dir_all(&tmp);
}

#[test]
fn trackclustertu_merges_overlapping_gene_blocks_in_anchored_bed12() {
    let tmp = unique_tmp_dir("trackclustertu_overlapping_gene_blocks");
    fs::create_dir_all(&tmp).unwrap();
    let reads_path = tmp.join("reads.bed");
    let genes_path = tmp.join("genes.bed");
    let out_bed12 = tmp.join("anchored.bed12");
    fs::write(&reads_path, "chr1\t100\t200\tr1\t0\t+\n").unwrap();
    fs::write(
        &genes_path,
        concat!(
            "chr1\t90\t170\tgeneA\t0\t+\n",
            "chr1\t150\t220\tgeneB\t0\t+\n",
        ),
    )
    .unwrap();

    let output = Command::new(env!("CARGO_BIN_EXE_trackclustertu"))
        .args([
            "cluster",
            "--tu-id-style",
            "sequential",
            "--in",
            reads_path.to_str().unwrap(),
            "--format",
            "bed6",
            "--annotation-bed",
            genes_path.to_str().unwrap(),
            "--out-tu-bed12",
            out_bed12.to_str().unwrap(),
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
        fs::read_to_string(out_bed12).unwrap(),
        "chr1\t100\t200\tTU000001\t0\t+\t100\t200\t0\t1\t100,\t0,\tr1,|1\tgeneA,geneB\n"
    );
    let _ = fs::remove_dir_all(tmp);
}
