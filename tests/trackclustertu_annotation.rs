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
            "--in",
            reads_path.to_str().unwrap(),
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
            "#contig\tstrand\ttu_id\ttu_start\ttu_end\tgene_id\tgene_start\tgene_end\toverlap_bp\n",
            "chr1\t+\tTU000001\t100\t200\tgeneA\t90\t150\t50\n",
            "chr1\t+\tTU000001\t100\t200\tgeneB\t160\t220\t40\n",
            "chr1\t+\tTU000002\t300\t400\tgeneC\t290\t310\t10\n",
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
        gene_count_text,
        concat!("gene_id,count\n", "geneA,2\n", "geneB,2\n", "geneC,1\n",)
    );

    let _ = fs::remove_dir_all(&tmp);
}
