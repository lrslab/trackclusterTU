use std::fs;
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
fn trackclustertu_gff_to_bed_subcommand_converts_gene_features() {
    let tmp = unique_tmp_dir("trackclustertu_gff_to_bed_test");
    fs::create_dir_all(&tmp).unwrap();

    let gff_path = tmp.join("genes.gff3");
    let out_bed = tmp.join("genes.bed");

    fs::write(
        &gff_path,
        concat!(
            "##gff-version 3\n",
            "chr1\tTest\tgene\t11\t50\t.\t+\t.\tID=id1;Name=geneA;gene=ga;gene_biotype=protein_coding\n",
            "chr1\tTest\tCDS\t11\t50\t.\t+\t0\tParent=id1\n",
            "chr1\tTest\tgene\t61\t100\t.\t-\t.\tID=id2;locus_tag=b0002\n",
            "chr2\tTest\tgene\t5\t20\t.\t+\t.\tID=id3\n",
        ),
    )
    .unwrap();

    let exe = env!("CARGO_BIN_EXE_trackclustertu");
    let output = Command::new(exe)
        .args([
            "gff-to-bed",
            "--annotation-gff",
            gff_path.to_str().unwrap(),
            "--out-bed",
            out_bed.to_str().unwrap(),
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
        fs::read_to_string(&out_bed).unwrap(),
        concat!(
            "chr1\t10\t50\tgeneA\t0\t+\n",
            "chr1\t60\t100\tb0002\t0\t-\n",
            "chr2\t4\t20\tid3\t0\t+\n",
        )
    );

    let _ = fs::remove_dir_all(&tmp);
}
