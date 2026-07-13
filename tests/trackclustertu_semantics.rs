use std::collections::BTreeMap;
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

fn run_cluster(args: &[&str]) -> std::process::Output {
    Command::new(env!("CARGO_BIN_EXE_trackclustertu"))
        .arg("cluster")
        .args(args)
        .output()
        .unwrap()
}

fn bed_ids(path: &std::path::Path) -> Vec<String> {
    fs::read_to_string(path)
        .unwrap()
        .lines()
        .map(|line| line.split('\t').nth(3).unwrap().to_owned())
        .collect()
}

#[test]
fn stable_ids_survive_filtering_and_sequential_mode_is_auditable() {
    let tmp = unique_tmp_dir("trackclustertu_stable_ids");
    fs::create_dir_all(&tmp).unwrap();
    let reads = tmp.join("reads.bed");
    fs::write(
        &reads,
        concat!(
            "chr1\t10\t20\tsingleton\t0\t+\n",
            "chr1\t100\t200\tkept1\t0\t+\n",
            "chr1\t100\t200\tkept2\t0\t+\n",
        ),
    )
    .unwrap();

    let all_dir = tmp.join("all");
    let filtered_dir = tmp.join("filtered");
    let sequential_dir = tmp.join("sequential");
    let all = run_cluster(&[
        "--in",
        reads.to_str().unwrap(),
        "--format",
        "bed6",
        "--out-dir",
        all_dir.to_str().unwrap(),
    ]);
    let filtered = run_cluster(&[
        "--in",
        reads.to_str().unwrap(),
        "--format",
        "bed6",
        "--min-tu-count",
        "2",
        "--out-dir",
        filtered_dir.to_str().unwrap(),
    ]);
    let sequential = run_cluster(&[
        "--in",
        reads.to_str().unwrap(),
        "--format",
        "bed6",
        "--min-tu-count",
        "2",
        "--tu-id-style",
        "sequential",
        "--out-dir",
        sequential_dir.to_str().unwrap(),
    ]);
    for output in [&all, &filtered, &sequential] {
        assert!(
            output.status.success(),
            "stderr:\n{}",
            String::from_utf8_lossy(&output.stderr)
        );
    }

    let all_ids = bed_ids(&all_dir.join("tus.bed"));
    let filtered_ids = bed_ids(&filtered_dir.join("tus.bed"));
    assert_eq!(
        all_ids,
        vec![
            "TUg_63687231_10_20_p".to_owned(),
            "TUg_63687231_100_200_p".to_owned(),
        ]
    );
    assert_eq!(filtered_ids, vec![all_ids[1].clone()]);
    assert_eq!(bed_ids(&sequential_dir.join("tus.bed")), vec!["TU000001"]);

    let mapping = fs::read_to_string(sequential_dir.join("tu_id_map.tsv")).unwrap();
    assert!(mapping.starts_with("#trackclustertu_tu_id_map_schema=v1\n"));
    assert!(mapping.contains("TU000001\tTUg_63687231_100_200_p\tTU000001\tchr1\t+\t100\t200\n"));
    let _ = fs::remove_dir_all(tmp);
}

#[test]
fn semantics_cover_directionality_multilabels_antisense_and_gff3_escaping() {
    let tmp = unique_tmp_dir("trackclustertu_semantics");
    fs::create_dir_all(&tmp).unwrap();
    let reads = tmp.join("reads.bed");
    let genes = tmp.join("genes.bed");
    let out_dir = tmp.join("out");
    fs::write(
        &reads,
        concat!(
            "chr1\t90\t210\tp_full\t0\t+\n",
            "chr1\t110\t210\tp_alt_start\t0\t+\n",
            "chr1\t90\t190\tp_alt_term\t0\t+\n",
            "chr1\t120\t180\tp_processed\t0\t+\n",
            "chr1\t90\t310\tp_readthrough\t0\t+\n",
            "chr1\t410\t490\tp_antisense\t0\t+\n",
            "chr1\t600\t650\tp_intergenic\t0\t+\n",
            "chr1\t690\t810\tm_full\t0\t-\n",
            "chr1\t690\t790\tm_alt_start\t0\t-\n",
            "chr1\t710\t810\tm_alt_term\t0\t-\n",
            "chr1\t720\t780\tm_processed\t0\t-\n",
            "chr1\t690\t910\tm_readthrough\t0\t-\n",
            "chr 2\t990\t1110\tweird\t0\t+\n",
        ),
    )
    .unwrap();
    fs::write(
        &genes,
        concat!(
            "chr1\t100\t200\tgeneA\t0\t+\n",
            "chr1\t210\t300\tgeneB\t0\t+\n",
            "chr1\t250\t280\tanti_readthrough\t0\t-\n",
            "chr1\t400\t500\tanti_only\t0\t-\n",
            "chr1\t700\t800\tgeneC\t0\t-\n",
            "chr1\t820\t900\tgeneD\t0\t-\n",
            "chr 2\t1000\t1100\tgene;weird=1\t0\t+\n",
        ),
    )
    .unwrap();

    let output = run_cluster(&[
        "--in",
        reads.to_str().unwrap(),
        "--format",
        "bed6",
        "--annotation-bed",
        genes.to_str().unwrap(),
        "--span-jaccard-threshold",
        "1",
        "--overlap-over-longer-threshold",
        "1",
        "--skip-score2-attachment",
        "--out-dir",
        out_dir.to_str().unwrap(),
    ]);
    assert!(
        output.status.success(),
        "stderr:\n{}",
        String::from_utf8_lossy(&output.stderr)
    );

    let semantics_text = fs::read_to_string(out_dir.join("tu_semantics.tsv")).unwrap();
    assert!(semantics_text.contains("#classification_order=alternative_start,alternative_termination,readthrough,antisense,intergenic,probable_processing_product\n"));
    let mut semantics: BTreeMap<(String, String, String), (String, String)> = BTreeMap::new();
    for line in semantics_text
        .lines()
        .filter(|line| !line.starts_with('#'))
        .skip(1)
    {
        let fields: Vec<&str> = line.split('\t').collect();
        semantics.insert(
            (
                fields[1].to_owned(),
                fields[2].to_owned(),
                format!("{}-{}", fields[3], fields[4]),
            ),
            (fields[5].to_owned(), fields[9].to_owned()),
        );
    }
    assert_eq!(
        semantics[&("chr1".into(), "+".into(), "90-210".into())].1,
        "alternative_start,alternative_termination"
    );
    assert_eq!(
        semantics[&("chr1".into(), "+".into(), "120-180".into())].1,
        "probable_processing_product"
    );
    assert_eq!(
        semantics[&("chr1".into(), "+".into(), "90-310".into())],
        ("geneA,geneB".into(), "readthrough,antisense".into())
    );
    assert_eq!(
        semantics[&("chr1".into(), "+".into(), "410-490".into())].1,
        "antisense"
    );
    assert_eq!(
        semantics[&("chr1".into(), "+".into(), "600-650".into())].1,
        "intergenic"
    );
    assert_eq!(
        semantics[&("chr1".into(), "-".into(), "690-910".into())],
        ("geneD,geneC".into(), "readthrough".into())
    );
    assert_eq!(
        semantics[&("chr1".into(), "-".into(), "690-810".into())].1,
        "alternative_start,alternative_termination"
    );
    assert_eq!(
        semantics[&("chr1".into(), "-".into(), "720-780".into())].1,
        "probable_processing_product"
    );
    assert_eq!(
        semantics[&("chr 2".into(), "+".into(), "990-1110".into())].0,
        "gene%3Bweird%3D1"
    );

    let relationships = fs::read_to_string(out_dir.join("tu_gene.tsv")).unwrap();
    assert!(relationships.contains("\tantisense\t1\tanti_readthrough\t"));
    assert!(relationships.contains("\tantisense\t1\tanti_only\t"));

    let gff3 = fs::read_to_string(out_dir.join("tus.gff3")).unwrap();
    assert!(gff3.starts_with("##gff-version 3\n"));
    assert!(gff3.contains("chr%202\ttrackclustertu\ttranscript\t991\t1110"));
    assert!(gff3.contains("gene_context=gene%3Bweird%3D1"));
    assert!(gff3.contains("gene_id=gene%3Bweird%3D1"));
    assert!(!gff3.contains("gene;weird=1"));

    let gene_counts = fs::read_to_string(out_dir.join("gene_count.csv")).unwrap();
    assert!(gene_counts.starts_with(
        "#count_semantics=nonexclusive_same_strand_tu_overlap;each_TU_assignment_is_counted_once_for_every_qualifying_gene;gene_totals_may_exceed_read_and_TU_totals;antisense_relationships_are_not_counted\n"
    ));
    assert!(gene_counts.lines().any(|line| line.starts_with("geneA,5,")));
    assert!(gene_counts.lines().any(|line| line.starts_with("geneB,1,")));
    assert!(gene_counts
        .lines()
        .any(|line| line.starts_with("anti_only,0,")));
    let _ = fs::remove_dir_all(tmp);
}

#[test]
fn overlap_policy_filters_relationships_and_invalid_fractions_publish_nothing() {
    let tmp = unique_tmp_dir("trackclustertu_gene_policy");
    fs::create_dir_all(&tmp).unwrap();
    let reads = tmp.join("reads.bed");
    let genes = tmp.join("genes.bed");
    fs::write(&reads, "chr1\t0\t100\tr1\t0\t+\n").unwrap();
    fs::write(&genes, "chr1\t90\t200\ttiny_overlap\t0\t+\n").unwrap();
    let out_dir = tmp.join("filtered");
    let filtered = run_cluster(&[
        "--in",
        reads.to_str().unwrap(),
        "--format",
        "bed6",
        "--annotation-bed",
        genes.to_str().unwrap(),
        "--gene-min-overlap-bp",
        "20",
        "--gene-min-tu-fraction",
        "0.05",
        "--gene-min-gene-fraction",
        "0.05",
        "--out-dir",
        out_dir.to_str().unwrap(),
    ]);
    assert!(filtered.status.success());
    let semantics = fs::read_to_string(out_dir.join("tu_semantics.tsv")).unwrap();
    let data_row = semantics
        .lines()
        .filter(|line| !line.starts_with('#'))
        .nth(1)
        .unwrap();
    assert!(data_row.ends_with("\tintergenic"));
    assert!(!fs::read_to_string(out_dir.join("tu_gene.tsv"))
        .unwrap()
        .lines()
        .any(|line| !line.starts_with('#')));

    let protected_output = tmp.join("protected.bed");
    fs::write(&protected_output, "previous-success\n").unwrap();
    let invalid = run_cluster(&[
        "--in",
        reads.to_str().unwrap(),
        "--format",
        "bed6",
        "--annotation-bed",
        genes.to_str().unwrap(),
        "--gene-min-tu-fraction",
        "1.1",
        "--out-tu",
        protected_output.to_str().unwrap(),
    ]);
    assert!(!invalid.status.success());
    assert!(String::from_utf8_lossy(&invalid.stderr)
        .contains("--gene-min-tu-fraction must be finite and between 0 and 1"));
    assert_eq!(
        fs::read_to_string(&protected_output).unwrap(),
        "previous-success\n"
    );
    let _ = fs::remove_dir_all(tmp);
}
