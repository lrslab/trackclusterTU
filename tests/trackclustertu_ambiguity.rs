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

fn run_cluster(input: &std::path::Path, out_dir: &std::path::Path) -> std::process::Output {
    Command::new(env!("CARGO_BIN_EXE_trackclustertu"))
        .args([
            "cluster",
            "--tu-id-style",
            "sequential",
            "--in",
            input.to_str().unwrap(),
            "--format",
            "bed6",
            "--span-jaccard-threshold",
            "0.95",
            "--overlap-over-longer-threshold",
            "0.90",
            "--three-prime-tolerance-bp",
            "5",
            "--ambiguity-margin",
            "0.02",
            "--fractional-assignment",
            "--min-tu-count",
            "2",
            "--out-dir",
            out_dir.to_str().unwrap(),
        ])
        .output()
        .unwrap()
}

#[test]
fn reports_ambiguity_fractions_unassigned_reads_and_endpoint_evidence_deterministically() {
    let tmp = unique_tmp_dir("trackclustertu_ambiguity");
    fs::create_dir_all(&tmp).unwrap();
    let canonical = tmp.join("canonical.bed");
    let permuted = tmp.join("permuted.bed");
    let records = [
        "chr1\t0\t100\ta1\t0\t+",
        "chr1\t0\t100\ta2\t0\t+",
        "chr1\t10\t110\tb1\t0\t+",
        "chr1\t10\t110\tb2\t0\t+",
        "chr1\t5\t105\tquery\t0\t+",
        "chr1\t500\t600\tfiltered\t0\t+",
    ];
    fs::write(&canonical, records.join("\n") + "\n").unwrap();
    fs::write(
        &permuted,
        records.iter().rev().copied().collect::<Vec<_>>().join("\n") + "\n",
    )
    .unwrap();

    let first_dir = tmp.join("first");
    let second_dir = tmp.join("second");
    let first = run_cluster(&canonical, &first_dir);
    let second = run_cluster(&permuted, &second_dir);
    assert!(
        first.status.success(),
        "stderr:\n{}",
        String::from_utf8_lossy(&first.stderr)
    );
    assert!(
        second.status.success(),
        "stderr:\n{}",
        String::from_utf8_lossy(&second.stderr)
    );

    let membership = fs::read_to_string(first_dir.join("membership.tsv")).unwrap();
    assert_eq!(
        membership,
        fs::read_to_string(second_dir.join("membership.tsv")).unwrap()
    );
    assert!(membership.starts_with("#trackclustertu_membership_schema=v2\n#columns="));

    let query_fields: Vec<&str> = membership
        .lines()
        .find(|line| line.starts_with("query\t"))
        .unwrap()
        .split('\t')
        .collect();
    assert_eq!(query_fields[1], ".");
    assert_eq!(query_fields[5], "ambiguous");
    assert_eq!(query_fields[6], "TU000001");
    assert_eq!(query_fields[7], "TU000002");
    assert_eq!(query_fields[16], "0.000000");
    let fractions: Vec<(&str, f64)> = query_fields[18]
        .split(',')
        .map(|item| {
            let (tu_id, weight) = item.split_once(':').unwrap();
            (tu_id, weight.parse().unwrap())
        })
        .collect();
    assert_eq!(fractions.len(), 2);
    assert_eq!(
        fractions.iter().map(|(_, weight)| *weight).sum::<f64>(),
        1.0
    );

    let filtered_fields: Vec<&str> = membership
        .lines()
        .find(|line| line.starts_with("filtered\t"))
        .unwrap()
        .split('\t')
        .collect();
    assert_eq!(filtered_fields[1], ".");
    assert_eq!(filtered_fields[5], "unassigned");
    assert_eq!(filtered_fields[18], ".");

    assert_eq!(
        fs::read_to_string(first_dir.join("tu_count.csv")).unwrap(),
        concat!(
            "tu_id,count,unique_count,full_length_evidence_count,total_count,fractional_count\n",
            "TU000001,2,2,0,2,2.5\n",
            "TU000002,2,2,0,2,2.5\n",
        )
    );

    let endpoint_stats = fs::read_to_string(first_dir.join("tu_endpoint_stats.tsv")).unwrap();
    assert!(endpoint_stats.starts_with(
        "#trackclustertu_tu_endpoint_stats_schema=v1\ntu_id\tsupport\tfive_prime_consensus"
    ));
    assert!(endpoint_stats
        .lines()
        .any(|line| line == "TU000001\t3\t0\t100\t0\t5\t5\t100\t105\t5\ta1"));
    assert!(endpoint_stats
        .lines()
        .any(|line| line == "TU000002\t2\t10\t110\t10\t10\t0\t110\t110\t0\tb1"));

    let _ = fs::remove_dir_all(tmp);
}
