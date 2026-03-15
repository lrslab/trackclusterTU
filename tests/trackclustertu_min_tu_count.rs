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
fn trackclustertu_filters_low_support_tus() {
    let tmp = unique_tmp_dir("trackclustertu_min_tu_count_test");
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
            "--min-tu-count",
            "2",
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
            "chr1\t300\t400\tTU000002\t0\t+\n",
        )
    );

    let membership_text = fs::read_to_string(&out_membership).unwrap();
    assert_eq!(
        membership_text,
        concat!(
            "r1\tTU000001\t1.000000\t1.000000\n",
            "r2\tTU000001\t0.980198\t0.990000\n",
            "r3\tTU000001\t0.600000\t1.000000\n",
            "r4\tTU000002\t1.000000\t1.000000\n",
            "r5\tTU000002\t0.980198\t0.990000\n",
            "r6\tTU000002\t0.400000\t1.000000\n",
        )
    );

    let count_text = fs::read_to_string(&out_tu_count).unwrap();
    assert_eq!(
        count_text,
        concat!("tu_id,count\n", "TU000001,3\n", "TU000002,3\n",)
    );

    let _ = fs::remove_dir_all(&tmp);
}
