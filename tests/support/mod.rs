#![allow(dead_code)]

pub fn membership_v1_projection(text: &str) -> String {
    let mut projected = String::new();
    for line in text.lines().filter(|line| !line.starts_with('#')) {
        let fields: Vec<&str> = line.split('\t').collect();
        assert!(
            fields.len() >= 4,
            "membership row has fewer than four columns"
        );
        projected.push_str(&fields[..4].join("\t"));
        projected.push('\n');
    }
    projected
}

pub fn count_v1_projection(text: &str) -> String {
    let mut projected = String::new();
    for line in text.lines().filter(|line| !line.starts_with('#')) {
        let fields: Vec<&str> = line.split(',').collect();
        assert!(fields.len() >= 2, "count row has fewer than two columns");
        projected.push_str(fields[0]);
        projected.push(',');
        projected.push_str(fields[1]);
        projected.push('\n');
    }
    projected
}

pub fn metric_matrix_total_projection(text: &str) -> String {
    let mut lines = text.lines().filter(|line| !line.starts_with('#'));
    let header: Vec<&str> = lines.next().expect("matrix header").split('\t').collect();
    let total_columns: Vec<usize> = header
        .iter()
        .enumerate()
        .filter_map(|(index, column)| column.ends_with(".total_count").then_some(index))
        .collect();
    assert!(
        !total_columns.is_empty(),
        "matrix has no total_count columns"
    );

    let mut projected = String::new();
    projected.push_str(header[0]);
    for &index in &total_columns {
        projected.push('\t');
        projected.push_str(header[index].trim_end_matches(".total_count"));
    }
    projected.push('\n');

    for line in lines {
        let fields: Vec<&str> = line.split('\t').collect();
        projected.push_str(fields[0]);
        for &index in &total_columns {
            projected.push('\t');
            projected.push_str(fields[index]);
        }
        projected.push('\n');
    }
    projected
}

pub fn metric_long_total_projection(text: &str) -> String {
    let mut lines = text.lines();
    let header: Vec<&str> = lines
        .next()
        .expect("long-table header")
        .split('\t')
        .collect();
    assert_eq!(header.get(4), Some(&"total_count"));
    let mut projected = format!("{}\t{}\tcount\n", header[0], header[1]);
    for line in lines {
        let fields: Vec<&str> = line.split('\t').collect();
        projected.push_str(&format!("{}\t{}\t{}\n", fields[0], fields[1], fields[4]));
    }
    projected
}
