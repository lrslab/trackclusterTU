# Performance (baseline benchmarks)

This project uses Criterion benchmarks to track performance of the hot paths:

- interval overlap + similarity scoring
- list-based intersection/union
- TU clustering (`cluster_tus`) on synthetic reads
- boundary-mode detection at high-depth loci with thousands of distinct endpoints

## How to run

```bash
cargo bench
```

Or run a specific bench:

```bash
cargo bench --bench intersection
cargo bench --bench clustering
cargo bench --bench endpoint_modes
```

Criterion stores historical results under `target/criterion/`.

## Benchmark setup

### `benches/intersection.rs`

- `single_interval/*` runs 10k `(Interval, Interval)` pairs per iteration.
- `interval_lists/*` uses two sorted, non-overlapping lists of fixed-size intervals.

### `benches/clustering.rs`

- Synthetic reads on one `(contig=strand)` partition.
- Each “seed” produces many near-identical full-length reads (high span Jaccard) plus one shorter truncation (exercises the overlap-over-longer second-pass checks).
- Thresholds: span Jaccard `0.95`, overlap over longer `0.99`.

### `benches/endpoint_modes.rs`

- High-depth loci with 1,000, 2,500, and 5,000 distinct endpoint coordinates.
- Three molecules occur at every coordinate and a centered 12 bp mode window is used.
- The indexed implementation has O(k² log k) worst-case complexity for k distinct endpoints, replacing the previous cubic repeated rescans.
- Record wall time and peak RSS alongside Criterion estimates when establishing a release baseline, and retain the exact command, hardware, and raw logs with that benchmark record.

## Results (2026-01-09)

Environment:

- `rustc 1.91.0`, `cargo 1.91.0`
- `Linux 5.15.0-160-generic`, `x86_64` (96 CPUs reported)

Selected Criterion estimates:

- `single_interval/overlap_len`: ~15.9 µs
- `single_interval/score1_interval`: ~44.0 µs (legacy Rust API name for span Jaccard)
- `single_interval/score2_interval`: ~84.6 µs (legacy Rust API name for overlap over longer)

- `interval_lists/intersection_len/1000`: ~5.15 µs
- `interval_lists/union_len/1000`: ~6.11 µs
- `interval_lists/intersection_len/5000`: ~26.6 µs
- `interval_lists/union_len/5000`: ~31.7 µs

- `tu_clustering/cluster_tus/10000`: ~8.54 ms
- `tu_clustering/cluster_tus/100000`: ~89.0 ms

These are intended as a **baseline**; if performance regresses, check `target/criterion/` diffs and profile before optimizing.

## Endpoint-mode baseline (2026-07-10)

Environment:

- `rustc 1.90.0`
- macOS 15.7.1, Apple M1 Max (10 logical CPUs), 64 GiB RAM
- Criterion sample size 10; Plotters backend

Initial Criterion estimates:

- `endpoint_modes/1000`: 517.75–530.03 µs
- `endpoint_modes/2500`: 2.7025–2.7129 ms
- `endpoint_modes/5000`: 10.070–10.244 ms

A complete Criterion process measured with `/usr/bin/time -l` used a maximum resident set size of 34,324,480 bytes (32.7 MiB). A second concurrently loaded 5,000-endpoint sample was noisy (10.792–18.286 ms), so release comparisons must use one warm-up, three isolated measured repeats, retained raw logs, and medians from the same hardware.
