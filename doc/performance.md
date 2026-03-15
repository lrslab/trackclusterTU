# Performance (baseline benchmarks)

This project uses Criterion benchmarks to track performance of the hot paths:

- interval overlap + similarity scoring
- list-based intersection/union
- TU clustering (`cluster_tus`) on synthetic reads

## How to run

```bash
cargo bench
```

Or run a specific bench:

```bash
cargo bench --bench intersection
cargo bench --bench clustering
```

Criterion stores historical results under `target/criterion/`.

## Benchmark setup

### `benches/intersection.rs`

- `single_interval/*` runs 10k `(Interval, Interval)` pairs per iteration.
- `interval_lists/*` uses two sorted, non-overlapping lists of fixed-size intervals.

### `benches/clustering.rs`

- Synthetic reads on one `(contig=strand)` partition.
- Each “seed” produces many near-identical full-length reads (high `score1`) plus one contained truncation (exercises the `score2` attachment pass).
- Thresholds: `score1=0.95`, `score2=0.99`.

## Results (2026-01-09)

Environment:

- `rustc 1.91.0`, `cargo 1.91.0`
- `Linux 5.15.0-160-generic`, `x86_64` (96 CPUs reported)

Selected Criterion estimates:

- `single_interval/overlap_len`: ~15.9 µs
- `single_interval/score1_interval`: ~44.0 µs
- `single_interval/score2_interval`: ~84.6 µs

- `interval_lists/intersection_len/1000`: ~5.15 µs
- `interval_lists/union_len/1000`: ~6.11 µs
- `interval_lists/intersection_len/5000`: ~26.6 µs
- `interval_lists/union_len/5000`: ~31.7 µs

- `tu_clustering/cluster_tus/10000`: ~8.54 ms
- `tu_clustering/cluster_tus/100000`: ~89.0 ms

These are intended as a **baseline**; if performance regresses, check `target/criterion/` diffs and profile before optimizing.
