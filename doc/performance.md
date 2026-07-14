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

That local build directory is not a durable release record. The raw logs and
exact invocations behind the historical figures below are not tracked in this
repository, so treat them as descriptive reference points rather than release
gates. For a new release comparison, record the commit, exact command, thread
settings, hardware, operating system, and toolchain, and retain the raw
Criterion and wall-time/RSS logs with the external release record.

## Benchmark setup

### `benches/intersection.rs`

- `single_interval/*` runs 10k `(Interval, Interval)` pairs per iteration.
- `interval_lists/*` uses two sorted, non-overlapping lists of fixed-size intervals.

### `benches/clustering.rs`

- `tu_clustering` spreads reads across disconnected loci, with 100 reads per
  seed. Each seed produces near-identical full-length reads (high span Jaccard)
  plus one shorter truncation (exercises the overlap-over-longer second-pass
  checks).
- `tu_clustering_dense_hot_locus` keeps 5,000 or 10,000 reads in one
  overlap-connected locus. It covers both duplicate-heavy alignments and 21
  variants of endpoint jitter. Compare the 5,000- and 10,000-read estimates:
  roughly fourfold growth when reads double is a warning that per-family work
  has become quadratic.
- `tu_assignment_sparse_index` assigns 10,000 or 20,000 reads against 10,000
  separated TUs on one reference/strand partition. It guards the spatial index
  used to avoid the former all-reads-by-all-TUs comparison.
- The two clustering suites use span Jaccard `0.95` and overlap over longer
  `0.99`. The sparse-assignment suite uses `0.95` and the CLI-default `0.80`.

### `benches/endpoint_modes.rs`

- High-depth loci with 1,000, 2,500, and 5,000 distinct endpoint coordinates.
- Three molecules occur at every coordinate and a centered 12 bp mode window is used.
- The indexed implementation has O(k² log k) worst-case complexity for k distinct endpoints, replacing the previous cubic repeated rescans.
- Record wall time and peak RSS alongside Criterion estimates when establishing
  a release baseline, following the reproducibility requirements above.

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

A complete Criterion process measured with `/usr/bin/time -l` used a maximum
resident set size of 34,324,480 bytes (32.7 MiB). A second concurrently loaded
5,000-endpoint sample was noisy (10.792–18.286 ms). Because the exact invocation
and raw logs for these historical samples are not tracked here, use them only as
context; a release comparison requires one warm-up, three isolated measured
repeats, retained raw logs, and medians from the same hardware.

## Clustering and assignment scalability baseline (2026-07-14)

Environment: the same Apple M1 Max system above (`rustc 1.90.0`, 10 logical
CPUs, 64 GiB RAM). Criterion sample size was 10.

- `identical/5000`: 1.999–2.034 ms
- `identical/10000`: 3.978–4.013 ms
- `jittered/5000`: 2.785–2.809 ms
- `jittered/10000`: 5.555–5.602 ms
- `tu_assignment_sparse_index/reads/10000`: 4.177–4.273 ms
- `tu_assignment_sparse_index/reads/20000`: 5.233–5.584 ms

Within these two synthetic cases, whose endpoint diversity stays fixed at one or
21 variants, doubling reads from 5,000 to 10,000 took about 2x. This is not a
general worst-case linearity claim. The sparse-assignment comparison likewise
keeps the TU count fixed at 10,000 while doubling reads per TU.

The v0.2.0 implementation took about 0.98 s at 8,000 duplicate-heavy reads,
3.72 s at 16,000, and 15.63 s at 32,000 in a direct one-thread CLI measurement,
showing an approximately quadratic curve for that historical input. Its exact
command, generated input, and raw log are not tracked in this repository, so the
comparison is contextual evidence rather than a reproducible release gate.
