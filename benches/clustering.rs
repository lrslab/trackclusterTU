use criterion::{black_box, criterion_group, criterion_main, BenchmarkId, Criterion, Throughput};

use trackclustertu::model::{Coord, Interval, Strand};
use trackclustertu::tu::{assign_reads_to_tus, cluster_tus, ReadRecord, Tu, TuClusteringOptions};

fn interval(start: u32, end: u32) -> Interval {
    Interval::new(Coord::new(start), Coord::new(end)).unwrap()
}

fn make_reads(total_reads: usize, reads_per_seed: usize) -> Vec<ReadRecord> {
    assert!(reads_per_seed >= 2);
    assert!(total_reads.is_multiple_of(reads_per_seed));

    let n_seeds = total_reads / reads_per_seed;
    let mut reads: Vec<ReadRecord> = Vec::with_capacity(total_reads);

    for seed in 0..n_seeds {
        let base_start = (seed as u32) * 10_000;
        let base_end = base_start + 1000;

        // Full-length reads (high score1 similarity).
        for i in 0..(reads_per_seed - 1) {
            let jitter = (i as u32) % 20;
            let start = base_start + jitter;
            let end = base_end + jitter;
            reads.push(ReadRecord {
                contig: "chr1".to_owned(),
                strand: Strand::Plus,
                interval: interval(start, end),
                id: format!("r{seed}_{i}"),
            });
        }

        // One shorter truncation to exercise the length-penalized score2 second-pass checks.
        reads.push(ReadRecord {
            contig: "chr1".to_owned(),
            strand: Strand::Plus,
            interval: interval(base_start + 100, base_end - 100),
            id: format!("r{seed}_trunc"),
        });
    }

    reads
}

/// Build one deeply covered, overlap-connected locus.
///
/// The ordinary benchmark above deliberately limits each disconnected seed to 100 reads. That
/// is representative of sparse inputs, but it does not exercise the high-depth case seen in
/// pooled bacterial direct-RNA data, where one connected region can contain many thousands of
/// nearly identical alignments. `endpoint_variants = 1` is the duplicate-heavy case; larger
/// values add realistic endpoint jitter while keeping every read in the same score1 component.
fn make_dense_reads(total_reads: usize, endpoint_variants: usize) -> Vec<ReadRecord> {
    assert!(endpoint_variants >= 1);

    (0..total_reads)
        .map(|i| {
            let jitter = (i % endpoint_variants) as u32;
            ReadRecord {
                contig: "chr1".to_owned(),
                strand: Strand::Plus,
                interval: interval(100 + jitter, 1100 + jitter),
                id: format!("dense_{i:08}"),
            }
        })
        .collect()
}

fn make_sparse_assignment_case(tu_count: usize, reads_per_tu: usize) -> (Vec<ReadRecord>, Vec<Tu>) {
    let mut reads = Vec::with_capacity(tu_count * reads_per_tu);
    let mut tus = Vec::with_capacity(tu_count);
    for tu_index in 0..tu_count {
        let start = (tu_index as u32) * 2_000;
        let tu_interval = interval(start, start + 1_000);
        tus.push(Tu {
            id: format!("TU{tu_index:08}"),
            contig: "chr1".to_owned(),
            strand: Strand::Plus,
            interval: tu_interval,
            rep_read_index: reads.len(),
        });
        for read_index in 0..reads_per_tu {
            reads.push(ReadRecord {
                contig: "chr1".to_owned(),
                strand: Strand::Plus,
                interval: tu_interval,
                id: format!("assignment_{tu_index:08}_{read_index}"),
            });
        }
    }
    (reads, tus)
}

fn bench_cluster_tus(c: &mut Criterion) {
    let mut group = c.benchmark_group("tu_clustering");
    group.sample_size(10);

    // Roughly: with 100 reads/seed, 10k reads -> 100 seeds; 100k reads -> 1000 seeds.
    for (n_reads, reads_per_seed) in [(10_000usize, 100usize), (100_000usize, 100usize)] {
        let reads = make_reads(n_reads, reads_per_seed);

        group.bench_with_input(
            BenchmarkId::new("cluster_tus", n_reads),
            &n_reads,
            |b, _| {
                b.iter(|| {
                    let result = cluster_tus(black_box(&reads), 0.95, 0.99).unwrap();
                    black_box(result.tus().len())
                })
            },
        );
    }

    group.finish();
}

fn bench_dense_hot_locus(c: &mut Criterion) {
    let mut group = c.benchmark_group("tu_clustering_dense_hot_locus");
    group.sample_size(10);

    // Two input sizes make quadratic regressions visible as roughly 4x growth when the read
    // count doubles. Twenty-one endpoint variants also cover a dense, jittered locus without
    // turning it into the disconnected 100-read components used by `bench_cluster_tus`.
    for (label, n_reads, endpoint_variants) in [
        ("identical", 5_000usize, 1usize),
        ("identical", 10_000usize, 1usize),
        ("jittered", 5_000usize, 21usize),
        ("jittered", 10_000usize, 21usize),
    ] {
        let reads = make_dense_reads(n_reads, endpoint_variants);
        group.throughput(Throughput::Elements(n_reads as u64));

        group.bench_with_input(BenchmarkId::new(label, n_reads), &n_reads, |b, _| {
            b.iter(|| {
                let result = cluster_tus(black_box(&reads), 0.95, 0.99).unwrap();
                black_box(result.tus().len())
            })
        });
    }

    group.finish();
}

fn bench_sparse_assignment_index(c: &mut Criterion) {
    let mut group = c.benchmark_group("tu_assignment_sparse_index");
    group.sample_size(10);

    for reads_per_tu in [1usize, 2usize] {
        let (reads, tus) = make_sparse_assignment_case(10_000, reads_per_tu);
        group.throughput(Throughput::Elements(reads.len() as u64));
        group.bench_with_input(
            BenchmarkId::new("reads", reads.len()),
            &reads.len(),
            |b, _| {
                b.iter(|| {
                    let assignments = assign_reads_to_tus(
                        black_box(&reads),
                        black_box(&tus),
                        0.95,
                        0.80,
                        TuClusteringOptions::default(),
                        0.02,
                        false,
                    )
                    .unwrap();
                    black_box(assignments.len())
                })
            },
        );
    }

    group.finish();
}

criterion_group!(
    benches,
    bench_cluster_tus,
    bench_dense_hot_locus,
    bench_sparse_assignment_index
);
criterion_main!(benches);
