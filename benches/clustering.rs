use criterion::{black_box, criterion_group, criterion_main, BenchmarkId, Criterion};

use trackcluster_rs::model::{Coord, Interval, Strand};
use trackcluster_rs::tu::{cluster_tus, ReadRecord};

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

        // One contained truncation to exercise the score2 attachment pass.
        reads.push(ReadRecord {
            contig: "chr1".to_owned(),
            strand: Strand::Plus,
            interval: interval(base_start + 100, base_end - 100),
            id: format!("r{seed}_trunc"),
        });
    }

    reads
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
                    black_box(result.tus.len())
                })
            },
        );
    }

    group.finish();
}

criterion_group!(benches, bench_cluster_tus);
criterion_main!(benches);
