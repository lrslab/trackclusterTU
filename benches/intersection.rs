use criterion::{black_box, criterion_group, criterion_main, BenchmarkId, Criterion};

use trackcluster_rs::interval::{intersection_len, union_len};
use trackcluster_rs::model::{Coord, Interval};
use trackcluster_rs::score::{score1_interval, score2_interval};

fn interval(start: u32, end: u32) -> Interval {
    Interval::new(Coord::new(start), Coord::new(end)).unwrap()
}

fn make_intervals(n: usize, offset: u32, len: u32, gap: u32) -> Vec<Interval> {
    (0..n)
        .map(|i| {
            let start = offset + (i as u32) * (len + gap);
            interval(start, start + len)
        })
        .collect()
}

fn bench_single_interval(c: &mut Criterion) {
    let mut group = c.benchmark_group("single_interval");

    let pairs: Vec<(Interval, Interval)> = (0..10_000)
        .map(|i| {
            let start = i * 10;
            let a = interval(start, start + 1000);
            let b = interval(start + (i % 20), start + 1000 + (i % 20));
            (a, b)
        })
        .collect();

    group.bench_function("overlap_len", |b| {
        b.iter(|| {
            let mut total = 0u64;
            for (a, other) in &pairs {
                total += a.overlap_len(*other) as u64;
            }
            black_box(total)
        })
    });

    group.bench_function("score1_interval", |b| {
        b.iter(|| {
            let mut total = 0.0f64;
            for (a, other) in &pairs {
                total += score1_interval(*a, *other);
            }
            black_box(total)
        })
    });

    group.bench_function("score2_interval", |b| {
        b.iter(|| {
            let mut total = 0.0f64;
            for (a, other) in &pairs {
                total += score2_interval(*a, *other);
            }
            black_box(total)
        })
    });

    group.finish();
}

fn bench_interval_lists(c: &mut Criterion) {
    let mut group = c.benchmark_group("interval_lists");

    for n in [10usize, 100, 1000, 5000] {
        let a = make_intervals(n, 0, 10, 10);
        let b = make_intervals(n, 5, 10, 10);

        group.bench_with_input(BenchmarkId::new("intersection_len", n), &n, |bch, _| {
            bch.iter(|| black_box(intersection_len(black_box(&a), black_box(&b))))
        });

        group.bench_with_input(BenchmarkId::new("union_len", n), &n, |bch, _| {
            bch.iter(|| black_box(union_len(black_box(&a), black_box(&b))))
        });
    }

    group.finish();
}

criterion_group!(benches, bench_single_interval, bench_interval_lists);
criterion_main!(benches);
