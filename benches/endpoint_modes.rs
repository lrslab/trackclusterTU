use criterion::{black_box, criterion_group, criterion_main, BenchmarkId, Criterion};
use trackclustertu::endpoint::{detect_endpoint_modes, TieDirection};

fn benchmark_endpoint_modes(c: &mut Criterion) {
    let mut group = c.benchmark_group("endpoint_modes");
    group.sample_size(10);

    for distinct_endpoints in [1_000usize, 2_500, 5_000] {
        // Three reads per distinct endpoint model a high-depth locus while a
        // 12 bp window exercises indexed range sums and bulk removal.
        let observations: Vec<u32> = (0..distinct_endpoints as u32)
            .flat_map(|coordinate| [coordinate, coordinate, coordinate])
            .collect();
        group.bench_with_input(
            BenchmarkId::from_parameter(distinct_endpoints),
            &observations,
            |b, observations| {
                b.iter(|| {
                    black_box(detect_endpoint_modes(
                        black_box(observations.iter().copied()),
                        12,
                        TieDirection::PreferLower,
                    ))
                });
            },
        );
    }

    group.finish();
}

criterion_group!(benches, benchmark_endpoint_modes);
criterion_main!(benches);
