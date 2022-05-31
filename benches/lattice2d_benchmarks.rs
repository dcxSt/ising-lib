use criterion::{criterion_group, criterion_main, Criterion};
use ising_lib::lattice2d::Lattice2d;
use ising_lib::monte_carlo_measurement::{MonteCarlo, MonteCarloParams};

// Criterion links:
// https://bheisler.github.io/criterion.rs/criterion/
// https://crates.io/crates/criterion

// NOTE: keep results consistent, always run set lattice size 50x50
fn bench_flip_100_spins(c: &mut Criterion) {
    let mut lattice = Lattice2d::new_basic([50, 50]);
    c.bench_function("try flip 100 spins", move |b| {
        b.iter(|| {
            lattice.update_n(100);
        })
    });
}

// Monte Carlo benchmarks
fn bench_sample_energy(c: &mut Criterion) {
    let mut lattice = Lattice2d::new_basic([50, 50]);
    let params = MonteCarloParams {
        n_runs: 500,
        flips_to_skip: 50, // 1500_000,
        samples_per_run: 5,
        flips_to_skip_between_samples: 10,
    };
    c.bench_function("sample energy", move |b| {
        b.iter(|| {
            let _erg_samples: Vec<Vec<f64>> = lattice.sample_energy(&params);
        })
    });
}

fn bench_sample_energy_parallel(c: &mut Criterion) {
    let mut lattice = Lattice2d::new_basic([50, 50]);
    let params = MonteCarloParams {
        n_runs: 500,
        flips_to_skip: 50, // 1500_000,
        samples_per_run: 5,
        flips_to_skip_between_samples: 10,
    };
    c.bench_function("sample energy parallel", move |b| {
        b.iter(|| {
            let _erg_samples = lattice.sample_energy_parallel(&params);
        })
    });
}

criterion_group!(
    benches,
    bench_flip_100_spins,
    bench_sample_energy,
    bench_sample_energy_parallel
);
criterion_main!(benches);
