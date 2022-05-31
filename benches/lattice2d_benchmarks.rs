use criterion::{criterion_group, criterion_main, Criterion};
use ising_lib::lattice2d::{Lattice2d};

// Criterion links:
// https://bheisler.github.io/criterion.rs/criterion/ 
// https://crates.io/crates/criterion 

// NOTE: keep results consistent, always run set lattice size 50x50
fn bench_flip_100_spins(c: &mut Criterion) {
    let mut lattice = Lattice2d::new_basic([50,50]);
    c.bench_function("try flip 100 spins", move |b| {
        b.iter(|| {
            lattice.update_n(100);
        })
    });
}


criterion_group!(benches, bench_flip_100_spins);
criterion_main!(benches);

