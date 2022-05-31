use ising_lib::lattice2d::{InitType, Lattice2d, SpinType, UpdateRule}; 
use ising_lib::monte_carlo_measurement::{MonteCarlo, MonteCarloParams};

#[test]
fn integration_test_test() {
    let params = MonteCarloParams {
        n_runs: 5,
        flips_to_skip: 100,
        samples_per_run: 3,
        flips_to_skip_between_samples: 10,
    };
    let mut lattice = Lattice2d::new(
        [8,9],
        UpdateRule::Metropolis,
        SpinType::SpinHalf,
        InitType::Random,
        1.0f64, // j
        0.1f64, // h, static external field term
        2.4f64, // 1/TkB
    );
    let _nn_samples:Vec<Vec<f64>> = lattice.sample_neighbor_correlations(&params);
}
