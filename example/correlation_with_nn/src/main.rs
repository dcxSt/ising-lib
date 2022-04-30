use ising_lib::lattice2d::{InitType, Lattice2d, SpinType, UpdateRule};
use ising_lib::monte_carlo_measurement::{MonteCarlo, MonteCarloParams};

/// Mont Carlo sample of correlation with nearest neighbour for many
/// temperatures.
fn main() {
    for temp in (15..=100).step_by(5) {
        // beta is 1/T (times boltzman constant, but we're ignoring that)
        let beta: f64 = 100.0 / temp as f64;

        // initiate the lattice
        let mut lattice = Lattice2d::new(
            [100, 100],
            UpdateRule::Metropolis,
            SpinType::SpinHalf,
            InitType::Random,
            0.0,
            0.0,
            beta,
        );

        let params = MonteCarloParams {
            n_per_sample: 25usize,
            lattice_dims: [100, 100],
            flips_to_skip: 5000,
            samples_per_measurement: 5,
            flips_to_skip_intra_measurment: 1000,
        };
        // let est_corr_nn =
    }
}
