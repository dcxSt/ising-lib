use ising_lib::lattice2d::{InitType, Lattice2d, SpinType, UpdateRule};
use ising_lib::monte_carlo_measurement::{MonteCarlo, MonteCarloParams};
use std::io::Write;

/// Mont Carlo sample of correlation with nearest neighbour for many
/// temperatures.
fn main() {
    let mut mu_vars = vec![];
    println!("Starting simulation");
    for temp in tqdm_rs::Tqdm::new((15..=100).step_by(5)) {
        print!("T={}, ",temp);
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
            n_runs: 500,
            flips_to_skip: 20000,
            samples_per_run: 50,
            flips_to_skip_intra_run: 1000,
        };

        // let mu_var = lattice.sample_neighbor_correlations(params);
        // mu_vars.push([mu_var[0],mu_var[1],temp as f64]);
        let mu_var = lattice.sample_average_magnetization(params);
        mu_vars.push([mu_var[0],mu_var[1],temp as f64])
    }
    println!("Done");
    let mut file = std::fs::File::create("data.txt").expect("Create failed");
    for idx in 0..mu_vars.len() {
        file.write_all(format!("{},{},{}\n",mu_vars[idx][0],mu_vars[idx][1],mu_vars[idx][2]).as_bytes()).expect("Write failed");
    }
    println!("Data written to file!");
}











