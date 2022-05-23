use ising_lib::lattice2d::{InitType, Lattice2d, SpinType, UpdateRule};
use ising_lib::monte_carlo_measurement::{MonteCarlo, MonteCarloParams};
// use std::io::Write;
use std::io;
use std::error::Error;
// use csv;
use std::env;
use std::fs::File;
use std::io::Write;

/// Mont Carlo sample of correlation with nearest neighbour for many
/// temperatures.
fn main() -> Result<(), Box<dyn Error>> {
    // TODO: find better way to iter through float temps so we don't have this
    // hacky scaling by a factor of 10
    let RANGE_LOW:usize = 5;
    let RANGE_HIGH:usize = 80;
    let STEP_BY:usize = 1;
    // let mut mu_vars = vec![];
    let mut samples = vec![];
    println!("Starting simulation");

    let params = MonteCarloParams {
        n_runs: 35,
        flips_to_skip: 1500_000,
        samples_per_run: 1,
        flips_to_skip_intra_run: 0,
    };

    for temp in tqdm_rs::Tqdm::new((RANGE_LOW..=RANGE_HIGH).step_by(STEP_BY)) {
        print!("T={}, ",temp);
        // beta is 1/T (times boltzman constant, but we're ignoring that)
        let beta: f64 = 10.0 / temp as f64;
        println!("beta={}",beta);

        // initiate the lattice
        let mut lattice = Lattice2d::new(
            [25, 25],
            UpdateRule::Metropolis,
            SpinType::SpinHalf,
            InitType::Random,
            1.0f64, // j interaction constant
            0.0f64, // h static field term
            beta, // 1/Tkb
        );

        let sample_points:Vec<f64> = lattice.sample_neighbor_correlations(&params);
        samples.push(sample_points)
        // let mu_var = lattice.sample_neighbor_correlations(params);
        // mu_vars.push([mu_var[0],mu_var[1],temp as f64]);
        // let mu_var = lattice.sample_average_magnetization(params);
        // mu_vars.push([mu_var[0],mu_var[1],temp as f64])
    }
    println!("Done computing. Writing to file...");

    let mut file = File::create("./data.csv").unwrap();
    for (idx,temp) in (RANGE_LOW..=RANGE_HIGH).step_by(STEP_BY).enumerate() {
        write!(&mut file, "{},", temp as f64 / 10.0).unwrap();
        for i in 0..(params.n_runs-1) {
            write!(&mut file, "{},", samples[idx][i]).unwrap();
        }
        writeln!(&mut file, "{}", samples[idx][params.n_runs-1]).unwrap();
    }

    Ok(())
}




