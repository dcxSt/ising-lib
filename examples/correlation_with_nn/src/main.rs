use ising_lib::lattice2d::{InitType, Lattice2d, SpinType, UpdateRule};
use ising_lib::monte_carlo_measurement::{MonteCarlo, MonteCarloParams};
use std::error::Error;
use std::fs::File;
use std::io::Write;

/// Mont Carlo sample of correlation with nearest neighbour for many
/// temperatures.
fn main() -> Result<(), Box<dyn Error>> {
    // TODO: find better way to iter through float temps so we don't have this
    // hacky scaling by a factor of 10
    const RANGE_LOW: f64 = 0.5;
    const RANGE_HIGH: f64 = 6.0;
    const STEP_BY: f64 = 0.2;
    let mut samples = vec![];
    println!("Starting simulation...");

    let params = MonteCarloParams {
        n_runs: 25,
        flips_to_skip: 300_000, // 1500_000,
        samples_per_run: 10,
        flips_to_skip_between_samples: 30_000,
    };

    for temp in tqdm_rs::Tqdm::new(
        (((10.0 * RANGE_LOW) as usize)..=((10.0 * RANGE_HIGH) as usize))
            .step_by((10.0 * STEP_BY) as usize)
            .map(|t| t as f64 / 10.0),
    ) {
        print!("T={}, ", temp);
        // beta is 1/T (times boltzman constant, but we're ignoring that)
        let beta: f64 = 1.0 / temp as f64;
        println!("beta={}", beta);

        // Initiate a lattice
        let mut lattice = Lattice2d::new(
            [25, 25],
            UpdateRule::Metropolis,
            SpinType::SpinHalf,
            InitType::Random,
            1.0f64, // j interaction constant
            0.0f64, // h static field term
            beta,   // 1/Tkb
        );

        samples.push(lattice.sample_neighbor_correlations(&params));
    }
    println!("Done computing. Writing to file...");

    let mut file = File::create("./data.csv").unwrap();
    for (idx, temp) in (((10.0 * RANGE_LOW) as usize)..=((10.0 * RANGE_HIGH) as usize))
        .step_by((10.0 * STEP_BY) as usize)
        .map(|t| t as f64 / 10.0)
        .enumerate()
    {
        write!(&mut file, "{},", temp as f64 / 10.0).unwrap();
        for i in 0..(params.n_runs) {
            for j in 0..(params.samples_per_run) {
                write!(&mut file, "{}", samples[idx][i][j]).unwrap();
                if (i != params.n_runs-1) | (j  != params.samples_per_run - 1) {
                    write!(&mut file, ",").unwrap();
                }
            }
        }
        writeln!(&mut file, "").unwrap();
    }
    Ok(())
}
