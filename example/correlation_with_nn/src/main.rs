use ising_lib::lattice2d::{InitType, Lattice2d, SpinType, UpdateRule};
use ising_lib::monte_carlo_measurement::{MonteCarlo, MonteCarloParams};
// use std::io::Write;
use std::io;
use std::error::Error;
use csv;

/// Mont Carlo sample of correlation with nearest neighbour for many
/// temperatures.
fn main() -> Result<(), Box<dyn Error>> {
    // let mut mu_vars = vec![];
    let mut samples = vec![];
    println!("Starting simulation");

    let params = MonteCarloParams {
        n_runs: 25,
        flips_to_skip: 500_000,
        samples_per_run: 1,
        flips_to_skip_intra_run: 0,
    };

    for temp in tqdm_rs::Tqdm::new((15..=100).step_by(5)) {
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
    println!("Done");
    // TODO: write to csv instead of printing
    // let mut wtr = csv::Writer::from_writer(io::stdout());
    for (idx,temp) in (15..=100).step_by(5).enumerate() {
    //     file.write_all(format!("{},{},{}\n",mu_vars[idx][0],mu_vars[idx][1],mu_vars[idx][2]).as_bytes()).expect("Write failed");
        print!("{},",temp as f64 / 10.0);
        for i in 0..params.n_runs {
            print!("{},",samples[idx][i])
        }
        println!();
        // wtr.write_record(vec![samples[idx]])?;
    }
    // wtr.flush()?;
    // println!("Data written to file!");

    Ok(())
}




