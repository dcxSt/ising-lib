//! The MonteCarlo type samples quantities on networks that are 
//! estimated over multiple runs, (such as the expected value of the 
//! magnetization squared) 

use crate::lattice2d::*;
use crate::measurement::Measurement; 
use std::thread;

/// Parameters for monte carlo sampling
pub struct MonteCarloParams {
    pub n_runs: usize,                        // number of dry runs
    pub flips_to_skip: usize,                 // skip flips for system to cool
    pub samples_per_run: usize,               // number of samples to make in each run
    pub flips_to_skip_between_samples: usize, // number of flips to skip between each sample from the same run
}

/// The measurement trait samples quantities across lattices and graphs
pub trait MonteCarlo {
    /// Calculates and returns basic metrics by Monto Carlo sampling
    /// - Energy fluctuations
    /// - Avg Magnetic Susceptibility
    fn sample_energy_parallel(&mut self, params: &MonteCarloParams) -> Vec<Vec<f64>>;
    fn sample_energy(&mut self, params: &MonteCarloParams) -> Vec<Vec<f64>>;
    fn sample_neighbor_correlations_parallel(&mut self, params: &MonteCarloParams) -> Vec<Vec<f64>>;
    fn sample_neighbor_correlations(&mut self, params: &MonteCarloParams) -> Vec<Vec<f64>>;
    fn sample_magnetization_parallel(&mut self, params: &MonteCarloParams) -> Vec<Vec<f64>>;
    fn sample_magnetization(&mut self, params: &MonteCarloParams) -> Vec<Vec<f64>>;
    // TODO: implement below function
    // fn sample_estimate_all_metrics(&self , params:MonteCarloParams) -> Vec;
}

/// Implements the measurement trait for the Lattice2d type
impl MonteCarlo for Lattice2d {
    // fn sample_energy_fluctuations(&self , params:MonteCarloParams) -> [f64;3] {
    //     // initiate energy vector
    //     // for 0..n_samples
    //         // randomly init lattice
    //         // skip nskip (~1000 * n_sites) time steps
    //         // measure and store energy of system
    //     // take the mean and variance of this vector
    //     // (this is the avg energy and avg energy fluctuations)
    //     //
    //     // Return the avg energy, the fluctuations, and the uncertainty
    //     //
    //     // Q: what is the uncertainty on the energy fluctuation?
    //     // This is some basic stats from experimental methods...
    //     // goes like 1 / sqrt n_samples
    // }
    // TODO: implement the following
    // (doc) Monte Carlo estimation for average Magnetic Susceptibility
    // (doc) Returns the estimate and uncertainty 1 sigma
    // fn sample_magnetic_fluctuations(&self) -> [f64;3] {
    //     // initiate magnetism vector
    //     // for 0..n_samples
    //         // randomly init lattice
    //         // skip some time steps to cool system
    //         // measure and store magetic field
    //     // take the mean and variance of this vector
    //     //
    //     // Return avg magnetization, fluctuations, uncertainty
    // }
    /// Monte Carlo sample of energy
    /// Returns a vec of energy samples, of length 
    /// params.n_runs * params.samples_per_run
    fn sample_energy(&mut self, params: &MonteCarloParams) -> Vec<Vec<f64>> {
        let mut energy = vec![vec![0.0; params.samples_per_run]; params.n_runs];
        for i in 0..params.n_runs {
            self.reset_spins();
            // Time evolve the system to cool (or heat) it
            self.update_n(params.flips_to_skip);
            for j in 0..params.samples_per_run {
                // Time evolve the system a bit
                self.update_n(params.flips_to_skip_between_samples);
                energy[i][j] = self.measure_energy();
            }
        }
        return energy;
    }

    /// Monte Carlo sample of energy in parallel
    /// Returns a vec of energy samples, of length 
    /// params.n_runs * params.samples_per_run
    fn sample_energy_parallel(&mut self, params: &MonteCarloParams) -> Vec<Vec<f64>> {
        // let mut energy = vec![vec![0.0; params.samples_per_run]; params.n_runs]; // all samples
        let mut energy = vec![];
        let mut fetch_handle = vec![];
        for _ in 0..params.n_runs {
            // Create a clone: inits new lattice with same input params
            let mut lattice_copy = self.clone();
            let flips_to_skip = params.flips_to_skip;
            let flips_to_skip_between_samples = params.flips_to_skip_between_samples;
            let samples_per_run = params.samples_per_run;
            // Time evolve the system to cool (or heat) it
            fetch_handle.push(thread::spawn(move || -> Vec<f64> {
                lattice_copy.update_n(flips_to_skip);
                let mut erg_samples = vec![];
                for _ in 0..samples_per_run {
                    // Time evolve the system a bit
                    lattice_copy.update_n(flips_to_skip_between_samples);
                    erg_samples.push(lattice_copy.measure_energy());
                }
                erg_samples
            }));
        }
        // fetch_handle.into_iter().map(|c| energy.push(c.join().unwrap()));

        for thread in fetch_handle.into_iter() {
            energy.push(thread.join().unwrap());
        }
        return energy;
    }



    /// Monte Carlo estimate of nearest neighbor correlations
    /// Returns a vec of mean samples, of length params.n_runs
    fn sample_neighbor_correlations(&mut self, params: &MonteCarloParams) -> Vec<Vec<f64>> {
        // initiate nearest-neigbour correlation vector (empty Vec)
        let mut nn_corr = vec![vec![0.0; params.samples_per_run]; params.n_runs]; 
        for i in 0..params.n_runs {
            self.reset_spins();
            // Time evolve the system to cool (or heat) it
            self.update_n(params.flips_to_skip);
            for j in 0..params.samples_per_run {
                // Time evolve the system a bit
                self.update_n(params.flips_to_skip_between_samples);
                nn_corr[i][j] = self.get_dot_spin_neighbours() as f64 / self.n_sites as f64 / 4.0;
                // dividing by 4.0 scales it between -1 and +1, since 4 neighbours
            }
        }
        nn_corr
    }

    /// Monte Carlo estimate of nearest neighbor correlations
    /// Returns a vec of mean samples, of length params.n_runs
    fn sample_neighbor_correlations_parallel(&mut self, params: &MonteCarloParams) -> Vec<Vec<f64>> {
        // initiate nearest neighbour correlation vector (empty Vec)
        let mut nn_corr = vec![];
        let mut fetch_handle = vec![];
        for _ in 0..params.n_runs {
            // Create a clone: inits new lattice with same input params
            let mut lattice_copy = self.clone();
            let flips_to_skip = params.flips_to_skip;
            let flips_to_skip_between_samples = params.flips_to_skip_between_samples;
            let samples_per_run = params.samples_per_run;
            // Time evolve system to cool (or heat) it
            fetch_handle.push(thread::spawn(move || -> Vec<f64> {
                lattice_copy.update_n(flips_to_skip);
                let mut nn_samples = vec![];
                for _ in 0..samples_per_run {
                    // Time evolve the system a bit
                    lattice_copy.update_n(flips_to_skip_between_samples);
                    nn_samples.push(lattice_copy.get_dot_spin_neighbours() as f64 / lattice_copy.n_sites as f64 / 4.0);
                }
                nn_samples
            }));
        }
        for thread in fetch_handle.into_iter() {
            nn_corr.push(thread.join().unwrap());
        }
        return nn_corr;
    }


    fn sample_magnetization_parallel(&mut self, params: &MonteCarloParams) -> Vec<Vec<f64>> {
        let mut fetch_handle = vec![];
        for _ in 0..params.n_runs {
            // Create a clone: inits a new lattice with same input params
            let mut lattice_copy = self.clone();
            let flips_to_skip = params.flips_to_skip;
            let flips_to_skip_between_samples = params.flips_to_skip_between_samples;
            let samples_per_run = params.samples_per_run;
            // Time evolve system to cool (or heat) it
            fetch_handle.push(thread::spawn(move || -> Vec<f64> {
                lattice_copy.update_n(flips_to_skip);
                let mut mag_samples = vec![];
                for _ in 0..samples_per_run {
                    // Time evolve the system a bit
                    lattice_copy.update_n(flips_to_skip_between_samples);
                    mag_samples.push(lattice_copy.get_dot_spin_neighbours() as f64 / lattice_copy.n_sites as f64 / 4.0);
                }
                mag_samples
            }));
        }
        let mut magnetization = vec![];
        for thread in fetch_handle.into_iter() {
            magnetization.push(thread.join().unwrap());
        }
        return magnetization;
    }

    /// Monte Carlo sample the magnetization
    /// Returns a vec of mean samples, of length params.n_runs
    fn sample_magnetization(&mut self, params: &MonteCarloParams) -> Vec<Vec<f64>> {
        let mut mag_samples = vec![vec![0.0; params.samples_per_run]; params.n_runs];
        // Instead, init new lattice object and use multithreding for parallel processing
        for i in 0..params.n_runs {
            self.reset_spins();
            // Time evolve the system to cool (or heat) it
            self.update_n(params.flips_to_skip);
            for j in 0..params.samples_per_run {
                // Time evolve the system a bit
                self.update_n(params.flips_to_skip_between_samples);
                mag_samples[i][j] = self.get_spin_mean();
            }
        }
        mag_samples
    }


    // TODO: implement the following
    // (doc) Monte Carlo estimation for spacial correlations after system is settled
    // (doc) Returns the estimate and uncertainty 1 sigma
    // fn sample_spacial_correlations(&self , params:MonteCarloParams) -> [f64;3] {
    //     // initiate nearest neighbo
    // }
    // TODO: implement the following
    // (doc) Monte Carlo estimation for temporal correlations after system is settled
    // (doc) Returns the estaimate and uncertainty 1 sigma
    // fn sample_temporal_correlations(&self) -> (f64 , f64);
    // TODO: implement the following
    // (doc) Monte Carlo estimagion of all metrics
    // (doc) Returns a Vec (or dict?) of all the metrics
    // fn sample_estimate_all_metrics(&self , params:MonteCarloParams) -> Vec;
}

#[cfg(test)]
mod test {
    use super::*;

    #[test]
    fn test_sample_energy() {
        let params = MonteCarloParams {
            n_runs: 5,
            flips_to_skip: 1_000, // 1500_000,
            samples_per_run: 10,
            flips_to_skip_between_samples: 100,
        };
        let beta: f64 = 2.4; // roughly critical temp
                             // Initiate a lattice
        let mut lattice = Lattice2d::new(
            [9, 9],
            UpdateRule::Metropolis,
            SpinType::SpinHalf,
            InitType::Random,
            1.0f64, // j interaction constant
            0.0f64, // h static field term
            beta,   // 1/TkB
        );
        let energy_samples: Vec<Vec<f64>> = lattice.sample_energy(&params);
        assert_eq!(energy_samples.len(), params.n_runs);
        assert_eq!(energy_samples[0].len(), params.samples_per_run);
    }

    #[test]
    fn test_sample_energy_parallel() {
        let params = MonteCarloParams {
            n_runs: 5,
            flips_to_skip: 1_000,
            samples_per_run: 10,
            flips_to_skip_between_samples: 100,
        };
        let beta: f64 = 2.4;
        let mut lattice = Lattice2d::new(
            [9,9],
            UpdateRule::Metropolis,
            SpinType::SpinHalf,
            InitType::Random,
            1.0f64, // j interaction constant
            0.0f64, // h static field term
            beta,   // 1/TkB
        );
        let energy_samples: Vec<Vec<f64>> = lattice.sample_energy_parallel(&params);
        assert_eq!(energy_samples.len(), params.n_runs);
        assert_eq!(energy_samples[0].len(), params.samples_per_run);
    }

    #[test]
    fn test_sample_neighbor_correlations() {
        let params = MonteCarloParams {
            n_runs: 5,
            flips_to_skip: 1_000, // 1500_000,
            samples_per_run: 10,
            flips_to_skip_between_samples: 100,
        };
        let beta: f64 = 2.4; // roughly critical temp
                             // Initiate a lattice
        let mut lattice = Lattice2d::new(
            [9, 9],
            UpdateRule::Metropolis,
            SpinType::SpinHalf,
            InitType::Random,
            1.0f64, // j interaction constant
            0.1f64, // h static field term
            beta,   // 1/TkB
        );
        let nn_corr: Vec<Vec<f64>> = lattice.sample_neighbor_correlations(&params);
        assert_eq!(nn_corr.len(), params.n_runs);
        assert_eq!(nn_corr[0].len(), params.samples_per_run);
    }

    #[test]
    fn test_neighbor_correlations_parallel() {
        let params = MonteCarloParams {
            n_runs: 5,
            flips_to_skip: 1_000,
            samples_per_run: 10,
            flips_to_skip_between_samples: 100,
        };
        let beta: f64 = 2.4;
        let mut lattice = Lattice2d::new(
            [9,9],
            UpdateRule::Metropolis,
            SpinType::SpinHalf,
            InitType::Random,
            1.0f64, // j interaction constant
            0.0f64, // h static field term
            beta,   // 1/TkB
        );
        let nn_corr: Vec<Vec<f64>> = lattice.sample_neighbor_correlations_parallel(&params);
        assert_eq!(nn_corr.len(), params.n_runs);
        assert_eq!(nn_corr[0].len(), params.samples_per_run);
    }
}


// DISCUSSION
// Thought: For temporal correlations, the metric to measure this would be
// mixing time or something else?
//
// (1) The mixing time could be given as an estimate of the time that it takes
// for the probability that a site changes spin to be 1/2 within that window;
// kind of like a half-life. I wonder if there are ways to do this analytically,
// i.e. if you can derive this value from the correlations with immediate neighbors.
//
// (2) However we must be careful because this is not the only metric, if instead
// we measured the half life to be the time it takes for your spin to change
// *at least* once, this would give us a much longer-time estimate; also that
// metric makes more sense to me intuitively.
//
// (3) A third possible metric could meausre "how much time does it take for the
// probability of my spin being different from where it started to be 1/4"
// (will this exponentially decay towards 1/2, certainly that's the assymtote)
//
// I think (3) is the most legit option.
//
// (3) could be implemented in multiple ways,
// (3.1) a computationally expensive way
// could be to measure the correlation (square of diff) of the state at time
// t(0) with all subsequent times until t(cnst * n), where n is the number of
// sites. (perhaps this constant cnst could also depend on temperature);
// then curve fit an exponential decay and voila
//
// (3.2) a cheaper way could be to measure only the correlation between the
// t(0)'th sample and the t(cnst * n * f(temp))'s sample; I think this is the way
// to go...
//
// YOU NEED TO UNDERSTAND HOW CORRELATIONS IN SPACE AND TIME ARE MEASURED
// corr(x,y) = cov(x,y) / sqrt( var(x) var(y) )
//
// We should divide the correlation across time by the size of the lattice
// because one flip is attempted with every tick i.e. unit of time
// so larger lattices and graphs will have stronger correlations across time
// unless we divide by the number of sites.
//
// Should the correlations accuracy of correlations in time be handled by the user?
// Or should it be automatic. I think giving the user too much optionality
// might make this tool confusing to use...


