use crate::lattice2d::*;
use crate::measurement::Measurement; // *;

/// Parameters for monte carlo sampling
pub struct MonteCarloParams {
    pub n_runs: usize, // number of dry runs
    pub flips_to_skip: usize, // skip flips for system to cool 
    pub samples_per_run: usize, // number of samples to make in each run
    pub flips_to_skip_intra_run: usize, // number of flips to skip between each sample from the same run
}

/// The measurement trait measures quantities across different graphs
pub trait MonteCarlo {
    /// Calculates and returns basic metrics by monto carlo sampling
    /// - Energy fluctuations
    /// - Avg Magnetic Susceptibility
    // fn sample_metrics(&self) -> Array? Vector?, simple better probably immutable basic array [] is best; 
    // TODO: implement the following
    // (doc) Monte Carlo estimation for average Energy Fluctuations
    // (doc) Returns tuple with estimate and uncertainty 1 sigma 
    fn sample_energy(&mut self,params:&MonteCarloParams) -> Vec<f64>;
    fn sample_energy_mean_var(&mut self,params:MonteCarloParams) -> [f64;2];
    // fn sample_energy_fluctuations(&self , params:MonteCarloParams) -> (f64 , f64);
    // TODO: implement the following
    // (doc) Monte Carlo estimation for average Magnetic Susceptibility
    // (doc) Returns the estimate and uncertainty 1 sigma
    // fn sample_magnetic_fluctuations(&self) -> (f64 , f64); 
    // TODO: implement the following 
    // (doc) Monte Carlo estimation for spacial correlations after system is settled
    // (doc) Returns the estimate and uncertainty 1 sigma
    // fn sample_spacial_correlations(&self) -> (f64 , f64);
    fn sample_neighbor_correlations(&mut self,params:&MonteCarloParams) -> Vec<f64>; 
    fn sample_neighbor_correlations_mean_var(&mut self, params:MonteCarloParams) -> [f64;2];
    fn sample_average_magnetization(&mut self , params:MonteCarloParams) -> [f64;2];
    // idea: temporal corrleations, need to see if this matters and how people usually implement this
    // TODO: implement below function
    // fn sample_temporal_correlations(&self) -> (f64 , f64);
    // (doc) Monte Carlo estimation of all of the above
    // (doc) re-implements each of the above metrics, and returns them in a vector
    // perhaps a dictionary would be better? what does rust offer instead of 
    // dictinoaries?
    // TODO: implement below function
    // fn sample_estimate_all_metrics(&self , params:MonteCarloParams) -> Vec;
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
//
// We should divide the correlation across time by the size of the lattice
// because one flip is attempted with every tick i.e. unit of time 
// so larger lattices and graphs will have stronger correlations across time 
// unless we divide by the number of sites. 
//
// Should the correlations accuracy of correlations in time be handled by the user?
// Or should it be automatic. I think giving the user too much optionality 
// might make this tool confusing to use...


/// Implements the measurement trait for the Lattice2d type
impl MonteCarlo for Lattice2d {
    // fn sample_energy_fluctuations(&self , params:MonteCarloParams) -> [f64;3] {
    //     // initiate energy vector
    //     // for 0..n_samples
    //         // randomly init lattice
    //         // skip nskip (~1000 * n_sites) timme steps
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
    /// Returns a vec of energy samples, of length params.n_runs * params.samples_per_run
    fn sample_energy(&mut self,params:&MonteCarloParams) -> Vec<f64> {
        let mut energy = vec![0.0; params.n_runs * params.samples_per_run];
        for i in 0..params.n_runs {
            self.reset_spins();
            // Time evolve the system to cool (or heat) it
            for _ in 0..params.flips_to_skip { self.update(); }
            for j in 0..params.samples_per_run {
                // Time evolve the system a bit
                for _ in 0..params.flips_to_skip_intra_run { self.update(); }
                energy[i * params.samples_per_run + j] = self.measure_energy();
            }
        }
        return energy
    }
    /// Monte Carlo estimate of mean in energy and fluctuation in energy
    /// Returns the mean and variance of the energies sampled 
    fn sample_energy_mean_var(&mut self, params:MonteCarloParams) -> [f64;2] {
        let energy:Vec<f64> = self.sample_energy(&params);
        let mean_energy:f64 = energy.iter().sum::<f64>() / energy.len() as f64;
        let var_energy = energy.iter().map(|x| f64::powi(*x,2) - mean_energy.powi(2)).sum::<f64>() / ((params.n_runs * params.samples_per_run) as f64);
        return [mean_energy, var_energy]
    }
    /// Monte Carlo estimate of nearest neighbor correlations
    /// Returns a vec of mean samples, of length params.n_runs
    fn sample_neighbor_correlations(&mut self,params:&MonteCarloParams) -> Vec<f64> {
        assert!(params.n_runs >= 2); // There must be at least two runs, or we cannot make estimate of var in mean
        // initiate nn correlation vector (empty Vec)
        let mut nn_corr = vec![0.0; params.n_runs]; // nearest neighbor correlations

        for i in 0..params.n_runs {
            let mut nn_corr_run = vec![0.0; params.samples_per_run];
            self.reset_spins();
            // Time evolve the system to cool (or heat) it 
            for _ in 0..params.flips_to_skip { self.update(); }
            for j in 0..params.samples_per_run {
                // Time evolve the system a bit
                for _ in 0..params.flips_to_skip_intra_run { self.update(); }
                nn_corr_run[j] = self.get_dot_spin_neighbours() as f64 / self.n_sites as f64 / 4.0; 
                // dividing by 4.0 scales it between -1 and +1
            }
            nn_corr[i] = nn_corr_run.iter().sum::<f64>() / nn_corr_run.len() as f64;
        }
        return nn_corr
    }
    /// Monte Carlo estimate of nearest neighbor correlations
    /// Returns the mean and variance in the mean of the sample
    fn sample_neighbor_correlations_mean_var(&mut self, params:MonteCarloParams) -> [f64;2] {
        let nn_corr:Vec<f64> = self.sample_neighbor_correlations(&params);
        let mean_corr:f64 = nn_corr.iter().sum::<f64>() / nn_corr.len() as f64;
        let var:f64 = nn_corr.iter().map(|x| f64::powi(*x,2) - mean_corr.powi(2)).sum::<f64>() / params.n_runs as f64;
        let var_in_mean:f64 = var / (params.n_runs as f64 - 1.0);
        return [mean_corr , var_in_mean]
    }
    /// Monte Carlo estimate of avg magnetization
    /// Returns the estimate and uncertainty 1 var
    fn sample_average_magnetization(&mut self , params:MonteCarloParams) -> [f64;2] {
        assert!(params.n_runs >= 2); // There must be at least two runs, or we cannot make estimate of var in mean
        let mut abs_spin = vec![0.0; params.n_runs]; 
        for i in 0..params.n_runs {
            let mut abs_spin_run = vec![0.0; params.samples_per_run];
            self.reset_spins();
            // Time evolve the system to cool (or heat) it
            for _ in 0..params.flips_to_skip { self.update(); }
            for j in 0..params.samples_per_run {
                // Time evolve the system a bit
                for _ in 0..params.flips_to_skip_intra_run { self.update(); }
                abs_spin_run[j] = self.get_spin_expected_value().abs();
            }
            abs_spin[i] = abs_spin_run.iter().sum::<f64>() / params.samples_per_run as f64;
        }
        let mean_abs_spin:f64 = abs_spin.iter().sum::<f64>() / params.n_runs as f64;
        let var:f64 = abs_spin.iter().map(|x| f64::powi(*x,2) - mean_abs_spin.powi(2)).sum::<f64>() / params.n_runs as f64;
        let var_in_mean:f64 = var / (params.n_runs as f64 - 1.0);
        [mean_abs_spin, var_in_mean]
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
    fn test_sample_neighbor_correlations() {
        // TODO: implement this test
    }
}



