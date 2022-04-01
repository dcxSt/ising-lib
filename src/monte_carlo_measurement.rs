use crate::lattice2d::*;
use crate::measurement::*;
use ndarray::prelude::*;

/// Parameters for monte carlo sampling
struct MonteCarloParams {
    n_per_sample: usize, // number of dry runs
    lattice_dims: [usize; 2],
    flips_to_skip: usize, // skip flips to allow for system to cool from random
    flips_per_measurement: usize,
}

/// The measurement trait measures quantities across different graphs
trait MonteCarlo {
    /// Calculates and returns basic metrics by monto carlo sampling
    /// - Energy fluctuations
    /// - Avg Magnetic Susceptibility
    // fn sample_metrics(&self) -> Array? Vector?, simple better probably immutable basic array [] is best; 
    /// Monte Carlo estimation for average Energy Fluctuations
    /// Returns tuple with estimate and uncertainty 1 sigma 
    fn sample_energy_fluctuations(&self) -> (f64 , f64);
    /// Monte Carlo estimation for average Magnetic Susceptibility
    /// Returns the estimate and uncertainty 1 sigma
    fn sample_magnetic_fluctuations(&self) -> (f64 , f64); 
    /// Monte Carlo estimation for spacial correlations after system is settled
    /// Returns the estimate and uncertainty 1 sigma
    fn sample_spacial_correlations(&self) -> (f64 , f64);
    // idea: temporal corrleations, need to see if this matters and how people usually implement this
    fn sample_temporal_correlations(&self) -> (f64 , f64);
    /// Monte Carlo estimation of all of the above
    /// re-implements each of the above metrics, and returns them in a vector
    // perhaps a dictionary would be better? what does rust offer instead of 
    // dictinoaries?
    fn sample_estimate_all_metrics(&self , params:MonteCarloParams) -> Vec;
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
    fn sample_energy_fluctuations(&self , params:MonteCarloParams) -> [f64;3] {
        // initiate energy vector
        // for 0..n_samples
            // randomly init lattice
            // skip nskip (~1000 * n_sites) timme steps
            // measure and store energy of system 
        // take the mean and variance of this vector
        // (this is the avg energy and avg energy fluctuations)
        //
        // Return the avg energy, the fluctuations, and the uncertainty
        //
        // Q: what is the uncertainty on the energy fluctuation?
        // This is some basic stats from experimental methods... 
        // goes like 1 / sqrt n_samples
    }
    /// Monte Carlo estimation for average Magnetic Susceptibility
    /// Returns the estimate and uncertainty 1 sigma
    fn sample_magnetic_fluctuations(&self) -> [f64;3] {
        // initiate magnetism vector
        // for 0..n_samples
            // randomly init lattice
            // skip some time steps to cool system
            // measure and store magetic field
        // take the mean and variance of this vector
        //
        // Return avg magnetization, fluctuations, uncertainty
    }
    /// Monte Carlo estimate of nearest neighbor correlations
    /// Returns the estimate and uncertainty 1 sigma
    fn sample_neighbor_correlations(&self , params:MonteCarloParams) -> [f64;2] {
        // initiate nn correlation vector (empty Vec)
        // for 0..n_samples
            // randomly init lattice
            // skip some time steps to cool system
            // measure the nearest neighbor correlation (same as interaction term)
    }
    /// Monte Carlo estimation for spacial correlations after system is settled
    /// Returns the estimate and uncertainty 1 sigma
    fn sample_spacial_correlations(&self , params:MonteCarloParams) -> [f64;3] {
        // initiate nearest neighbo
    }
    /// Monte Carlo estimation for temporal correlations after system is settled
    /// Returns the estaimate and uncertainty 1 sigma
    fn sample_temporal_correlations(&self) -> (f64 , f64);
    /// Monte Carlo estimagion of all metrics
    /// Returns a Vec (or dict?) of all the metrics
    fn sample_estimate_all_metrics(&self , params:MonteCarloParams) -> Vec;
}

#[cfg(test)]
mod test {
    use super::*;
}



