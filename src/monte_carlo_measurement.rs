use crate::lattice2d::*;
use crate::measurement::*;
use ndarray::prelude::*;

/// Parameters for monte carlo sampling
struct MonteCarloParams {
    n_samples: usize,
    beta_range: Array1::<f64>, // array of beta, must be above zero
    lattice_dims: [usize; 2],
    flips_to_skip: usize,
    measurements_per_beta: usize,
    flips_per_measurement: usize,
}

/// The measurement trait measures quantities across different graphs
trait MonteCarlo {
    /// Calculates and returns basic metrics by monto carlo sampling
    /// - Energy fluctuations
    /// - Avg Magnetic Susceptibility
    // fn sample_metrics(&self) -> Array? Vector?, simple better probably immutable basic array [] is best; 
    /// Monte Carlo estimation for average Energy Fluctuations
    fn sample_energy_fluctuations(&self) -> f64;
    /// Monte Carlo estimation for average Magnetic Susceptibility
    fn sample_magnetic_suceptibility(&self) -> f64;
}

// Implement the measurement trait for the Lattice2d type
impl MonteCarlo for Lattice2d {
    fn sample_energy_fluctuations(&self, params:MonteCarloParams) -> f64 {
        
    }
}

#[cfg(test)]
mod test {
    use super::*;
}



