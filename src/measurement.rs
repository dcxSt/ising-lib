//! The Measurement type calculates quantities on systems that do not require 
//! averaging over multiple samples

use ndarray::Array2;
use crate::lattice2d::*;

/// The measurement trait measures quantities across different graphs
trait Measurement {
    fn get_spin_sum(&self) -> i32;              // get the sum of the spin values
    fn get_spin_expected_value(&self) -> f64;   // get the expected value of the spin, i.e. the magnetization per spin 
    fn _convolve_2d_circ_neighbours(mat:Array2<i32>) -> Array2<i32>; // convolves mat with filt with circular boundary conditions
    fn get_dot_spin_neighbours(&self) -> i32;   // get dot-product of each spin with the sum of it's neighbours
    fn measure_energy(&self) -> f64;            // get total energy of system
    fn measure_energy_per_spin(&self) -> f64;   // get energy divided by number of sites
    // IDEAS
    // spacial correlation function
        // correlation with immediate and 2nd degree neighbors... more?
        // is there a fast way to implement this?
}

/// Implement the measurement trait for the Lattice2d type
impl Measurement for Lattice2d {
    /// method returns sum of spins in lattice
    /// ∑ s_i
    fn get_spin_sum(&self) -> i32 {
        self.nodes.iter()
            .fold(0 , |acc, &x| acc + x)
    }
    /// method returns expected value of lattice
    /// ∑ s_i / n
    fn get_spin_expected_value(&self) -> f64 {
        self.get_spin_sum() as f64 / (self.n_sites as f64)
    }
    /// Convolves the 2d array mat, with a filter array filt
    /// Assumes periodic (/circular) boundary conditions
    /// This is kindof a bad hacky solution, we'll see how well it performs...
    fn _convolve_2d_circ_neighbours(mat:Array2<i32>) -> Array2<i32> {
        let mut result = Array2::<i32>::zeros(mat.shape());
        for i in 0..mat.shape()[0] {
            for j in 0..mat.raw_dim()[1] {
                match {
                    i == 0 && j == 0 => {
                        // result[[i,j]] = mat[]
                    }
                    i == 0 && j == mat.shape()[1] - 1 => {
                        // dummy
                    }
                    i == mat.shape()[0] - 1 && j == 0 => {
                        // dummy
                    }
                    i == mat.shape()[0] - 1 && j == mat.shape()[1] - 1 => {
                        // dummy 
                    }
                    _ => {
                        // dummy
                    }
                }
            }
        }
    }
    /// method returns dot of spins with their neighbors
    /// ∑ (s_i * s_j)   summing over all i,j pairs of neighbors
    fn get_dot_spin_neighbours(&self) -> i32 {
        // circular boudary convolution with neighbor filter
        // 0 1 0
        // 1 0 1
        // 0 1 0 
        // dot product of result with all_sites
        // (There may be room for optimization here... possibly a 2x speed 
        // up... at the expense of readable code.)
        0 // dummy
    }
    /// Return the energy of the lattice
    ///
    /// ```text
    /// E = -J * ∑(s_i * s_j) - H * ∑ s_i 
    /// ```
    fn measure_energy(&self) -> f64 {
        let spin_sum = self.get_spin_sum() as f64; // calculate H term
        let spin_neighbours_dot = self.get_dot_spin_neighbours() as f64; // J term
        // Q: should we take precautions in case of overflow errors here 
        // when converting from i32 to f64 ? 
        return - self.j * spin_neighbours_dot - self.h * spin_sum 
    }
    /// Returns the energy per spin
    fn measure_energy_per_spin(&self) -> f64 {
        self.measure_energy() / self.n_sites as f64 
    }
}

#[cfg(test)]
mod test {
    use super::*;

    #[test]
    fn test_spin_expected_value() {
        println!("\n----------------------\nTesting Lattice2d, spin_expected_value");
        let mut lattice = Lattice2d::new(
            [25,25],
            UpdateRule::Metropolis,
            SpinType::SpinHalf,
            InitType::Random,
            1.0,  // interaction constant, default 1.0
            0.0,  // external uniform magnetic field, default 0.0
            0.5, // beta = 1/(k_b * T), defaults to 0.43
        );
        // UNCOMMENT TWO LINES BELOW WITH --nocapture FLAG
        // use std::time::Duration;
        // std::thread::sleep(Duration::from_secs(3));
        for _ in 0..1 { // set 1 to 1000 for slideshow
            println!("< Spin > = {}", lattice.get_spin_expected_value());
            for _ in 0..1 { // set 1 to 100 for slideshow
                lattice.update();
            }
        }
    }
}



