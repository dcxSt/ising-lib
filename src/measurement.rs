use crate::lattice2d::*;

/// The measurement trait measures quantities across different graphs
trait Measurement {
    fn spin_expected_value(&self) -> f64; // get the expected value of the spin, i.e. the magnetization per spin 
    fn measure_energy(&self) -> f64; // energy of the lattice
}

// Implement the measurement trait for the Lattice2d type
impl Measurement for Lattice2d {
    fn spin_expected_value(&self) -> f64 {
        let mut sum = 0;
        for i in 0..self.dims[0] {
            for j in 0..self.dims[1] {
                sum += self.nodes[[i,j]];
            }
        }
        sum as f64 / ((self.dims[0] * self.dims[1]) as f64)
    }
    /// Return the energy of the lattice
    ///
    /// ```text
    /// E = -J * ∑(s_i * s_j) - H * ∑ s_i 
    /// ```
    fn measure_energy(&self) -> f64 {
        let spin_mean = self.spin_expected_value();
        // STOPPED CODING HERE
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
            println!("< Spin > = {}", lattice.spin_expected_value());
            for _ in 0..1 { // set 1 to 100 for slideshow
                lattice.update();
            }
        }
    }
}



