//! The Measurement type calculates quantities on networks that can be 
//! measured instantaneously, (such as the magnetization at a given 
//! point in time), i.e. quantities that do not require averaging over 
//! multiple samples.

use ndarray::prelude::*;
use crate::lattice2d::*;

/// The measurement trait measures quantities across different graphs.
/// Each method returns a 2d vector of dim (x,y) where x is the number
/// of simulated graphs run, and y is the number of samples taken with
/// each run. 
pub trait Measurement {
    fn get_spin_sum(&self) -> i32;      // get the sum of the spin values
    fn get_spin_mean(&self) -> f64;     // get the mean value of spins 
    fn _convolve_2d_circ_neighbours(mat:&Array2<i32>) -> Array2<i32>; // convolves mat with filt with circular boundary conditions
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

    /// method returns mean spin of lattice
    /// ∑ s_i / n
    fn get_spin_mean(&self) -> f64 {
        self.get_spin_sum() as f64 / (self.n_sites as f64)
    }

    /// Convolves the 2d array mat, with a filter array filt
    /// Assumes periodic (/circular) boundary conditions
    /// This can still be optimized
    fn _convolve_2d_circ_neighbours(mat:&Array2<i32>) -> Array2<i32> {
        // Fill the result matrix (result of convolution)
        let roll = |ix: usize, amt: i32, max: usize| {
            let max = max as i32;
            ((ix as i32 + amt + max) % max) as usize // +max, might be neg
        };
        let (width, height) = mat.dim();
        Array2::from_shape_fn((width, height), |ix| {
            mat[[ix.0,roll(ix.1,1,height)]] 
                + mat[[ix.0,roll(ix.1,-1,height)]]
                + mat[[roll(ix.0,1,width),ix.1]]
                + mat[[roll(ix.0,-1,width),ix.1]]
        })
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
        // up... at the expense of readable code?)
        let neighbors = Lattice2d::_convolve_2d_circ_neighbours(&self.nodes);
        assert_eq!(neighbors.shape() , self.nodes.shape());
        let mut dot_spin:i32 = 0;
        for (x,y) in self.nodes.iter().zip(neighbors) {
            dot_spin += x * y;
        }
        dot_spin 
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
    fn test_spin_mean() {
        println!("\n----------------------\nTesting Lattice2d, spin_expected_value");
        let mut lattice = Lattice2d::new(
            [25,25],
            UpdateRule::Metropolis,
            SpinType::SpinHalf,
            InitType::Random,
            1.0, // interaction constant, default 1.0
            0.0, // external uniform magnetic field, default 0.0
            0.5, // beta = 1/(k_b * T), defaults to 0.43
        );
        // UNCOMMENT TWO LINES BELOW WITH --nocapture FLAG
        // use std::time::Duration;
        // std::thread::sleep(Duration::from_secs(3));
        for _ in 0..1 { // set 1 to 1000 for slideshow
            println!("< Spin > = {}", lattice.get_spin_mean());
            for _ in 0..1 { // set 1 to 100 for slideshow
                lattice.update();
            }
        }
    }

    #[test]
    fn test_convolve_2d_circ_neighbours() {
        let vec1 = vec![
            vec![0,0,0],
            vec![0,1,0],
            vec![0,0,0],
        ];
        let vec1_conv = vec![ // we expect vec1 to convolve into this
            vec![0,1,0],
            vec![1,0,1],
            vec![0,1,0],
        ];
        let vec2 = vec![
            vec![0,0,0,1],
            vec![0,0,0,0],
            vec![0,0,0,0],
            vec![0,0,0,0],
        ];
        // we expect the convolution operator to turn vec2 into vec2_conv
        let vec2_conv = vec![
            vec![1,0,1,0],
            vec![0,0,0,1],
            vec![0,0,0,0],
            vec![0,0,0,1],
        ];
        let mut arr1 = Array2::<i32>::default((3,3));
        for (i, mut row) in arr1.axis_iter_mut(Axis(0)).enumerate() {
            for (j, col) in row.iter_mut().enumerate() {
                *col = vec1[i][j];
            }
        }
        let mut arr2 = Array2::<i32>::default((4,4));
        for (i, mut row) in arr2.axis_iter_mut(Axis(0)).enumerate() {
            for (j, col) in row.iter_mut().enumerate() {
                *col = vec2[i][j];
            }
        }

        let result1 = Lattice2d::_convolve_2d_circ_neighbours(&arr1);
        let result2 = Lattice2d::_convolve_2d_circ_neighbours(&arr2);

        for i in 0..3 {
            for j in 0..3 {
                assert_eq!(vec1[i][j] , arr1[[i,j]]); // just test intuition about array referencing/dereferencing
                assert_eq!(result1[[i,j]] , vec1_conv[i][j]);
                assert_eq!(result2[[i,j]] , vec2_conv[i][j]);
            }
        }
    }

    #[test]
    fn test_get_dot_spin_neighbours() {
        let lattice = Lattice2d::new([2,3],
            UpdateRule::Metropolis,
            SpinType::SpinHalf,
            InitType::AllUp,
            1.0, // interaction constant, default 1.0
            0.0, // external uniform magnetic field, default 0.0
            0.5, // beta = 1/(k_b * T), defaults to 0.43
        );

        // The spins have all been initialized to up
        assert_eq!(lattice.get_dot_spin_neighbours() , 6 * 4);
    }

}



