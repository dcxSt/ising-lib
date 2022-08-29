//! The 2D Spin Lattice Type. This is the first (and, so far, the only)
//! type of network (/graph) on which Monte-Carlo sampling can be done. 
//! You can use it with the built-in measurement and monte-carlo types
//! (refer to examples), or you can just use the lattice. 

use ndarray::prelude::*;
use rand::prelude::SliceRandom;
use rand::Rng;

/// Update rule options for Lattice 2d (Metropolis only for now)
#[derive(Clone, Copy)]
pub enum UpdateRule {
    Metropolis,
    Glauber,
}

// TODO: implement a hamiltonian type, for different models?
// consider renaming SpinType to Model Type, then we can call them
// things like SpinHalfFerromagnet or Sznajd or XY
// - SpinThreeHalf
// - XY
/// Types of spin system (SpinHalf only for now)
#[derive(Clone, Copy)]
pub enum SpinType {
    SpinHalf,
}

/// Initial condition (Random or AllUp)
#[derive(Clone, Copy)]
pub enum InitType {
    Random,
    AllUp,
}

/// A type encapsulating the 2d spin lattice
/// and basic operations performed on it
///
/// The lattice behaves like a torus - spins on
/// opposite edges are considered each other's
/// neighbours
///
/// The 2D lattice type
pub struct Lattice2d {
    pub dims: [usize; 2],
    pub n_sites: i32,       // the number of spin 1/2 sites == dims[0] * dims[1]
    pub nodes: Array2<i32>, // this language generalizes better to other graphs
    update_rule: UpdateRule,
    pub spin_type: SpinType,
    pub init_type: InitType,
    pub j: f64,    // interaction constant, default 1.0
    pub h: f64,    // external uniform magnetic field, default 0.0
    pub beta: f64, // beta = 1/(k_b * T), defaults to 0.43
}

impl Clone for Lattice2d {
    fn clone(&self) -> Self {
        Self::new(
            self.dims,
            self.update_rule,
            self.spin_type,
            self.init_type,
            self.j,
            self.h,
            self.beta,
        )
    }
}

/// Implement basic methods for the 2d lattice type
impl Lattice2d {
    /// Create a new lattice of given dims with specific implementation details
    pub fn new(
        dims: [usize; 2],
        update_rule: UpdateRule,
        spin_type: SpinType,
        init_type: InitType,
        j: f64,
        h: f64,
        beta: f64,
    ) -> Self {
        let nodes: Array2<i32> = Lattice2d::init_spins(&init_type, &dims);

        let (width, height) = nodes.dim();

        Lattice2d {
            dims: [width, height],
            n_sites: width as i32 * height as i32,
            nodes: nodes, // should it be called notes or sites?
            update_rule: update_rule,
            spin_type: spin_type,
            init_type: init_type,
            j: j,
            h: h,
            beta: beta,
        }
    }

    /// Create a new lattice of given dims with randomly generated spins
    pub fn new_basic(dims: [usize; 2]) -> Self {
        Self::new(
            dims,
            UpdateRule::Metropolis,
            SpinType::SpinHalf,
            InitType::Random,
            1.0f64,
            0.0f64,
            0.43f64,
        )
    }

    /// initiates the sites to some config (often random) as specified by init_type
    fn init_spins(init_type: &InitType, dims: &[usize; 2]) -> Array2<i32> {
        let nodes: Array2<i32>;
        match init_type {
            InitType::Random => {
                let mut rng = rand::thread_rng();
                nodes = Array2::from_shape_fn(*dims, |_| *[-1, 1].choose(&mut rng).unwrap());
            }
            InitType::AllUp => {
                nodes = Array2::<i32>::ones(*dims);
            }
        }
        return nodes;
    }

    /// resets the sites to some config (often random) as specified by init_type
    pub fn reset_spins(&mut self) {
        self.nodes = Lattice2d::init_spins(&self.init_type, &self.dims);
    }

    /// Gets the difference in energy from flipping the spin at [idx0,idx1]
    #[allow(non_snake_case)] // just for this function
    fn get_dE(&self, idx0: usize, idx1: usize) -> f64 {
        let neighbour_spin_sum: i32 = self.nodes[[idx0, (idx1 + 1) % self.dims[1]]]
            + self.nodes[[
                idx0,
                match idx1 {
                    0 => self.dims[1] - 1,
                    _ => idx1 - 1,
                },
            ]]
            + self.nodes[[(idx0 + 1) % self.dims[0], idx1]]
            + self.nodes[[
                match idx0 {
                    0 => self.dims[0] - 1,
                    _ => idx0 - 1,
                },
                idx1,
            ]];

        // two times dot prod of spin w/ it's neighbours
        // this is the energy required to flip
        2.0 * self.j * ((neighbour_spin_sum * self.nodes[[idx0, idx1]]) as f64)
            + self.h * (self.nodes[[idx0, idx1]] as f64)
    }

    /// Update the lattice by one timestep, (=one potential flip)
    pub fn update(&mut self) {
        match self.update_rule {
            UpdateRule::Metropolis => {
                // pick a random index
                let mut rng = rand::thread_rng();
                let idx0: usize = rng.gen::<usize>() % self.dims[0];
                let idx1: usize = rng.gen::<usize>() % self.dims[1];
                // determine weather to flip or not to flip
                #[allow(non_snake_case)]
                let dE: f64 = self.get_dE(idx0, idx1);
                if dE > 0.0 {
                    let p: f64 = rng.gen::<f64>(); // random f64 between 0 and 1
                    if p < (-self.beta * dE).exp() {
                        self.nodes[[idx0, idx1]] *= -1;
                    }
                } else {
                    self.nodes[[idx0, idx1]] *= -1; // something more complicated for spin 3/2
                }
            }
            UpdateRule::Glauber => {
                // not yet implemented
                panic!("Warning Glauber Rule not yet implemented");
            }
        }
    }

    /// Update the lattice by n timesteps
    pub fn update_n(&mut self, n: usize) {
        for _ in 0..n {
            self.update();
        }
    }

    /// Display lattice in terminal
    pub fn disp_terminal(&self) {
        let mut string = "----------------\n".to_owned();
        for idx0 in 0..self.dims[0] {
            string = string + "|";
            for idx1 in 0..self.dims[1] {
                match self.nodes[[idx0, idx1]] {
                    -1 => {
                        string = string + " ";
                    }
                    1 => {
                        string = string + "#";
                    }
                    _ => {
                        panic!("Ising lattice is an array of -1s and 1s");
                    }
                }
            }
            string = string + "|\n";
        }
        string = string + "---------------------";
        println!("{}", string);
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_lattice_new_basic() {
        let lattice = Lattice2d::new_basic([5, 10]);
        println!("{}", lattice.nodes); // cargo test -- --nocapture
        assert_eq!(lattice.nodes.raw_dim(), lattice.dims);
        assert_eq!(
            lattice.n_sites,
            lattice.dims[0] as i32 * lattice.dims[1] as i32
        )
    }

    #[test]
    fn test_lattice_new() {
        let lattice = Lattice2d::new(
            [7, 8],
            UpdateRule::Metropolis,
            SpinType::SpinHalf,
            InitType::Random,
            1.0,  // interaction constant, default 1.0
            0.0,  // external uniform magnetic field, default 0.0
            0.43, // beta = 1/(k_b * T), defaults to 0.43
        );
        assert_eq!(lattice.nodes.raw_dim(), lattice.dims);
        assert_eq!(
            lattice.n_sites,
            lattice.dims[0] as i32 * lattice.dims[1] as i32
        )
    }

    #[test]
    fn test_clone_lattice() {
        let lattice = Lattice2d::new_basic([3, 3]);
        let _newlat = lattice.clone();
    }

    #[test]
    fn test_update_n() {
        let mut lattice = Lattice2d::new_basic([5, 10]);
        lattice.update_n(10);
    }

    #[test]
    #[allow(non_snake_case)]
    fn test_get_dE() {
        let lattice = Lattice2d::new_basic([5, 6]);
        let i0: usize = 0;
        let i1: usize = 3;
        let _dE: f64 = lattice.get_dE(i0, i1);
    }

    #[test]
    fn test_init_spins() {
        let nodes: Array2<i32> =
            Lattice2d::init_spins(&InitType::Random, &[4 as usize, 5 as usize]);
        let (width, height) = nodes.dim();
        assert_eq!(width, 4usize);
        assert_eq!(height, 5usize);
        assert!(nodes[[3, 4]] == 1 || nodes[[3, 4]] == -1);
        assert!(nodes[[0, 0]] == 1 || nodes[[0, 0]] == -1);

        let nodes: Array2<i32> = Lattice2d::init_spins(&InitType::AllUp, &[2 as usize, 3 as usize]);
        let (width, height) = nodes.dim();
        assert_eq!(width, 2usize);
        assert_eq!(height, 3usize);
        assert_eq!(nodes[[1, 1]], 1i32);
    }

    #[test]
    fn test_reset_spins() {
        let mut lattice = Lattice2d::new_basic([5, 10]);
        lattice.reset_spins(); // all we test for here is runtime errors
    }

    #[test]
    fn test_update_disp() {
        let mut lattice = Lattice2d::new_basic([5, 5]);
        lattice.disp_terminal();
        // this tests the update function 300 times, it should only take an instant
        for _ in 0..3 {
            // set 3 to 3000 for slideshow
            for _ in 0..1 {
                // set 1 to 100 for slideshow
                lattice.update();
            }
            lattice.disp_terminal();
        }
    }
}




