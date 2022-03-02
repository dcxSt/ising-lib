//! The 2D Spin Lattice Type

// use ndarray::{prelude::*, NdIndex};
use ndarray::prelude::*;
// use rand::prelude::*;
use rand::prelude::SliceRandom;
use rand::Rng;

/// update rule for Lattice 2d
pub enum UpdateRule {
    Metropolis,
    Glauber,
}

/// types of spin system
pub enum SpinType {
    SpinHalf,
    SpinThreeHalf,
    XY,
}

/// initial condition
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
    dims: [usize; 2],
    n_sites: i32,       // the number of spin 1/2 sites == dims[0] * dims[1]
    nodes: Array2<i32>, // this language generalizes better to other graphs
    update_rule: UpdateRule,
    spin_type: SpinType,
    init_type: InitType,
    j: f64,    // interaction constant, default 1.0
    h: f64,    // external uniform magnetic field, default 0.0
    beta: f64, // beta = 1/(k_b * T), defaults to 0.43
}

/// Implement basic methods for the 2d lattice type
impl Lattice2d {
    /// Create a new lattice of given dims with randomly generated spins.
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
        // TODO: implement initialization for different spin types
        let nodes: Array2<i32>;
        let mut rng = rand::thread_rng();
        nodes = Array2::from_shape_fn(dims, |_| *[-1, 1].choose(&mut rng).unwrap());

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

    /// Gets the difference in energy from flipping a spin
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
        // Calculation with H, not yet implemented
        2.0 * self.j * ((neighbour_spin_sum * self.nodes[[idx0, idx1]]) as f64)
    }

    /// Update the lattice by one timestep
    pub fn update(&mut self) {
        match self.update_rule {
            UpdateRule::Metropolis => {
                // pick a random index
                let mut rng = rand::thread_rng();
                let idx0: usize = rng.gen::<usize>() % self.dims[0];
                let idx1: usize = rng.gen::<usize>() % self.dims[1];
                // determine weather to flip or not to flip
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
                println!("Warning Glauber Rule not yet implemented");
            }
        }
    }

    /// Display lattice in terminal
    pub fn disp_terminal(&self) {
        println!("----------------");
        for idx0 in 0..self.dims[0] {
            print!("|");
            for idx1 in 0..self.dims[1] {
                match self.nodes[[idx0, idx1]] {
                    -1 => {
                        print!(" ");
                    }
                    1 => {
                        print!("*");
                    }
                    _ => {
                        panic!("Ising lattice is an array of -1s and 1s");
                    }
                }
            }
            print!("|\n");
        }
        println!("\n----------------");
    }
}

pub struct Graph {
    n_sites: i32,
    nodes: Array1<i32>,
    edges: Array2<i32>, // less memory efficient, more readable
    update_rule: UpdateRule,
    spin_type: SpinType,
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn playground() {
        assert!(true);
    }

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
    fn test_get_dE() {
        let lattice = Lattice2d::new_basic([5, 10]);
        let i0: usize = 0;
        let i1: usize = 9;
        let dE: f64 = lattice.get_dE(i0, i1);
    }

    #[test]
    fn test_update() {
        let mut lattice = Lattice2d::new_basic([25, 75]);
        lattice.disp_terminal();
        for i in 0..3000 {
            for j in 0..100 {
                lattice.update();
            }
            lattice.disp_terminal();
        }
    }
}
