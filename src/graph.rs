//! Graph spin 1/2 Type


use ndarray::prelude::*;
use rand::prelude::SliceRandom;
use rand::Rng;

/// update rule for Graph
pub enum UpdateRule {
    Metropolis,
}

/// types of spin system
pub enum SpinType {
    SpinHalf,
}

/// initial spin condition
pub enum InitType {
    Random, // uniformly random
    AllUp,
}

/// the type of edge between nodes
pub enum EdgeType {
    BinaryRandom { p:f64 }, // random edges {0 to 1} between nodes i,j with prob p
    UnifRandom, // uniformly random numbers between zero and 1
}

/// A type encapsulating an Ising model on
/// a graph and basic operations performed on it

pub struct Graph {
    pub n_sites: u32, // = nodes.len()
    pub nodes: Array<i32>, // An array of nodes
    pub edges: Array2<usize>, // 
    pub update_rule: UpdateRule,
    pub edge_type: EdgeType,
    pub spin_type: SpinType,
    pub init_type: InitType,
    pub j: f64, // interaction constant, default 1.0 
    pub h: f64, // external uniform magnetic field, default 0.0
    pub beta: f64, // beta = 1/(kb * T), defaults to 0.43
}

/// Implement basic methods on Graph type
impl Graph {
    /// Create a new Graph of given size with random edges
    pub fn new_basic(n_sites:u32, p:f64) -> Self {
        assert!(p>=0.0 && p<= 1.0);
        Self::new(
            n_sites,
            UpdateRule::Metropolis,
            EdgeType::BinaryRandom {p:p},
            SpinType::SpinHalf,
            InitType::Random,
            1.0f64,
            0.0f64,
            0.43f64,
        )
    }

    pub fn new(
        n_sites: u32,
        update_rule: UpdateRule,
        edge_type: EdgeType,
        spin_type: SpinType,
        init_type: InitType,
        j: f64,
        h: f64,
        beta: f64,
    ) -> Self {
        // TODO: implement init for different edge types
        // TODO: implement init for different spin types
        let nodes: Array<i32>;
        nodes = Array::from_shape_fn([n_sites], |_| *[-1,1].choose(&mut rng).unwrap());
        let mut edges: Array2<usize>;
        match edge_type {
            EdgeType::BinaryRandom { p } => {
                for i in 0..(n_sites.pow(2)) {
                    
                }
            }
        }
    }
}
