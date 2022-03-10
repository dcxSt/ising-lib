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
    BinaryRandom { prob: f64 }, // random edges {0 to 1} between nodes i,j with prob p
    UnifRandom,                 // uniformly random numbers between zero and 1
}

/// A type encapsulating an Ising model on
/// a graph and basic operations performed on it

pub struct Graph {
    pub n_sites: u32,       // = nodes.len()
    pub nodes: Array1<i32>, // An array of nodes
    pub edges: Array2<f64>, // matrix
    pub update_rule: UpdateRule,
    pub edge_type: EdgeType,
    pub spin_type: SpinType,
    pub init_type: InitType,
    pub j: f64,    // interaction constant, default 1.0
    pub h: f64,    // external uniform magnetic field, default 0.0
    pub beta: f64, // beta = 1/(kb * T), defaults to 0.43
}

/// Implement basic methods on Graph type
impl Graph {
    /// Create a new Graph of given size with random edges
    pub fn new_basic(n_sites: u32, prob: f64) -> Self {
        assert!(prob >= 0.0 && prob <= 1.0);
        Self::new(
            n_sites,
            UpdateRule::Metropolis,
            EdgeType::BinaryRandom { prob: prob },
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
        // TODO: implement init for different spin types
        let nodes: Array1<i32>;
        let mut rng = rand::thread_rng();
        match init_type {
            InitType::Random => {
                nodes = Array::from_iter((0..n_sites).map(|_| *[-1, 1].choose(&mut rng).unwrap()));
            }
            InitType::AllUp => {
                nodes = Array::from_iter((0..n_sites).map(|_| 1));
            }
        }
        // TODO: implement init for different edge types
        let edges; // : Array2<f64> = Array2<f64>::zeros((n_sites, n_sites));
        match edge_type {
            EdgeType::BinaryRandom { prob } => {
                // Probabilistically fill the edge matrix with ones with prob p, and zeros with prob 1-p
                edges = Array2::from_shape_fn([n_sites as usize , n_sites as usize], |_| -> f64 {if rng.gen::<f64>() < prob {1.0} else {0.0}}); // *[0.0,1.0].choose(&mut rng, prob).unwrap());
            }
            _ => {
                panic!("Not yet implemented edge type, try EdgeType::BinaryRandom instead")
            }
        };
        Graph {
            n_sites:n_sites,
            nodes:nodes,
            edges:edges,
            update_rule:update_rule,
            edge_type:edge_type,
            spin_type:spin_type,
            init_type:init_type,
            j:j,
            h:h,
            beta:beta,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_graph_new() {
        let _graph: Graph = Graph::new(
            10u32,
            UpdateRule::Metropolis,
            EdgeType::BinaryRandom { prob: 0.2 },
            SpinType::SpinHalf,
            InitType::Random,
            1.0f64,
            0.0f64,
            0.4f64,
        );
    }
}
