//! Everything you need to run Ising model simulations on different 
//! networks (although, so far, only the classic 2d grid ising model 
//! has been implemented as the Ising2d type). 
//! Despite its simplicity, the simulation allows us to observe an 
//! interesting physical phenomenon - phase transition.
//! Refer to the github repository for [examples](https://github.com/micouy/ising_lib). 


pub mod lattice2d;
pub mod measurement;
pub mod monte_carlo_measurement;
// pub mod prelude; // TODO: do this
// pub mod graph; // TODO: implement this

