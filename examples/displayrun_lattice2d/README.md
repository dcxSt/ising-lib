# Example dry-run of Ising Lattice 2d
Visually displays in the terminal what's going on under the hood of the simulation. To run it just `cd` into this directory and run 

`cargo run`

# Demo

From the exmple's root directory (the same directory this README is in), build the release version with the `cargo build --release` command, then run it with `./target/release/displayrun_lattice2d`



https://user-images.githubusercontent.com/21654151/156935828-114c918a-d309-42ed-81c6-7f76f75c0f62.mov

:warning: The config file `Cargo.toml` uses `ising_lib` as a local dependency. Once version `^1.0.0` of `ising_lib` is published to [crates.io](https://crates.io/crates/ising_lib) the following config options
```toml
[dependencies.ising_lib]
path = "../../"
```
should be changed to 
```toml
[dependencies]
ising_lib = "1.0.0"
```
