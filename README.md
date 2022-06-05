# Rust Ising Library

Ising Lib is a tool to aid researchers and curious people perform Ising model simulations on lattices and other graphs. It aims to provide a broad range of fast Ising simulation tools.

The ising model is popular because it is very simple to describe mathematically, and has a broad range of applications. Originally it was concieved to model spontaneous magnitization of ferromagnets due to cooling, but became very popular within the physics community in the mid 1900s for modelling many types of phase changes. ([Click here to learn about universality classes](https://www.wikiwand.com/en/Phase_transition#/Critical_exponents_and_universality_classes))

Thanks to increasing computational capabilities, today we are able to run larger and more Ising model simulations than ever before! Recently Ising models have started gaining traction within the social sciences in modelling herding behaviour and population dynamics. 


### Demo 

Visit [this website](https://ising-2d-lattice.netlify.app/) for an interactive demo. [Hint: try clicking on the `+T`, `-T`, and `random` buttons. When the system cools to below a critical temperature it 'quenches' or spontaneously magnetizes—this is an example of a phase change.]

Below is a toy in-terminal visualization demo ([examples/displayrun_lattice2d](https://github.com/dcxSt/ising-lib/tree/main/examples/displayrun_lattice2d)) to give you a look at what's going on under the hood of the simulation.

https://user-images.githubusercontent.com/21654151/156935828-114c918a-d309-42ed-81c6-7f76f75c0f62.mov

Below is a plot of the nearest-neighbour correlations with temperature. The transparent-green dots are data-points (i.e. samples), the blue crosses is the mean (i.e. the best guess estimate) and the red errorbar is 1 sigma of uncertainty in the mean. The below plot was generated for a 25x25 lattice grid and took about 3 minutes to run on my mac air (M1) (no gpu).

![Plot of Nearest Neighbor correlation against Temperature](https://github.com/dcxSt/ising-lib/blob/main/examples/correlation_with_nn/data/plot_76temps_nn_corr.png?raw=true)



### Implementation Details

This library makes it easy for you to estimate properties of a class of markov chains. They are called Ising models because they were concieved of by Dr. Ising in order to model ferromagnets. They became popular in the 60s (fact-check this) because it turns out they have many other applications (universality classes). 

**Types**

This library provides a type for each kind of graph you may want to use. Currently implemented are:
- **Lattice2D**, the typical spin-half (for now) lattice
TODO:
- **Lattice3D**
- **Graph** (maybe rename to GraphGeneral)

*Wishlist: once we have a fully functional library that implements Lattice2D and Graph type for spin-half, we will first make the spin types more vercetile (e.g. with spin three-half or xy-model), then we will introduce more specialized graph types: Lattice1D (which is trivial to solve mathematically and will only be useful as an example), Lattice3D, LatticeND, and other types of graph.*

**Traits**

The **Measurement** trait measures quantities across each graph type. The quantities associated to this trait are those which can be measured instantaneously, such as the average spin or the energy of the lattice. (*link to docs here*)

The **MonteCarlo** trait probabilistically estimates quantities associated with the system considered as an ensemble by averaging across many runs. 

You can use these traits in the same way regardless of what the underlying graph structure is. I.e. with the same methods and associated functions. This way, once you see one example implementation across one type of graph, you've seen them all. 

### TODO
- [x] Implement threading in MonteCarlo so that everything can run in [parallel](https://www.programming-idioms.org/cheatsheet/Rust)
  - [x] Implement deep clone for Lattice2d 
- [ ] Implement MonteCarlo trait, three metrics: energy, neighbor correlations, magnetization, in parallel. 
- [ ] Implement benchmarks for lattice 2d, including monte carlo, parallel processing etc.
- [ ] Implement [Sznajd model](https://www.wikiwand.com/en/Sznajd_model) hamiltonian for lattice 2d.
- [ ] Complete MonteCarlo trait for lattice2d
- [ ] Generate docs, make them pretty and informative
- [ ] Ship the lib
- [x] clean up `_convolve_2d_circ_neighbours` in measurements
- [ ] Implement 3d lattice
- [ ] Implement random graph and general graph type.
- [ ] Example calculate magnetic susceptibility
- [ ] Example calculate energy fluctuation

### Ideas
Monte Carlo Routines
- Correlations over larger spatial distances
- Correlations with n'th neighbour (n to the right)
- Correlations with (k,n)'th neighbour (n right, k up)
- Restructure so that you can sample multiple metrics each time
- Temporal correlations: put some though into how 'time' will scale. If you want to evaluate how much time is going by and compare different size grids, we need to scale the number of times we attempt a flip by nsize (number of sites).

### Citations

[1] [Wonseok Oh,  Sangyong Jeon, Membership Herding and Network Stability in the Open Source Community: The Ising Perspective, Management Science, Vol. 53, No. 7, 1086-1101 (2007)](https://www.jstor.org/stable/20122271). Draws from data from two open source projects. Concepts include: membership hearding, Open source communities. 

[2] [Ising Model of User Behavior Decision in Network Rumor Propagation](https://www.hindawi.com/journals/ddns/2018/5207475/). Discusses a use of ising models for rumour propagation. Uses 2d ising lattice for simulations. Interaction forces have three terms, the 'micropart' or 'self identity attribute', the 'middle part' or 'user-user interaction' (i.e. interaction constant `J`), the 'macroscopic part' or 'social enviornment's influence' (i.e. ext mag field `H`). Concepts include Von Neuman entropy, Game theory. 

[3] [The Ising Model: Brief Introduction and It's Application](https://www.intechopen.com/chapters/71210) is a short intro to Ising models, including some basic applications. 

[4] [Hadrien ising spins](https://github.com/HadrienG2/ising-spins)

---

Note: Ising lib is going to ship version 1.0.0 this Sunday. The maintainer of the `ising_lib` crate has kindly agreed to transfer ownership of his crate to my project. 

