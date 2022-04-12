# Rust Ising Library

Ising Lib is a tool to aid researchers perform ising model simulations on lattices and other graphs. It aims to provide researchers and curious people with a broad range of tools optimized for parallel processing and gpus. 

The ising model is popular because it is very simple to describe mathematically, and has a broad range of applications. Originally it was concieved to model spontaneous magnitization of ferromagnets due to cooling, but became very popular within the physics community in the mid 1900s for modelling many types of phase changes. ([Click here to learn about universality classes](https://www.wikiwand.com/en/Phase_transition#/Critical_exponents_and_universality_classes))

Thanks to increasing computational capabilities, today we are able to run larger and more Ising model simulations than ever before! Recently Ising models have started gaining traction within the social sciences in modelling herding behaviour and population dynamics. 


### Demo 

Visit [this website](https://ising-2d-lattice.netlify.app/) for an interactive demo. [Hint: try clicking on the `+T`, `-T`, and `random` buttons. When the system cools to below a critical temperature it 'quenches' or spontaneously magnetizes.]

https://user-images.githubusercontent.com/21654151/156935828-114c918a-d309-42ed-81c6-7f76f75c0f62.mov


### Implementation Details

This library makes it easy for you to estimate properties of a class of markov chains. They are called Ising models because they were concieved of by Dr. Ising in order to model ferromagnets. They became popular in the 60s (fact-check this) because it turns out they have many other applications (universality classes). 

**Types**

This library provides a type for each kind of graph you may want to use. Currently implemented are:
- **Lattice2D**, the typical spin-half (for now) lattice
- **Graph** (maybe rename to GraphGeneral)

*Wishlist: once we have a fully functional library that implements Lattice2D and Graph type for spin-half, we will first make the spin types more vercetile (e.g. with spin three-half or xy-model), then we will introduce more specialized graph types: Lattice1D (which is trivial to solve mathematically and will only be useful as an example), Lattice3D, LatticeND, and other types of graph.*

**Traits**

The **Measurement** trait measures quantities across each graph type. The quantities associated to this trait are those which can be measured instantaneously, such as the average spin or the energy of the lattice. (*link to docs here*)

The **MonteCarlo** trait probabilistically estimates quantities associated with the system considered as an ensemble by averaging across many runs. 

You can use these traits in the same way regardless of what the underlying graph structure is. I.e. with the same methods and associated functions. This way, once you see one example implementation across one type of graph, you've seen them all. 



### Citations

[1] [Wonseok Oh,  Sangyong Jeon, Membership Herding and Network Stability in the Open Source Community: The Ising Perspective, Management Science, Vol. 53, No. 7, 1086-1101 (2007)](https://www.jstor.org/stable/20122271). Draws from data from two open source projects. Concepts include: membership hearding, Open source communities. 

[2] [Ising Model of User Behavior Decision in Network Rumor Propagation](https://www.hindawi.com/journals/ddns/2018/5207475/). Discusses a use of ising models for rumour propagation. Uses 2d ising lattice for simulations. Interaction forces have three terms, the 'micropart' or 'self identity attribute', the 'middle part' or 'user-user interaction' (i.e. interaction constant `J`), the 'macroscopic part' or 'social enviornment's influence' (i.e. ext mag field `H`). Concepts include Von Neuman entropy, Game theory. 

[3] [The Ising Model: Brief Introduction and It's Application](https://www.intechopen.com/chapters/71210) is a short intro to Ising models, including some basic applications. 


---

Note: Ising lib is going to ship version 1.0.0 this Sunday. The maintainer of the `ising_lib` crate has kindly agreed to transfer ownership of his crate to my project. 

