# Rust Ising Library

Ising Lib is a tool to aid researchers perform ising model simulations on lattices and other graphs. 

The ising model is popular because it is very simple to describe mathematically, and has a broad range of applications. Originally it was concieved to model spontaneous magnitization of ferromagnets due to cooling, but became very popular within the physics community in the mid 1900s for modelling many types of phase changes. ([Click here to learn about universality classes](https://www.wikiwand.com/en/Phase_transition#/Critical_exponents_and_universality_classes))

Thanks to increased computational capabilities, we are able to run larger and more Ising model simulations than ever before! Recently ising models have started gaining traction within the social sciences in modelling herding behaviour and population dynamics. 


### Demo 

Visit [this website](https://ising-2d-lattice.netlify.app/) for an interactive demo. [Hint: try clicking on the `+T`, `-T`, and `random` buttons. When the system cools to below a critical temperature it 'quenches' or spontaneously magnetizes.]

https://user-images.githubusercontent.com/21654151/156935828-114c918a-d309-42ed-81c6-7f76f75c0f62.mov


### Citations

[1] [Wonseok Oh,  Sangyong Jeon, Membership Herding and Network Stability in the Open Source Community: The Ising Perspective, Management Science, Vol. 53, No. 7, 1086-1101 (2007)](https://www.jstor.org/stable/20122271). Draws from data from two open source projects. Concepts include: membership hearding, Open source communities. 

[2] [Ising Model of User Behavior Decision in Network Rumor Propagation](https://www.hindawi.com/journals/ddns/2018/5207475/). Discusses a use of ising models for rumour propagation. Uses 2d ising lattice for simulations. Interaction forces have three terms, the 'micropart' or 'self identity attribute', the 'middle part' or 'user-user interaction' (i.e. interaction constant `J`), the 'macroscopic part' or 'social enviornment's influence' (i.e. ext mag field `H`). Concepts include Von Neuman entropy, Game theory. 

[3] [The Ising Model: Brief Introduction and It's Application](https://www.intechopen.com/chapters/71210) is a short survey of the history of ising models and some of it's basic applications. 


---

Note: Ising lib is going to ship version 1.0.0 this Sunday. The maintainer of the `ising_lib` crate has kindly agreed to transfer ownership of his crate to my project. 

