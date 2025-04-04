# Binocular Rivalry Research

This repository contains the core implementation behind the primary task discussed in [New insights into binocular rivalry from the reconstruction of evolving percepts using model network dynamics](https://www.frontiersin.org/journals/computational-neuroscience/articles/10.3389/fncom.2023.1137015/full).

## Background

What is binocular rivalry? When two eyes perceive different stimuli, they alternate irregularly in a phenomenon termed binocular rivalry. 

The paper experimented with a new framework for characterizing the interplay between visual system connectivity, non-linear network dynamics, and stimulus encoding in binocular rivalry. We investigated binocular rivalry using a two-layer neuron network model where each layer competes based on distinct stimuli. We reconstructed the perceptions and found that the dominant image is prominently represented in the network. We observed slower perceptual shifts with specific parameters and found support for the theory that autism might stem from reduced brain inhibition.

### Note

The file [`main.cpp`](./main.cpp) includes the C++ code responsible for simulating rivalry between two networks with adaption. It serves as the foundation for the results and experiments presented in the paper.
