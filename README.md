# QADsimulation
The dimension of a quantum state is traditionally seen as the number of superposed distinguishable states in a given basis. We propose an absolute, i.e.~basis-independent, notion of dimensionality for ensembles of quantum states. It is based on whether a quantum ensemble can be simulated with states confined to arbitrary lower-dimensional subspaces and classical postprocessing. In order to determine the absolute dimension of quantum ensembles, we develop both analytical witness criteria and a semidefinite programming criterion based on the ensemble's information capacity. Furthermore, we construct explicit simulation models for arbitrary ensembles of pure quantum states subject to white noise, and in natural cases we prove their optimality. Also, efficient numerical methods are provided for simulating generic ensembles. Finally, we discuss the role of absolute dimensionality in high-dimensional quantum information processing [1].

In this repository, we provide the function needed to evaluate the SDP for the QAD simulation.

# List of the files

- QADsimulation: this function evaluates the maximum visibility for the QAD = r simulation of a quantum ensemble of states subject to white noise;
- test_function: this is the example used in the paper to compute the results of Table I.

# References

[1]: A. Bernal, G. Cobucci, M.J. Renner, A. Tavakoli (2024). Absolute dimensionality of quantum ensembles (https://arxiv.org/abs/2409.01752)


