# JuWeightedEnsemble
Julia Implementation of the Weighted Ensemble Algorithm

# Overview

This Julia package provides tools for using the Weighted Ensemble algorithm for estimating rare events.  The essential idea of this approach is to use an ensemble of particles evolving in an unbiased way under the underlying dynamics which are episodically selected, in the sense of a genetic algorithm, and redistributed.  The total number of particles remains fixed (a user parameter), but particles that are in the more important portions of state space _for the quantity of interest (QoI)_ are resampled more frequently.  Additionally, each particle has a weight, which is used to ensure that that we have an unbiased estimator of our QoI.

To guide the selection process, the state space is partitioned up into a collection of bins.

The key functions and data structures that must be provided by the user are:

* An initial ensemble, stored in the `Ensemble` data structure, with initial particle positions and weights.
* A list of of bins, stored in the `Bins` data structure, together with a `bin_id` function that uniquely maps positions in state to one of the bins.  
* A `mutation` function which evolves particles under the underlying, unbiased, dynamic
* Value vectors, `value_vectors`, and a resampling function (`Systematic`, `Residual`, and `Stratified` are included) to guide the selection process.

The value vectors can be constructed by building up a bin to bin coarse Markov model, estimating the QoI when restricted to the bins, and then using the included `build_value_vectors`.  A coarse Markov model can be constructed with the included `build_coarse_transition_matrix`

Parallel versions of the construction of the coarse transition matrix and the actual WE are also included.  These distribute the work of the mutation step, usually the most costly, and an inherently independent computation, across workers.  

More details on the algorithm are described in the references below.    In particular, this code uses optimal allocation of the particles to minimize selection variance.

# Examples

* The double well potential in dimension one, where we estimate the probability of moving from the point -1 in the left basin to the right basin in a specified amount of time.  The underlying dynamic in this version is the discretized overdamped Langevin equation.

* The Muller potential, where we estimate the probability of moving from one minimum to another in a specified time.  This version uses the overdamped Langevin equation.

# Caveats

* This version only includes optimal allocation
* The mutation step in the provided examples make use of [`JuBasicMD`](https://github.com/gideonsimpson/JuBasicMD)

# TO DO

* Support uniform allocation, if only for comparison
* Provide examples of steady state (reaction rate) problems
* Provide additional examples

# References

1. _Analysis and optimization of weighted ensemble sampling_, D. Aristoff, ESAIM: M2AN, 52(4), 1219 - 1238, 2018. [(Link)](https://www.esaim-m2an.org/articles/m2an/abs/2018/04/m2an160145/m2an160145.html)

2. _Optimizing weighted ensemble sampling of steady states_, D. Aristoff and D.M. Zuckerman, arXiv:1806.00860. [(Link)](https://arxiv.org/abs/1806.00860)
3. _An ergodic theorem for weighted ensemble_, D. Aristoff, arXiv:1906.00856. [(Link)](https://arxiv.org/abs/1906.00856)
