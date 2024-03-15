# WeightedEnsemble

[![](https://img.shields.io/badge/docs-dev-blue.svg)](https://gideonsimpson.github.io/WeightedEnsemble.jl/dev)

Julia Implementation of the Weighted Ensemble Algorithm

# Overview

This Julia package provides tools for using the Weighted Ensemble algorithm for
estimating rare events.  The essential idea of this approach is to use an
ensemble of particles evolving in an unbiased way under the underlying dynamics
which are episodically selected, in the sense of a genetic algorithm, and
redistributed.  The total number of particles remains fixed (a user parameter),
but particles that are in the more important portions of state space _for the
quantity of interest (QoI)_ are resampled more frequently.  Additionally, each
particle has a weight, which is used to ensure that that we have an unbiased
estimator of our QoI.  To guide the selection process, the state space is
partitioned up into a collection of bins.

This package can be added to a Julia environment with the command:
```
(@v1.XYZ) pkg> add https://github.com/gideonsimpson/WeightedEnsemble.jl
```
We expect to add it to the general registry in the future.


The key functions and data structures that must be provided by the user are:
* An initial ensemble, stored in the `Ensemble` data structure, with initial
  particle positions and weights.
* A list of of bins, stored in the `Bins` data structure, together with a
  `rebin!` function that associates particle positions with bins.
* A `mutation` function which evolves particles under the underlying, unbiased,
  dynamic.
* A `selection!` function, which selects which particles to produce offspring
  before performing the `mutation`.  Two selection schemes are currently
  included:
    * `uniform_selection!` - This uniformly samples from the particles, ensuring
      that there is at least one particle spawned from a bin which is non-empty.
    * `targeted_selection!` - This allocates particles proprtionally to a user
      specified function, `G(p, E, B, t)`, for bin `p` at time `t`.
    * `optimal_selection!` - This uses a coarse scale model and a quantity of
      interest to allocate particles in a way that minimizes mutation variance.
      To use it, it is necessary to construct a "value function" which
      approximates the mutation variance of the problem.  The value function can
      be constructed by first building a coarse grained transition matrix for
      the model, form which "value_vectors", based on the QoI, can be
      constructed.  Tools for building up the value vectors and the coarse model
      are included with  `build_value_vectors` and
      `build_coarse_transition_matrix`.
    * `trivial_selection!` - This is included, for convenience, such that each
      parent has exactly one offspring.  It is useful for benchmarking against a direct simualtion with equivalent resources.

Parallel versions of the construction of the coarse transition matrix and the
actual WE are also included.  These distribute the work of the mutation step,
usually the most costly, and an inherently independent computation, across
workers. 

# Caveats and To Dos
* The code is currently implemented for problems where the underlying problem is
  time homogeneous.  Thus, the `mutation` function should only take the current
  state as its argument.  However, both the `selection!` and `rebin!` functions
  will take the algorithmic time as an argument.  This is relevant for computing
  certain QoI and for adaptive binning strategies.
* When using `Distributed` parallelism, the parallelization is in the use of
  `pmap` for the mutation step.  A long term goal is to support the ensemble and
  bins datastructures as `SharedArrays` or `DistributedArrays`.
* Multithreading is now implemented using `trun_we` type commands.
* Make sure the number of threads/processes is set correctly before running the multiprocessing/multithreading examples
* Modify code to handle and provide examples of steady state (reaction rate) problems
* Provide additional examples
* Update examples


# Collaborators
* D. Aristoff
* J. Copperman
* L. F. Doherty
* F. G. Jones
* R. J. Webber
* D. M. Zuckerman


# Acknowledgements
This work was supported in part by the US National Science Foundation Grants
DMS-1818716 and DMS-2111278.

# References
1. [_Analysis and optimization of weighted ensemble sampling_, D. Aristoff, ESAIM: M2AN, 52(4), 1219 - 1238, 2018](https://www.esaim-m2an.org/articles/m2an/abs/2018/04/m2an160145/m2an160145.html) 
2. [_Optimizing weighted ensemble sampling of steady states_, D. Aristoff and D.M. Zuckerman, MMS, 18(2), 2020.](https://epubs.siam.org/doi/10.1137/18M1212100)
3. [_An ergodic theorem for weighted ensemble_, D. Aristoff, J. Appl. Probab. 59, 152–166 (2022)](https://doi.org/10.1017/jpr.2021.38)
4. [_Parallel replica dynamics method for bistable stochastic reaction networks: Simulation and sensitivity analysis_, T. Wang and P. Plecháč, J. Chem. Phys., 147, 234110, 2017.](https://doi.org/10.1063/1.5017955)
5. [_A splitting method to reduce MCMC variance_, R.J. Webber, D. Aristoff, G.  Simpson, arXiv:2011.13899.](https://arxiv.org/abs/2011.13899)
6. [_Weighted ensemble: Recent mathematical developments_, D. Aristoff, J. Copperman, G. Simpson, R.J. Webber,  D.M. Zuckerman, J. Chem. Phys., 158, 014108, 2023](https://pubs.aip.org/aip/jcp/article/158/1/014108/2867485/Weighted-ensemble-Recent-mathematical-developments)
