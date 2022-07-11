# WeightedEnsemble
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
estimator of our QoI.

To guide the selection process, the state space is partitioned up into a
collection of bins.

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
      parent has exactly one offspring.
* In the new version, what were formerly interchangeable terms (`allocation` and
  `selection`) are now distinct concepts.  `allocation` involves deciding how to
  allocate a given number of particles amongst the bins, but does not actually
  resample the particles.  After it has been decided how to allocate the
  particle budget amongst all the bins, call `offspring_allocation!`, which
  decides how to allocate the particle budget within the bins.  Lastly, `repopulate!` is the function which will actually copy particles over.  

Parallel versions of the construction of the coarse transition matrix and the
actual WE are also included.  These distribute the work of the mutation step,
usually the most costly, and an inherently independent computation, across
workers. 

More details on the algorithm are described in ``WE_description.pdf``, which may
be found in the `doc` folder.

This package can be added to a Julia environment with the command:
```
(@v1.XYZ) pkg> add https://github.com/gideonsimpson/WeightedEnsemble.jl

```
# Examples

* The double well potential in dimension one, where we estimate the probability of
moving from the point -1 in the left basin to the right basin in a specified
amount of time.  The underlying dynamic in this version is the  discretized
overdamped Langevin equation.

* The Muller potential, where we estimate the probability of moving from one
minimum to another in a specified time.  This version uses the overdamped
Langevin equation.

* The Schlögl reaction rate network.  This version has parameters chosen so as
  to be bistable with resevoirs for species A and B.  Parameters were taken from
  Wang and Plecháč (2017).

# Notes

* The code is currently implemented for problems where the underlying problem is
  time homogeneous.  Thus, the `mutation` function should only take the current
  state as its argument.  However, both the `selection!` and `rebin!` functions
  will take the algorithmic time as an argument.  This is relevant for computing
  certain QoI and for adaptive binning strategies.
* When using `Distributed` parallelism, the parallelization is in the use of
  `pmap` for the mutation step.  A long term goal is to support the ensemble and
  bins datastructures as `SharedArrays` or `DistributedArrays`.
* Multithreading is now implemented using `trun_we` type commands.

# Caveats
* Make sure the number of threads/processes is set correctly before running the multiprocessing/multithreading examples
* The mutation step in the provided examples make use of [`BasicMD`](https://github.com/gideonsimpson/BasicMD.jl)
* One of the examples makes use of [`TestLandscapes`](https://github.com/gideonsimpson/TestLandscapes.jl)
* The mutation step in the Schlögl example code makes use of [`DiffEqBiological`](https://github.com/SciML/DiffEqBiological.jl)

# TO DO

* Modify code to handle and provide examples of steady state (reaction rate) problems
* Provide additional examples
* Update examples

# Acknowledgements
This work was supported in part by the US National Science Foundation Grants DMS-1818716 and DMS-2111278.

# References

1. _Analysis and optimization of weighted ensemble sampling_, D. Aristoff, ESAIM: M2AN, 52(4), 1219 - 1238, 2018. [(Link)](https://www.esaim-m2an.org/articles/m2an/abs/2018/04/m2an160145/m2an160145.html)
2. _Optimizing weighted ensemble sampling of steady states_, D. Aristoff and D.M. Zuckerman, arXiv:1806.00860. [(Link)](https://arxiv.org/abs/1806.00860)
3. _An ergodic theorem for weighted ensemble_, D. Aristoff, arXiv:1906.00856. [(Link)](https://arxiv.org/abs/1906.00856)
4. _Parallel replica dynamics method for bistable stochastic reaction networks: Simulation and sensitivity analysis_, T. Wang and P. Plecháč, J. Chem. Phys., 147, 234110, 2017. [(Link)](https://doi.org/10.1063/1.5017955)
5. _A splitting method to reduce MCMC variance_, R.J. Webber, D. Aristoff, G.  Simpson, arXiv:2011.13899. [(Link)](https://arxiv.org/abs/2011.13899)