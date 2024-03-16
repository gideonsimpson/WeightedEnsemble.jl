# WeightedEnsemble.jl

```@contents
Pages = ["index.md"]
Depth = 2
```


## Algorithm Overview
```@contents
Pages = ["index.md"]
Depth = 3:4
```

Weighted Ensemble (WE) is a variance reduction strategy that can be used to,
amongst other things, estimate rare events.  More generally, it is designed to estimate
```math
\mathbb{E}_{\mu}[f(X)] = \mu(f)
```
for some Quantity of Interest (QoI), ``f``, with respect to some target measure,
``\mu``.   One of the major advantages of WE  over other variance reduction
strategies is that it makes rather minimal assumptions on the Markov underlying
process.  Indeed, there is no assumption of reversibility of the process, nor is
there any need to evaluate or otherwise estimate a likelihood.  It is sufficient
to be able to simulate the associated process, i.e., sample a Markov kernel,
``K``, for which ``\mu`` is the stationary distribution.  A key application of
this approach is in the computation of mean first passage times (MFPTs) via the
[Hill relation](https://statisticalbiophysicsblog.org/?p=8).

This implementation of WE makes use of a fixed number, ``N``, particles, ``\xi^{i}``, each carrying weight ``\omega^{i}``.  At algorithmic time ``t``, we denote
```math
\mu_t^{N}(dx) = \sum_{i=1}^N \omega_t^{i}\delta_{\xi_t^{i}}(dx)
```
The weights are assumed to be positive and sum to one.  The particle ensemble is
stored in an [`Ensemble`](@ref) data structure. The ensemble is evolved,
successively, in two steps, _selection_ and _mutation_.  As the partilces and
weights evolve in time, we add a ``t`` subscript to indicate that these quantities change, ``\{(\omega_t^{i},\xi_t^{i}\}``.

This implementation has been used in [aristoff_weighted_2023, webber_splitting_2020](@cite).  It is heavily influenced by the algoirthmic
description found in [aristoff_optimizing_2020](@cite).  

### Selection
During the selection step, particles are resampled from the empirical
distribution in such a way so as to maintain unbiasedness _and_ reduce variance
with respect to some quantity of interest.  The ensemble after selection is
represented by
```math
\hat{\mu}_t^{N}(dx) = \sum_{i=1}^N \hat{\omega}_t^{i}\delta_{\hat{\xi}_t^{i}}(dx)
```
See [Selection Algorithms](@ref) and related documents for additional details on
this step.

### Bins
To understand more of how the selection step works, it is essential to introduce
the concept of bins.  Indeed, this implementation of WE is _bin based_, in the
following sense.  The indices of the particles are partitioned into bins, ``u_1, u_2,\ldots, u_{M}``, where ``M`` is the number of partitions.  Often, these bins are determined by a disjoint partition of the state space, ``\mathcal{X}`` as
```math
\mathcal{X} = \bigcup_{i=1}^M \mathcal{B}_i, \quad j \in u_i \leftrightarrow \xi_t^j\in \mathcal{B}_i
```
Though the presentation here assumes that the bins are statically (independent
of ``t``), these need not be the case in general.  Associated with each bin is the particle count at time ``t``, ``n_t^i=|u_i|``, and the bin weight, ``\nu_t^i``, given by
```math
\nu^i_t = \sum_{j\in u_i}\omega^j.
```
The bins and their features are recorded in a [`Bins`](@ref) data structure.

To perform the selection step, one first performs the allocation steps:
* First, we determine how many particles should reside in each of the bins, ``C_t(u)``, after selection.  This is implemented to satisfy total particle conservation
```math
\sum_{i=1}^{M}C_t(u_i) = N.
```
A simple, though suboptimal, choice is uniform allocation, in which ``C_t(u_i)\approx 1/M``.

* Second, we determine how many offspring each particle should have.  This is implemented to constrain the bin count from the first step,
```math
C_t(u_i)=\sum_{j\in u_i}C_t^j.
```
This is typically done using multinomial resampling, with
```math
(C_t^{j_k}) \sim \mathrm{Multi.}\left(C_t(u_i),\frac{\omega_t^{j_1}}{\nu_t^i},\frac{\omega_t^{j_2}}{\nu_t^i},\ldots \right)
```
where the ``j_k \in u_i``.  

Having allocated the number of offspring to each particle, we then repopulate our ensemble as:
* Offspring of particle ``\xi_{t}^j`` are assigned this as their position;
* Offspring of particles in bin ``u_i`` are all assigned weight ``\nu_t^i/C_t(u_i)``.

There are some additional subtlesties that are implemented to ensure the
allocation and selection with bins behaves as expected.

### Mutation
During the mutation step, the particles evolve freely over one algorithmic time unit according to a user defined Markov kernel, ``K``, with ``\xi_{t+1}^{i}\sim K(\hat{\xi}_t^{i}, dx)``, and ``\omega_{t+1}^{i}= \hat{\omega}_t^{i}``.  We then have the updated ensemble,
```math
\mu_{t+1}^{N}(dx) = \sum_{i=1}^N \omega_{t+1}^{i}\delta_{\xi_{t+1}^{i}}(dx)
```
See [Mutation](@ref mutation1) for additional information on implemenation.

### Unbiasedness 
Repeating the above alternating selection/mutation steps, the key feature of WE
is as follows. Denote by ``\mu`` the invariant measure associated with the Markov kernel ``K`` such that ``\mu K=\mu``.  Additionally, suppose ``\omega_0^{i}=1/N`` and ``\xi_0^{i}\sim \mu_0`` be i.i.d.  Lastly, denote by ``\mu_t = \mu_0 K^t``.  Then for any observable ``f``,
```math
\mu_t(f)=\mathbb{E}_{\mu_0}[f(X_t)] = \mathbb{E}[\mu_t^{N}(f)]= \mathbb{E}[\hat{\mu}_t^{N}(f)],
```
where the latter two expectations are ensemble averages over the WE process.  Consequently, WE is _unbiased_, and it is providing $N$-particle approximations of the distribution, ``\mu_t``.


### Estimation
WE can thus be used to compute averages against the associated a target stationary measure, ``\mu``, in the usual way:
```math
\mu(f) = \lim_{t\to\infty} \mu_t(f)= \lim_{t\to \infty}\frac{1}{t}\sum_{s=0}^{t-1} \mu_s(f) = \lim_{t\to \infty}\frac{1}{t}\sum_{s=0}^{t-1} \mu^{N}_s(f).
```
This is what we will use to make estimates; by taking ``t`` sufficiently large,
```math
\mu(f) \approx \frac{1}{t}\sum_{s=0}^{t-1} \mu^{N}_s(f)
```

A key thing to note is that WE does __not__ accelerate mixing.  The mixing rate is controlled by the underlying Markov process.  What WE does do is, ideally, provide very low variance estimates of ``\mu_t(f)``.  

WE can also be used to estimate finite time non-equilibrium quantities.

### Implementation
To run WE, we must define a [`WEsampler`](@ref) `sampler`, which defines the above steps.  Indeed, to recap, one must code up steps for
* `selection!(E, B, t)` - This is an in place transform that modifies the ensemble, `E`, and the bins `B`.  It may be time dependent.
* `mutation!(ξ)` - This is an in place transform corresponds to the application of the Markov kernel to the state, ``\xi \sim K(\xi,\bullet)``.  Though this operation is vectorized, it does defined on the underlying state space, not the ensemble or bin structures.  It is assumed to be time homogeneous.
* `sampler.rebin!(E, B, t+1)` - This is applied after mutation, at time ``t+1``, and determines which particles are in which bins, updating the structures accordingly.  It also retabulates the weight, ``\nu_t^i`` of each bin.

Having defined a sampler structure and provided initial `E0` ensemble and `B0` bins structures, WE is then executed with the command
```
n_we_steps = 10; 
E_trajectory, B_trajectory = run_we(E0, B0, sampler, n_we_steps);
```
The data structures contain, amongst other things,:
* `E_trajctory` - ``\{(\omega_t^i, \xi_t^i\}``, ``\{(\\hat{omega}_t^i, \hat{\xi}_t^i\}``;
* `B_trajectory` - ``n_t^i``, ``\nu_t^i``.
These can be used to compute the QoI.

## Caveats
There is some partial parallelization implemented, both with threads and with
`Distributed`.  In both cases, the parallelization is only with respect to the
application of the `mutation!` step, presumed to be the most costly piece of the
algorithm.  It does not currently make use of distributed data structures or
parallelize the selection step.

## Collaborators
* D. Aristoff
* J. Copperman
* L. F. Doherty
* F. G. Jones
* R. J. Webber
* D. M. Zuckerman

## Acknowledgements
This work was supported in part by the US National Science Foundation Grants DMS-1818716 and DMS-2111278.

