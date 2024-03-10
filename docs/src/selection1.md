# Particle Selection Algorithms

```@contents
Pages = ["selection1.md"]
Depth = 2
```


## More about the Selection Step
The selection step, encoded in a user defined `selection!` function, is the
essential element of WE.  It is designed to ensure that the results of WE are
unbiased _and_ work to reduce the variance of some QoI.  

Generically, when we call `selection!(E, B, t)`, at algorthmic step `t`, we
advance from ``\{(\omega_t^i, \xi_t^i\}`` (before selection) to
``\{(\hat{\omega}_t^i, \hat{\xi}_t^i\}`` (after selection).  

For algorithmic consistency and numerical stability, the selection step
typically includes the following two conditions:
* Every non-empty bin __must__ have at least one offspring.  Thus, if ``u_p`` is
  occupied before selection (there exists ``i`` such that ``\xi_t^i \in u_p``),
  then after selection, there exists ``\hat{i}`` such that
  ``\hat{\xi}_t^{\hat{i}} \in u_p``.  Hence, if we have a total of ``N``
  particles, and ``K`` nonempty bins, we immediately allocate 1 particle to each
  of the nonempty bins, and then must allocate the remaining ``N-K`` particles.
* If the total mass in a bin falls beneath some threshold, ``\nu_t^p
  \leq \nu_{\min}``, then the number of offspring of that bin equals the number
  of particles currently in the bin:
  ```math
  n_t^p = \sum_{\xi_t^i \in \mathcal{B}_p} 1 = \sum_{\hat{\xi}_t^i \in \mathcal{B}_p}1
  ```
  This is neccessary to avoid floating point underflow issues.
  
These steps are both handled by `minimal_bin_allocation!(E, B)`.  After this
step is completed, the remaining  particles are allocated to the bins, then
particles are allocated within each bin, and, finally, [`repopulate!(E, B)`](@ref) is
called, which assigns ``\{(\hat{\omega}_t^i, \hat{\xi}_t^i\}``, storing them in
`E.ξ̂` and `E.ω̂`.

As several of the included selection routines show, one may include additional
arguments.  But when the `selection!` step is used to define a
[`WEsampler`](@ref) structure, it must be such that it only takes as arguments
`(E,B,t)`.  For instance, if we want to use `targeted_selection!`, we need to
define a target function, 
```
selection! = (E, B, t) -> targeted_selection!(E, B, G, t);
sampler = WEsampler(mutation!, selection!, rebin!); # mutation! and rebin! defined elsewhere
```

## Included Selection Methods

```@docs
uniform_selection!
```
Uniform selection is suboptimal (in the sense of variance reduction), but
provided there are a sufficient number of bins and particles, it can often
provide quite satisfactory results.  It has the advantage of being very
straightforward to implement.  

```@docs
optimal_selection!
```
Optimal selection is taken from Aristoff & Zuckerman (2020), where particles are allocated to bins according to 
```math
\frac{C_t(u)}{N} \propto \nu_t \sqrt{\sum_{\xi_t^i\mid i\in u}\frac{\omega_t^i}{\nu_t}v^2(\xi_t^i,t)}
```
For steady state problems 
```math
v^2(\xi, t) = v^2(\xi) =\mathrm{Var}_{K(\xi,\bullet)}(h),
```
where ``h`` solves the Poisson equation,
```math
(I-K)h = f - \mu(f).
```
This of course, must be estimated, and tools are provided for this in [Coarse Models](@ref)


```@docs
targeted_selection!
```
Targeted selection allows one steer bin allocation according to the rule
```math
\frac{C_t(u)}{N} \propto G(u)
```
where the target function, ``G``, is of the form `G(p, E, B, t)`, with bin index `p`.  


```@docs
static_selection!
```
Static selection allows one to predetermine the target number of particles in
each bin by providing a vector `n_static`.  It is assumed `n_static[i]>0` and
`sum(n_static)≤N`.  It may be that fewer than `N` particles are allocated.  In
this case, some of the offspring will be given zero mass.

```@docs
trivial_selection!
```
This a convenience tool built in which ``\{(\omega_t^i, \xi_t^i\}=\{(\hat{\omega}_t^i, \hat{\xi}_t^i\}``.  It can be useful when benchmarking against the non-interacting particle system case.