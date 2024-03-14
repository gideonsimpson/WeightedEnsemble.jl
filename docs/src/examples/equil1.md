# Equilibrium Sampling 

These exmaples involve sampling from an equilibirum Boltzmann type distribution
of the sort commonly found in molecular dynamics and Bayesian inverse problems.  


## Doublewell Potential
As a first example, we will consider the classical double well potential,
```math
V(x) = (x^2-1)^2
```
in dimension one, and its associated Boltzmann distribution,
```math
\mu(dx) \propto e^{-\beta V(x)}dx.
```
For a sufficiently low temperature (large ``\beta``), we expect the distribution
to concentrate near ±1.  

Consider the question of estimating ``\mathbb{P}_\mu(A)``, where ``A =
(-\delta, \delta)``, is a set in the neighborhood of the saddle at the origin.
At low temperature, this will be a rare event, ``\mathbb{P}_\mu(A)\ll 1``.  Indeed, even at `β=10` and `δ=0.1`, we find:
```@example 1
using QuadGK
using TestLandscapes
β = 10.;
δ = 0.1;
f(x) = exp(-β*SymmetricDoubleWell(x))
Z = quadgk(f, -Inf, Inf)[1];
p = quadgk(f, -δ, δ)[1]/Z;
@show p;
```
While not the rareest of events, this will take a bit of effort to estimate with high confidence.

We will explore this landscape using the MALA approximation of overdamped
Langevin dyanmics; thus, we build up the mutation step:
```@example 1
using ForwardDiff
using BasicMD
V(x) = SymmetricDoubleWell(x[1])
∇V! = (gradV, X) -> ForwardDiff.gradient!(gradV, V, X);
Δt = 1e-2;  # time step
sampler = MALA(V, ∇V!, β, Δt);

nΔt = 1; # number of fine time steps per coarse WE step

# define the mutation mapping
opts = MDOptions(n_iters=nΔt);
mutation! = x -> sample_trajectory!(x, sampler, options=opts); nothing
```
Next, we define our observable, corresponding to the indicator function on the
set ``A``:
```@example 1
fA(x) = Float64(-δ<x[1]<δ); nothing
```
Next, we select a set of bins based on a Voronoi tesselation of the real line:
```@example 1
using WeightedEnsemble
voronoi_pts = [[x_] for x_ in LinRange(-1.5, 1.5, 23)];
B0, bin_id, rebin! = setup_Voronoi_bins(voronoi_pts); nothing
```
before defining the selection step, we will also set the number of particles, WE
steps, and initialize an ensemble:
```@example 1
n_we_steps = 10^3;
n_particles = 10^3;

x0 = [-1.0]; # this is somewhat arbitrary
E0 = Dirac_to_Ensemble(x0, n_particles);
rebin!(E0, B0, 0); nothing # ensure weights are properly tabulated
```


### Uniform Selection
In the case that we use uniform selection, we proceed with:
```@example 1
# define selection function
uni_selection! = (E, B, t) -> uniform_selection!(E, B, t)
# set up sampler - trivial analysis step
uni_we_sampler = WEsampler(mutation!, uni_selection!, rebin!); nothing
```
Next, we can run and examine our results:
```@example 1
using Statistics
using Random # for reproducbility
Random.seed!(100);
# f_trajectory = run_we_observables(E0, B0, uni_we_sampler, n_we_steps, (fA,))[:]
# @show mean(f_trajectory);
```
Visualizing our computation:
<!-- ```
# @example 1
# using Plots
# plot(1:n_we_steps, cumsum(f_trajectory) ./(1:n_we_steps),label="Uniform")
# plot!(1:n_we_steps, p*ones(n_we_steps),label="Exact",ls=:dash, color=:black)
# xlabel!("Iterate")
# ylabel!("Running Average")
``` -->
This gives a fairly good estimate, but one might ask how it fares against having used the same amount of resources on the problem directly:
<!-- ```@example 1
# set up sampler - trivial analysis step
# trivial_we_sampler = WEsampler(mutation!, (E, B, t)->trivial_selection!(E), rebin!); nothing
# f_direct_trajectory = run_we_observables(E0, B0, trivial_we_sampler, n_we_steps, (fA,))[:];
# @show mean(f_direct_trajectory);
# plot(1:n_we_steps, cumsum(f_trajectory) ./(1:n_we_steps),label="Uniform")
# plot!(1:n_we_steps, cumsum(f_direct_trajectory) ./(1:n_we_steps),label="Direct")
# plot!(1:n_we_steps, p*ones(n_we_steps),label="Exact",ls=:dash, color=:black)
# xlabel!("Iterate")
# ylabel!("Running Average")
``` -->

### Optimal Selection
TBW


