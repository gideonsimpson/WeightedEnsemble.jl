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
While not the rarest of events, this will take a bit of effort to estimate with
high confidence.  This will be estimated using WE with the associated QoI, ``1_A(x)``.

We will explore this landscape using the MALA approximation of overdamped
Langevin dyanmics; thus, we build up the mutation step:
```@example 1
using ForwardDiff
using BasicMD
V(x) = SymmetricDoubleWell(x[1])
∇V! = (gradV, X) -> ForwardDiff.gradient!(gradV, V, X);
Δt = 1e-2;  # time step
sampler = MALA(V, ∇V!, β, Δt);

nΔt = 10; # number of fine time steps per coarse WE step

# define the mutation mapping
opts = MDOptions(n_iters=nΔt);
mutation! = x -> sample_trajectory!(x, sampler, options=opts); nothing
```
We define our observable, the indicator function on the
set ``A``:
```@example 1
fA(x) = Float64(-δ<x[1]<δ); nothing
```
Next, we select a set of bins based on a Voronoi tesselation of the real line:
```@example 1
using WeightedEnsemble
voronoi_pts = [[x_] for x_ in LinRange(-1.5, 1.5, 21)];
B0, bin_id, rebin! = setup_Voronoi_bins(voronoi_pts); nothing
```
before defining the selection step, we will also set the number of particles, WE
steps, and initialize an ensemble:
```@example 1
n_we_steps = 10^2;
n_particles = 10^3;

# the following is somwhat arbitrary, but it spreads particles through state space
E0 = Ensemble([[x_] for x_ in LinRange(-1., 1., n_particles)]);
rebin!(E0, B0, 0); nothing # ensure weights are properly tabulated
```

### Uniform Selection
In the case that we use uniform selection, we proceed with:
```@example 1
uni_selection! = (E, B, t) -> uniform_selection!(E, B, t)
uni_we_sampler = WEsampler(mutation!, uni_selection!, rebin!); nothing
```
Next, we can run and examine our results:
```@example 1
using Statistics
using Random # for reproducbility
Random.seed!(100);
f_trajectory = run_we_observables(E0, B0, uni_we_sampler, n_we_steps, (fA,))[:]
@show mean(f_trajectory);
```
Visualizing our computation:
```@example 1
using Plots
using Printf
plot(1:n_we_steps, f_trajectory , yaxis=(:log10, [1e-6, :auto]),lw=2,label="WE Est.")
plot!(1:n_we_steps, cumsum(f_trajectory) ./ (1:n_we_steps),lw=2, label="WE Time Avg.")
plot!(1:n_we_steps, p * ones(n_we_steps),lw=2,ls=:dash,color=:black, label="Exact")
xlabel!("Iterate")
```
This gives a fairly good estimate, and if we were to apply a burnin rule, discarding the first 25 over so iterates, the time average would be spot on, as shown by the following, naive, 95% confidence interval calculation.
```@example 1
using HypothesisTests
confint(OneSampleTTest(f_trajectory[26:n_we_steps]))
```

### Visualizing the WE
To get a good sense of why WE works, we can animate the empirical distribution
as a function of time:
```@example 1
Random.seed!(100);
E_trajectory, _ = run_we(E0, B0, uni_we_sampler, n_we_steps);

xplt = LinRange(-1.5, 1.5,501); 
hist_bins = LinRange(-1.5,1.5,31);

we_anim = @animate for (t, E) in enumerate(E_trajectory)
    histogram([ξ_[1] for ξ_ in E.ξ], bins=hist_bins, 
        weights=E.ω, norm=:pdf, alpha=0.5, label="Weighted Ensemble", legend=:top)
    histogram!([ξ_[1] for ξ_ in E.ξ], bins=hist_bins,
     norm=:pdf, alpha=0.5, label="Unweighted Ensemble")
    plot!(xplt, fA.(xplt), label="QoI",lw=2)
    plot!(xplt, f.(xplt)/Z,label="Analytic",lw=2)
    xlabel!("x")
    ylabel!("Density")
    xlims!(-1.5, 1.5)
    ylims!(0,2)
    title!(@sprintf("t = %d", t))
end
gif(we_anim, fps =6)
```
In the above animation, the `Weighted Ensemble` histogram corresponds to the distribution
```math
\sum_i \omega_t^{i} \delta_{\xi_t^{i}}
```
while the `Unweighted Ensemble` histogram corresponds to
```math
\sum_i \tfrac{1}{N} \delta_{\xi_t^{i}}
```
This shows that WE's selection step makes an effort to keep particles in the
support of the QoI, and by weighting these particles properly, we sample the
correct distribution.



### Comparison with Direct Computation
One might ask how it fares against having used the same amount of resources on
the problem directly:
```@example 1
trivial_we_sampler = WEsampler(mutation!, (E, B, t)->trivial_selection!(E), rebin!);
f_direct_trajectory = run_we_observables(E0, B0, trivial_we_sampler, n_we_steps, (fA,))[:];

plot(1:n_we_steps, f_trajectory ,yaxis=(:log10, [1e-6, :auto]),lw=2,label="WE Est.")
plot!(1:n_we_steps, cumsum(f_trajectory) ./ (1:n_we_steps),lw=2, label="WE Time Avg.")
plot!(1:n_we_steps, f_direct_trajectory , yscale=:log10,lw=2,label="Direct Est.")
plot!(1:n_we_steps, cumsum(f_direct_trajectory) ./ (1:n_we_steps),lw=2, label="Direct Time Avg.")
plot!(1:n_we_steps, p * ones(n_we_steps),lw=2,ls=:dash,color=:black, label="Exact")
xlabel!("Iterate")
```
The direct computation, which involves using `n_particles` independent trials,
each run for `n_we_steps` struggles to estimate the QoI.  Indeed, the reason the
line labeled `Direct Est.` appears to stop is because it is returning zeros
which are not showing up on the logarithmic scale.

### Optimal Selection
TBW


