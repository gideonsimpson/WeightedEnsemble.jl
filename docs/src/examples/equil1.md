# Equilibrium Sampling 

These exmaples involve sampling from an equilibirum Boltzmann type distribution
of the sort commonly found in molecular dynamics and Bayesian inverse problems.  


## [Doublewell Potential](@id double1)
```@contents
Pages = ["equil1.md"]
Depth = 3:4
```

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

### Defining the Mutation Step
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
### Defining the Observable
We define our observable, the indicator function on the
set ``A``:
```@example 1
fA(x) = Float64(-δ<x[1]<δ); nothing
```
### Defining the Bins
Next, we select a set of bins based on a Voronoi tesselation of the real line:
```@example 1
using WeightedEnsemble
voronoi_pts = [[x_] for x_ in LinRange(-1.5, 1.5, 13)];
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
f_uni_trajectory = run_we_observables(E0, B0, uni_we_sampler, n_we_steps, (fA,))[:]
@show mean(f_uni_trajectory);
```
Visualizing our computation:
```@example 1
using Plots
using Printf
plot(1:n_we_steps, f_uni_trajectory , yaxis=(:log10, [1e-6, :auto]),lw=2,label="WE Est.")
plot!(1:n_we_steps, cumsum(f_uni_trajectory) ./ (1:n_we_steps),lw=2, label="WE Time Avg.")
plot!(1:n_we_steps, p * ones(n_we_steps),lw=2,ls=:dash,color=:black, label="Exact")
xlabel!("Iterate")
```
This gives a fairly good estimate, and if we were to apply a burnin rule, discarding the first 25 over so iterates, the time average would be spot on, as shown by the following, naive, 95% confidence interval calculation.
```@example 1
using HypothesisTests
confint(OneSampleTTest(f_uni_trajectory[26:n_we_steps]))
```

### Visualizing WE with Uniform Selection
To get a good sense of why WE works, we can animate the empirical distribution
as a function of time:
```@example 1
Random.seed!(100);
E_uni_trajectory, _ = run_we(E0, B0, uni_we_sampler, n_we_steps);

xplt = LinRange(-1.5, 1.5,501); 
hist_bins = LinRange(-1.5,1.5,31);

we_anim = @animate for (t, E) in enumerate(E_uni_trajectory)
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



### Comparison of Uniform WE with Direct Computation
One might ask how it fares against having used the same amount of resources on
the problem directly:
```@example 1
Random.seed!(200)
trivial_we_sampler = WEsampler(mutation!, (E, B, t)->trivial_selection!(E), rebin!);
f_direct_trajectory = run_we_observables(E0, B0, trivial_we_sampler, n_we_steps, (fA,))[:];

plot(1:n_we_steps, f_uni_trajectory ,yaxis=(:log10, [1e-6, :auto]),lw=2,label="WE Est.")
plot!(1:n_we_steps, cumsum(f_uni_trajectory) ./ (1:n_we_steps),lw=2, label="WE Time Avg.")
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
To perform Optimal Selection, in the spirit of
[aristoff_optimizing_2020](@cite), we must first approximate
```math
v^2(x) = \mathrm{Var}_{K(x,\bullet)}(h)
```
where ``K`` is our Markov transition kernel, and ``h`` is the solution of the
Poisson problem
```math
(I - K)h = f - \mu(f)
```
An approximation is found by first estimating coarse, bin scale, transition, matrix ``\tilde{K}``, where
```math
\tilde{K}_{pq} \approx \mathbb{P}(p\to q),
```
the probability of a particle starting in bin ``p`` transitioning to bin ``q``.
Left unsaid is what hte initial distribution is for the particles; in practice,
we often use a Dirac (centered at the Voronoi cell center).  This matrix is
estimated using our coarse model tools (see [Coarse Models](@ref)):
```@example 1
x0_vals = deepcopy(voronoi_pts); # Dirac at the bin
n_bins = length(B0);
n_samples_per_bin = 10^3;

Random.seed!(5000)
K̃ = WeightedEnsemble.build_coarse_transition_matrix(mutation!, bin_id,
    x0_vals, n_bins, n_samples_per_bin);
```
Having constructed this, we now obtain the coarse, bin scale, approximation of
``v^2``.  This requires us to first construct a bin function approximation of
the observable, then compose it with the bin identification function, ``\tilde{f}_A``.
```@example 1
f̃A = fA.(voronoi_pts); # bin function for 
h̃, ṽ² = WeightedEnsemble.build_coarse_poisson(K̃, f̃A);
v² = (x, t) -> ṽ²[bin_id(x)]; nothing
```
We are now ready to run the optimized sampler:
```@example 1
optimal! = (E, B, t) -> optimal_selection!(E, B, v², t);
opt_we_sampler = WEsampler(mutation!, optimal!, rebin!);
Random.seed!(300)
f_opt_trajectory = run_we_observables(E0, B0, opt_we_sampler, n_we_steps, (fA,))[:];

plot(1:n_we_steps, f_uni_trajectory ,yaxis=(:log10, [1e-6, :auto]),lw=2,label="WE Est.")
plot!(1:n_we_steps, cumsum(f_uni_trajectory) ./ (1:n_we_steps),lw=2, label="WE Time Avg.")
plot!(1:n_we_steps, f_opt_trajectory , yscale=:log10,lw=2,label="Optimal WE Est.")
plot!(1:n_we_steps, cumsum(f_opt_trajectory) ./ (1:n_we_steps),lw=2, label="Optimal WE Time Avg.")
plot!(1:n_we_steps, p * ones(n_we_steps),lw=2,ls=:dash,color=:black, label="Exact")
xlabel!("Iterate")
```
The difference in this example is neglible to the eye, though when we check the variance (after burnin), we do see an improvement:
```@example 1
println(var(f_opt_trajectory[26:end]));
println(var(f_uni_trajectory[26:end]));
```
More substantial improvements can be found in other problems, such as those shown in [aristoff_weighted_2023](@cite).