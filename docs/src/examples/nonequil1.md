# Nonequilibrium Sampling 

These exmaples involve sampling from a nonequilbrium distribution.  They are
particularly focused on computing the mean firest passage time (MFPT) via the
Hill relation [hill_free_1989](@cite); see, also,
<https://statisticalbiophysicsblog.org/?p=275>.  

## The Hill Relation
Consider the case of a discrete in time process, ``X_t``, and two sets, ``A``
and ``B``.  We wish to compute the mean first passage time for the process,
starting in ``A`` to reach ``B``.   One way of doing this would be sample
initial conditions in ``A`` (this distribution, ``\mu_0``, must be specified),
and run them until they reach ``B``.  But if  the MFPT is long, it will be
computationally prohibitive to observe enough transitions to obtain a
statistically meaningful estimate.

Now, we study an augmented processs, ``Y_t``, where ``Y_t`` obeys the same
dynamics as ``X_t`` until it reaches ``B``, at which point it resets in ``A``.
Denote the steady state distribution of the ``Y_t`` process by ``\nu``  The the
Hill relation then says
```math
\frac{1}{\text{MFPT for $X_\cdot$}} = \mathbb{P}_\nu(B)
```
Consequently, large MFPT's correspond to low probability events; some form of
importance sampling will be needed.  This is precisely the regime where WE can
pay off.  We will simulate the augmented process and estimate
``\mathbb{P}_\nu(B)`` using the we with the observable ``1_B(y)``.

To connect this to the continuous time case, frequently of interest, let us assume that ``X_t`` follows overdamped Langevin dynamics
```math
dX_t = -\nabla V(X_t)dt + \sqrt{2\beta^{-1}}dW_t,
```
The MFPT for this problem is computed first by solving the Poisson equation,
```math
L\tau_B = 0, \quad \tau_B|_{\partial B} = 0,
```
and then computing
```math
\text{MFPT $A\to B$} = \mathbb{E}_{\mu_0}[\tau_B(X)].
```
In the event that ``A=\{x_0\}``, as it will be in the computational example, this is just ``\tau_B(x_0)``.  As this Poisson problem will typically be in high dimension, solving it via classical methods is not feasable. 

The associated ``Y_t``, obeys this the overdamped Langevin equation until it
hits ``\partial B``, at which point it instantaneously retarts in ``A``
acoording ot the initial distribution ``\mu_0(dx) = \gamma_0(x)dx``.  The
modified process obeys the nonlocal Fokker-Planck equation:
```math
\partial_t \rho = L^\ast ρ + \gamma_0(x)\int_{\partial B}J_{\rho}\cdot n dS.
```
In the above expression, ``\int_{\partial B}J_{\rho}\cdot n dS`` is the
integrated probability flux into the set ``B``.  This has an associated noequilibrium steady state distribution (NESS).

To set this up to work with WE, we time discretize our process, ``\tilde{X}_{k}\approx X_{t_k}`` in some fashion and work with the associated augmented process (recylcing upon reaching ``B``), ``\tilde{Y}_k``.  The estimate of the MFPT of the original process is then
```math
\text{MFPT}\approx \frac{\Delta t}{\tilde{\nu}(B)}.
```
Indeed, this is approximate even without statistical error because we have lost
information on the sub ``\Delta t`` time scale.  But this error will typically
be small in comparison to the aforementioned statistical error. 

## Doublewell Example
```@contents
Pages = ["nonequil1.md"]
Depth = 3:4
```
Following the [equilibrium example](@ref double1), let us first compute, via ODE
methods, the NESS along with the MFPT.  Here, we will take ``B = [b,\infty)``,
with `b=0.5` and ``\mu_0 =\delta_{x_0}``, with `x0 = -1`.  Thus, we are
estimating the MFPT for a trajecotry, starting in the left basin to ge into the
"heart" of the right basin.  As we do not compute on an infinite domain, we
impose a Neumann (reflecting) boundary condition at ``x=a<x_0`` with `a=-2` in
the example.  For low enough temperature, this is a good approximation.


### ODE Results
For later comparison, we illustrate the results that are obtained usine ODE
methods.  The NESS density can be obtained by integration:
```@example 1
using QuadGK
using TestLandscapes
using Plots

β = 10.0;
V_ode(x) = SymmetricDoubleWell(x);
x0 = -1.0;
a = -2;
b = 0.5;

C0 = quadgk(t -> exp(β * V_ode(t)), x0, b)[1];
function f_(x)
    if (x < x0)
        return C0 * exp(-β * V_ode(x))
    else
        return quadgk(t -> exp(β * V_ode(t)), x, b)[1] * exp(-β * V_ode(x))
    end
end
Z = quadgk(f_, a, b)[1];
xx = LinRange(a, b, 200);
ρ = f_.(xx) / Z;

plot(xx, ρ, label="Density",lw=2)
xlabel!("x")
title!("NESS")
```
For the MFPT, we use `BVPProblem` from `BoundaryValueDiffEq`:
```@example 1
using BoundaryValueDiffEq
using ForwardDiff

∇V_ode = x->ForwardDiff.derivative(V_ode, x);

function τ_rhs!(du, u, p, t)
    du[1] = u[2]
    du[2] = β * (∇V_ode(t) * u[2] - 1.0)
    du
end

function τ_bc!(res, u, p, t)
    res[1] = u[1][2];
    res[2] = u[end][1];
    res
end

tspan = (a, b)

τ_bvp = BVProblem(τ_rhs!, τ_bc!, [1., 0.], tspan)    
τ_soln = solve(τ_bvp, MIRK4(), dt = 0.05);

MFPT = τ_soln(x0[1])[1];
@show MFPT;
```

### Defining the Mutation Step
In this example we use an Euler-Maruyama integrator with the recycling condition, with the `BasicMD` package:
```@example 1
using BasicMD
V(x) = SymmetricDoubleWell(x[1])
∇V! = (gradV, X) -> ForwardDiff.gradient!(gradV, V, X);
Δt = 1e-2;  # time step
Δt_recycle = 1e-2
nΔt_recycle = Int(Δt_recycle / Δt); # number of time steps before applying recycler
nΔt_coarse = 1 * nΔt_recycle # number of time steps in a coarse step

sampler = EM(∇V!, β, Δt);

mutation_opts = MDOptions(n_iters=nΔt_coarse, n_save_iters=nΔt_coarse); 

function restartA!(state::BasicMD.EMState)
    if (state.x[1] > b)
        @. state.x = [x0]
        ∇V!(state.∇V, [x0])
    end
    state
end

constraints = Constraints(restartA!, trivial_constraint!, nΔt_recycle, nΔt_coarse);

mutation! = x -> sample_trajectory!(x, sampler,constraints, options=mutation_opts); nothing
```

### Defining the Bins
Bins will be defined via Voronoi cells:
```@example 1
using WeightedEnsemble
Δb = 0.2
x_voronoi = [[x - 0.5 * Δb] for x in -1.5:Δb:0.5+Δb]
@show n_bins = length(x_voronoi);
B0, bin_id, rebin! = setup_Voronoi_bins(x_voronoi); nothing
```

### Run WE with Uniform Selection
As we will see, uniform selection will be sufficient to obtain relatively high quality results:
```@example 1
using Random

n_particles = 10^3;
E0 = Ensemble([[x_] for x_ in LinRange(-1.5, 0.5, n_particles)])
rebin!(E0, B0, 0);

uni_we_sampler = WEsampler(mutation!, uniform_selection!, rebin!);

n_we_steps = 10^3;

Random.seed!(100);
E_uni_trajectory, B_uni_trajectory = run_we(E0, B0, uni_we_sampler, n_we_steps); nothing
```
First, we verify convergence to the target NESS distribution:
```@example 1
using Printf
hist_bins = LinRange(-2.0, 1.5, 41);

we_anim = @animate for (t, E) in enumerate(E_uni_trajectory)
    histogram([ξ_[1] for ξ_ in E.ξ], bins=hist_bins, 
        weights=E.ω, norm=:pdf, alpha=0.5, label="Weighted Ensemble",
        yaxis=(:log10, [1e-6, 1e1]), legend=:topright)
    plot!(xx, ρ, lw=2, label="ODE")
    xlabel!("x")
    ylabel!("Density")
    xlims!(-1.5, 1)
    title!(@sprintf("t = %d", t))
end
gif(we_anim, fps =30)
```
Next, we take a look at our quantity of interest, the estimate of
``\mathbb{P}_{\tilde{\nu}}(B)``.  Using the Hill relation, this should
correspond to an estimate of ``\Delta t/\text{MFPT}``, the value we already
computed by ODE methods:
```@example 1
using LinearAlgebra

f_uni_trajectory = zeros(n_we_steps);

fB = X -> Float64(X[1] > b); # define observable

for (j,E) in enumerate(E_uni_trajectory)
    f_uni_trajectory[j] = (E.ω ⋅ fB.(E.ξ));
end

plot(1:n_we_steps, f_uni_trajectory,  lw=2,yaxis=(:log10, [1e-8, :auto]),label="WE Est.")
plot!(1:n_we_steps, cumsum(f_uni_trajectory) ./(1:n_we_steps),lw=2, label="WE Time Avg.")
plot!(1:n_we_steps, Δt / MFPT * ones(n_we_steps),
    lw=2,label="Δt/MFPT", color=:black, ls=:dash)
xlabel!("Iterate")
```
We have a good estimate after 500 iterations.

### Comparison with Direct Computation
As in the equilibrium case, an attempt to estimate the MFPT via the Hill
relation, without using WE, will result in very poor results:
```@example 1
Random.seed!(200)
trivial_we_sampler = WEsampler(mutation!, (E, B, t)->trivial_selection!(E), rebin!);
f_direct_trajectory = run_we_observables(E0, B0, trivial_we_sampler, n_we_steps, (fB,))[:];

plot(1:n_we_steps, f_uni_trajectory,  lw=2,yaxis=(:log10, [1e-8, :auto]), label="WE Est.")
plot!(1:n_we_steps, cumsum(f_uni_trajectory) ./(1:n_we_steps),lw=2, label="WE Time Avg.")
plot!(1:n_we_steps, f_direct_trajectory, lw=2, label="Direct Est.")
plot!(1:n_we_steps, cumsum(f_direct_trajectory) ./(1:n_we_steps),lw=2, label="Direct Time Avg.")
plot!(1:n_we_steps, Δt / MFPT * ones(n_we_steps),
    lw=2,label="Δt/MFPT", color=:black, ls=:dash)
xlabel!("Iterate")

```