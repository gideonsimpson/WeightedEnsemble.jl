using Random
using LinearAlgebra

Δb = 0.2
n_particles = 10^5;
seed = 100;
β = 10.0    # inverse temperature
b = 0.5 # target set [b, ∞)
x0 = [-1.0] # starting point

T = 4 # terminal time


Δt = 1e-2  # time step
nΔt = Int(T / Δt)
ΔT_recycle = 1e-2
nΔt_recycle = Int(ΔT_recycle / Δt) # number of time steps before applying recycler
nΔt_coarse = 1 * nΔt_recycle # number of time steps in a coarse step
@show n_we_steps = nΔt ÷ (nΔt_coarse);

x_voronoi = [[x - 0.5 * Δb] for x in -1.7:Δb:0.5+Δb]
n_bins = length(x_voronoi)

@show T;
@show Δt;
@show β;
@show n_particles;
@show Δb;
@show n_bins;
@show x_voronoi;

include("setup1.jl");

uni_we_sampler = WEsampler(mutation!, uniform_selection!, rebin!);

Random.seed!(seed);
# E_trajectory, B_trajectory = run_we(E0, B0, uni_we_sampler, n_we_steps);
f_trajectory = run_we_observables(E0, B0, uni_we_sampler, n_we_steps - 1, (fB,))[:];
prepend!(f_trajectory, fB.(E0.ξ) ⋅ E0.ω);

