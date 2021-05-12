using Statistics
using HypothesisTests
using Printf
using Distributed

addprocs(4)

@everywhere using WeightedEnsemble

@everywhere include("doublewell_setup.jl");

@everywhere begin
    Δt = 1e-4;  # time step
    tmax = 100.0;    # terminal time
    nΔt = ceil(Int, tmax/Δt);
end

# number of coarse steps in WE
@everywhere n_we_steps = 1000;
# number of time steps during mutation step
@everywhere nΔt_coarse =nΔt ÷ n_we_steps;
# number of samples in coarse matrix
n_samples_per_bin = 10^2;
# ensemble size
n_particles = 10^2;

# define bin structure
@everywhere voronoi_pts = [[x] for x in LinRange(-1.2,1.2,21)];
@everywhere B₀, bin_id, rebin! = setup_Voronoi_bins(voronoi_pts);

# define the mutation mapping
@everywhere opts = MDOptions(n_iters=nΔt_coarse, n_save_iters = nΔt_coarse)
@everywhere mutation! = x-> sample_trajectory!(x, sampler, options=opts);
@everywhere function mutation(x)
    Xvals, _ = sample_trajectory(x, sampler, options=opts);
    return Xvals[end]
end
selection! = (E, B, t)-> WeightedEnsemble.uniform_allocation!(E, B);

f1 = x-> Float64(-0.1 < x[1] < 0.1); # define observables
f2 = x-> Float64( x[1] > 1.1);
obs = (f1, f2);

# set up ensemble
E₀ = WeightedEnsemble.Dirac_to_EnsembleWithBins(x₀, n_particles);
rebin!(E₀, B₀, 0);

Random.seed!(200)
obs_trajectory = prun_we_observables(E₀, B₀, mutation, selection!, rebin!, n_we_steps, obs);
