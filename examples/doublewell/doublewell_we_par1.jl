#=
Parallel WE estimation of the probability of for a diffusion with X(0) = a
satisfying X(T) ∈ (b, ∞) for the double well potential V(x) = (x²-1)².

Either run as julia -p 4 doublewell_we_par1.jl or uncomment addprocs below.
=#

using Distributed
using Statistics
using HypothesisTests
using Printf

nw = 4; # number of workers
addprocs(nw);

@everywhere using WeightedEnsemble

@everywhere include("doublewell_setup.jl");

# number of coarse steps in WE
n_we_steps = 10;
# number of time steps during mutation step
nΔt_coarse = nΔt ÷ n_we_steps;
# number of samples in coarse matrix
n_samples_per_bin = 10^2;
# ensemble size
n_particles = 10^2;

# define bin structure
@everywhere voronoi_pts = [[x] for x in LinRange(a-.1,b+.1,21)];
@everywhere B₀, bin_id, rebin! = setup_Voronoi_bins(voronoi_pts);

# define the mutation mapping
opts = MDOptions(n_iters=nΔt_coarse, n_save_iters = nΔt_coarse)
@everywhere function mutation(x)
    Xvals, _ = sample_trajectory(x, sampler, options=opts);
    return Xvals[end]
end
@everywhere mutation! = x-> sample_trajectory!(x, sampler, options=opts);

# build the coarse model for optimal allocation
Random.seed!(100);
x0_vals = copy(voronoi_pts);
n_bins = length(B₀);
K̃ = WeightedEnsemble.pbuild_coarse_transition_matrix(mutation!, bin_id, x0_vals, n_bins, n_samples_per_bin);

# define coarse observable as a bin function
f̃ = f.(voronoi_pts);
_,v²_vectors = WeightedEnsemble.build_coarse_vectors(n_we_steps,K̃,float.(f̃));
v² = (x,t)-> v²_vectors[t+1][bin_id(x)]
selection! = (E, B, t)-> WeightedEnsemble.optimal_allocation!(E, B, v², t)
we_sampler = DistributedWEsampler(mutation!, selection!, rebin!);

# set up ensemble
E₀ = Dirac_to_Ensemble(x₀, n_particles);
rebin!(E₀, B₀, 0);

# run
E = deepcopy(E₀);
B = deepcopy(B₀);
Random.seed!(200)
prun_we!(E, B, we_sampler, n_we_steps);
p_est = f.(E.ξ) ⋅ E.ω
@printf("WE Estimate = %g\n", p_est)
