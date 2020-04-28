#=
WE estimation of the probability of for a diffusion with X(0) = a satisfying
X(T) ∈ (b, ∞) for the double well potential V(x) = (x²-1)².
=#


using StatsBase
using HypothesisTests
using Printf
using NearestNeighbors

include("doublewell_setup.jl");
push!(LOAD_PATH,"../../src/");
using JuWeightedEnsemble

# number of coarse steps in WE
n_we_steps = 10;
# number of time steps during mutation step
nΔt_coarse = ceil(Int,nΔt/n_we_steps);
# number of samples in coarse matrix
n_samples_per_bin = 10^2;
# ensemble size
n_particles = 10^2;

# define bin structure
voronoi_pts = [[x] for x in LinRange(a-.1,b+.1,21)];
B₀ = JuWeightedEnsemble.Voronoi_to_Bins(voronoi_pts);
tree = KDTree(hcat(voronoi_pts...));

# define bin id mapping
bin_id = x-> JuWeightedEnsemble.Voronoi_bin_id(x,tree);
# define the rebinning function
function rebin!(E, B, t)
    @. E.b = bin_id(E.ξ);
    JuWeightedEnsemble.update_bin_weights!(B, E);
    E, B
end

# define the mutation mapping
opts = MDOptions(n_iters=nΔt_coarse, n_save_iters = nΔt_coarse)
function mutation(x)
    Xvals, _ = sample_trajectory(x, sampler, options=opts);
    return Xvals[end]
end
mutation! = x-> sample_trajectory!(x, sampler, options=opts);

# construct coarse model matrix
Random.seed!(100);
x0_vals = copy(voronoi_pts);
bin0_vals = bin_id.(voronoi_pts);
n_bins = length(B₀);
K̃ = JuWeightedEnsemble.build_coarse_transition_matrix(mutation!, bin_id, x0_vals,bin0_vals, n_bins, n_samples_per_bin);

# define coarse observable as a bin function
f̃ = f.(voronoi_pts);
_,v²_vectors = JuWeightedEnsemble.build_coarse_vectors(n_we_steps,K̃,float.(f̃));
v² = (x,t)-> v²_vectors[t+1][bin_id(x)]
# define selection function
selection! = (E, B, t)-> JuWeightedEnsemble.optimal_allocation_selection!(E, B, v², t)
# selection! = (E, B, t)-> JuWeightedEnsemble.uniform_selection!(E, B);

# set up ensemble
E₀ = JuWeightedEnsemble.Dirac_to_EnsembleWithBins(x₀, n_particles);
rebin!(E₀, B₀, 0);

# run
E = deepcopy(E₀);
B = deepcopy(B₀);
Random.seed!(200)
JuWeightedEnsemble.run_we!(E, B, mutation,selection!, rebin!, n_we_steps);
p_est = f.(E.ξ) ⋅ E.ω
@printf("WE Estimate = %g\n", p_est)
