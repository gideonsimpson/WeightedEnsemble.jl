
using Distributed
using StatsBase
using HypothesisTests
using Printf

# addprocs(4);

@everywhere using NearestNeighbors

@everywhere include("schlogl_setup.jl");
@everywhere push!(LOAD_PATH,"../../src/");
@everywhere using JuWeightedEnsemble


# number of coarse steps in WE
@everywhere n_we_steps = 20;
# number of time steps during mutation step
@everywhere T_coarse = Tmax/n_we_steps;
# number of samples in coarse matrix
n_samples_per_bin = 10^2;
# ensemble size
n_particles = 10^2;

# define bin structure
@everywhere voronoi_pts = [[x] for x in 5:5:100];
B₀ = JuWeightedEnsemble.Voronoi_to_Bins(voronoi_pts);
@everywhere tree = KDTree(hcat(float.(voronoi_pts)...));

# define bin id mapping
@everywhere bin_id = x-> JuWeightedEnsemble.Voronoi_bin_id(x,tree);
# define the rebinning function
function rebin!(E, B, t)
    @. E.b = bin_id(E.ξ);
    JuWeightedEnsemble.update_bin_weights!(B, E);
    E, B
end

@everywhere mutation = x -> integrate_schlogl(x, T_coarse);
@everywhere function mutation!(x)
    y = integrate_schlogl(x, T_coarse)
    @. x = y;
    x
end

Random.seed!(100);
x0_vals = copy(voronoi_pts);
bin0_vals = bin_id.(voronoi_pts);
n_bins = length(B₀);
K̃ = JuWeightedEnsemble.pbuild_coarse_transition_matrix(mutation!, bin_id, x0_vals,bin0_vals, n_bins, n_samples_per_bin);

# define coarse observable as a bin function
f̃ = f.(voronoi_pts);
_,v²_vectors = JuWeightedEnsemble.build_coarse_vectors(n_we_steps,K̃,float.(f̃));
v² = (x,t)-> v²_vectors[t+1][bin_id(x)]
# define selection function
selection! = (E, B, t)-> JuWeightedEnsemble.optimal_allocation_selection!(E, B, v², t)

# set up ensemble
E₀ = JuWeightedEnsemble.Dirac_to_EnsembleWithBins(x₀, n_particles);
rebin!(E₀, B₀, 0);

# run
E = deepcopy(E₀);
B = deepcopy(B₀);
Random.seed!(200)
JuWeightedEnsemble.prun_we!(E, B, mutation,selection!, rebin!, n_we_steps);
p_est = f.(E.ξ) ⋅ E.ω
@printf("WE Estimate = %g\n", p_est)
