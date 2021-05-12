#=
Serial WE estimation of the probability for the Schlögl chemical reaction network
with X(0)= x₀ to reach X(T) ∈ A, the target set.

This reaction network corresponds to:
A + 2S → 3S
3S → A + 2S
B → S
S → B

Parameters were chosen so that the system is bistable, as in:
https://doi.org/10.1063/1.5017955
=#

using Statistics
using HypothesisTests
using Printf
using NearestNeighbors
using LinearAlgebra

include("schlogl_setup.jl");
using WeightedEnsemble


# number of coarse steps in WE
n_we_steps = 20;
# number of time steps during mutation step
T_coarse = Tmax/n_we_steps;
# number of samples in coarse matrix
n_samples_per_bin = 10^2;
# ensemble size
n_particles = 10^2;

# define bin structure
voronoi_pts = [[x] for x in 5:5:100];
B₀ = WeightedEnsemble.Voronoi_to_Bins(voronoi_pts);
tree = KDTree(hcat(float.(voronoi_pts)...));

# define bin id mapping
bin_id = x-> WeightedEnsemble.Voronoi_bin_id(x,tree);
# define the rebinning function
function rebin!(E, B, t)
    @. E.b = bin_id(E.ξ);
    WeightedEnsemble.update_bin_weights!(B, E);
    E, B
end

function mutation!(x)
    y = integrate_schlogl(x, T_coarse)
    @. x = y;
    x
end


# construct coarse model matrix
Random.seed!(100);
x0_vals = copy(voronoi_pts);
bin0_vals = bin_id.(voronoi_pts);
n_bins = length(B₀);
K̃ = WeightedEnsemble.build_coarse_transition_matrix(mutation!, bin_id, x0_vals,bin0_vals, n_bins, n_samples_per_bin);

# define coarse observable as a bin function
f̃ = f.(voronoi_pts);
_,v²_vectors = WeightedEnsemble.build_coarse_vectors(n_we_steps,K̃,float.(f̃));
v² = (x,t)-> v²_vectors[t+1][bin_id(x)]
# define selection function
selection! = (E, B, t)-> WeightedEnsemble.optimal_allocation!(E, B, v², t)

# set up ensemble
E₀ = WeightedEnsemble.Dirac_to_EnsembleWithBins(x₀, n_particles);
rebin!(E₀, B₀, 0);

# run
E = deepcopy(E₀);
B = deepcopy(B₀);
WeightedEnsemble.trun_we!(E, B, mutation!, selection!, rebin!, n_we_steps);
p_est = f.(E.ξ) ⋅ E.ω
@printf("WE Estimate = %g\n", p_est)
