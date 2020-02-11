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
nΔt_coarse = Int(nΔt/n_we_steps);
# number of samples in coarse matrix
n_samples_per_bin = 10^2;
# ensemble size
n_particles = 10^2;

# define bin structure
voronoi_pts = [[x] for x in LinRange(a-.1,b+.1,21)];
B₀ = JuWeightedEnsemble.Voronoi_to_Bins(voronoi_pts);
tree = KDTree(hcat(voronoi_pts...));

# define the mutation mapping
mutation = x-> MALA(x, V, gradV!, β, Δt, nΔt_coarse, return_trajectory=false)[1];
mutation! = x-> MALA!(x, V, gradV!, β, Δt, nΔt_coarse);

# define bin id mapping
bin_id = x-> JuWeightedEnsemble.Voronoi_bin_id(x,tree);

# construct coarse model matrix
Random.seed!(100);
x0_vals = copy(voronoi_pts);
bin0_vals = bin_id.(voronoi_pts);
n_bins = length(B₀);
T = JuWeightedEnsemble.build_coarse_transition_matrix(mutation!, bin_id, x0_vals,bin0_vals, n_bins, n_samples_per_bin);

# define coarse observable as a bin function
F = f.(voronoi_pts);
value_vectors = JuWeightedEnsemble.build_value_vectors(n_we_steps,T,float.(F));

#  define selection function
selection! = (E, B, j)-> JuWeightedEnsemble.optimal_allocation_selection!(E,B,value_vectors,j)
# selection! = (E, B, j)-> JuWeightedEnsemble.uniform_selection!(E,B);

# set up ensemble
ξ₀ = [copy(x₀) for i = 1:n_particles];
ω₀ = 1.0/n_particles * ones(n_particles);

E₀ = Ensemble{Array{Float64,1}, Float64, Int}(copy(ξ₀),copy(ξ₀),copy(ω₀), copy(ω₀),
                            zeros(Int, n_particles),zeros(Int, n_particles));
@. E₀.bin =bin_id(E₀.ξ);
JuWeightedEnsemble.update_bin_weights!(B₀, E₀);

# run
E = deepcopy(E₀);
B = deepcopy(B₀);
Random.seed!(200)
JuWeightedEnsemble.run_we!(E, B, mutation,selection!, bin_id, n_we_steps);
p_est = f.(E.ξ) ⋅ E.ω
@printf("WE Estimate = %g\n", p_est)
