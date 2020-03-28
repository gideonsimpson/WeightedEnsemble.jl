#=
WE estimation of the probability of for a diffusion with X(0) = x₀ satisfying
X(T) ∈ B for the Muller potential.
=#


using StatsBase
using HypothesisTests
using Printf
using NearestNeighbors

include("muller_setup.jl");
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

# define bin structure using Voronoi
xc = LinRange(-1.5,1,7)
yc = LinRange(-0.5,2,7)
voronoi_pts = Array{Float64,1}[];
for x in xc, y in yc
    # only include points that are likely to be accessed
    if(V([x,y])<250)
        push!(voronoi_pts, [x,y])
    end
end
B₀ = JuWeightedEnsemble.Voronoi_to_Bins(voronoi_pts);
tree = KDTree(hcat(voronoi_pts...));

bin_id = x-> JuWeightedEnsemble.Voronoi_bin_id(x,tree);
# define the rebinning function
function rebin!(E, B, t)
    @. E.b = bin_id(E.ξ);
    JuWeightedEnsemble.update_bin_weights!(B, E);
    E, B
end

# mutation = X-> EM(X,∇V!,β, Δt, Int(nΔt/n_we_steps), return_trajectory=false);
# mutation! = X-> EM!(X,∇V!,β, Δt, Int(nΔt/n_we_steps));
mutation = X-> MALA(X,V, ∇V!,β, Δt, Int(nΔt/n_we_steps), return_trajectory=false)[1];
mutation! = X-> MALA!(X,V, ∇V!,β, Δt, Int(nΔt/n_we_steps));

# construct coarse model
Random.seed!(100);
x0_vals = copy(voronoi_pts);
bin0_vals = bin_id.(voronoi_pts);
n_bins = length(B₀);
T = JuWeightedEnsemble.build_coarse_transition_matrix(mutation!, bin_id, x0_vals,bin0_vals, n_bins, n_samples_per_bin);

# define coarse observable as a bin function
F = f.(voronoi_pts);
value_vectors = JuWeightedEnsemble.build_value_vectors(n_we_steps,T,float.(F));
h = (x,t)-> value_vectors[t][bin_id(x)];
selection! = (E, B, t)-> JuWeightedEnsemble.optimal_allocation_selection!(E,B,h,t);

# set up ensemble
ξ₀ = [copy(x₀) for i = 1:n_particles];
ω₀ = 1.0/n_particles * ones(n_particles);

E₀ = Ensemble{Array{Float64,1}, Float64, Int}(copy(ξ₀),copy(ξ₀),copy(ω₀), copy(ω₀),
                            zeros(Int, n_particles),zeros(Int, n_particles),zeros(Int, n_particles));

rebin!(E₀, B₀, 0);

# run
E = deepcopy(E₀);
B = deepcopy(B₀);
Random.seed!(200)
JuWeightedEnsemble.run_we!(E, B, mutation,selection!, rebin!, n_we_steps);
p_est = f.(E.ξ) ⋅ E.ω
@printf("WE Estimate = %g\n", p_est)
