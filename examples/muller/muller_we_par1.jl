#=
Parallel WE estimation of the probability of for a diffusion with X(0) = x₀
satisfying X(T) ∈ B for the Muller potential.
=#

using Distributed
using StatsBase
using HypothesisTests
using Printf

# addprocs(4);

@everywhere using NearestNeighbors

@everywhere include("muller_setup.jl");
@everywhere push!(LOAD_PATH,"../../src/");
@everywhere using WeightedEnsemble

# number of coarse steps in WE
n_we_steps = 10;
# number of time steps during mutation step
nΔt_coarse = ceil(Int, nΔt/n_we_steps);
# number of samples in coarse matrix
n_samples_per_bin = 10^2;
# ensemble size
n_particles = 10^2;

# define bin structure
@everywhere xc = LinRange(-1.5,1,7)
@everywhere yc = LinRange(-0.5,2,7)
@everywhere voronoi_pts = Array{Float64,1}[];
@everywhere for x in xc, y in yc
    if(V([x,y])<250)
        push!(voronoi_pts, [x,y])
    end
end
B₀ = WeightedEnsemble.Voronoi_to_Bins(voronoi_pts);
@everywhere tree = KDTree(hcat(voronoi_pts...));

# define bin id mapping
@everywhere bin_id = x-> WeightedEnsemble.Voronoi_bin_id(x,tree);
# define the rebinning function
function rebin!(E, B, t)
    @. E.b = bin_id(E.ξ);
    WeightedEnsemble.update_bin_weights!(B, E);
    E, B
end

# define the mutation mapping
opts = MDOptions(n_iters=nΔt_coarse, n_save_iters = nΔt_coarse)
@everywhere function mutation(x)
    Xvals, _ = sample_trajectory(x, sampler, options=opts);
    return Xvals[end]
end
@everywhere mutation! = x-> sample_trajectory!(x, sampler, options=opts);

Random.seed!(100);
x0_vals = copy(voronoi_pts);
bin0_vals = bin_id.(voronoi_pts);
n_bins = length(B₀);
K̃ = WeightedEnsemble.pbuild_coarse_transition_matrix(mutation!, bin_id, x0_vals,bin0_vals, n_bins, n_samples_per_bin);

# define coarse observable as a bin function
f̃ = f.(voronoi_pts);
_,v²_vectors = WeightedEnsemble.build_coarse_vectors(n_we_steps,K̃,float.(f̃));
v² = (x,t)-> v²_vectors[t+1][bin_id(x)]
selection! = (E, B, t)-> WeightedEnsemble.optimal_allocation_selection!(E, B, v², t)


# set up ensemble
E₀ = WeightedEnsemble.Dirac_to_EnsembleWithBins(x₀, n_particles);
rebin!(E₀, B₀, 0);

# run
E = deepcopy(E₀);
B = deepcopy(B₀);
Random.seed!(200)
WeightedEnsemble.prun_we!(E, B, mutation,selection!, rebin!, n_we_steps);
p_est = f.(E.ξ) ⋅ E.ω
@printf("WE Estimate = %g\n", p_est)
