#=
WE estimation of the probability of for a diffusion with X(0) = x₀ satisfying
X(T) ∈ B for the Muller potential.

Tests uniform allocation
=#

include("muller_setup.jl");

# number of coarse steps in WE
n_we_steps = 10;
# number of time steps during mutation step
nΔt_coarse = Int(nΔt / n_we_steps);
# number of samples in coarse matrix
n_samples_per_bin = 10^2;
# ensemble size
n_particles = 10^2;

# define bin structure using Voronoi
xc = LinRange(-1.5, 1, 7)
yc = LinRange(-0.5, 2, 7)
voronoi_pts = Array{Float64,1}[];
for x in xc, y in yc
    # only include points that are likely to be accessed
    if (Muller([x, y]) < 250)
        push!(voronoi_pts, [x, y])
    end
end
B₀, bin_id, rebin! = setup_Voronoi_bins(voronoi_pts);

opts = MDOptions(n_iters = nΔt_coarse, n_save_iters = nΔt_coarse)
mutation! = x -> sample_trajectory!(x, sampler, options = opts);

# define selection function
selection! = (E, B, t)-> uniform_selection!(E, B, t)

we_sampler = WEsampler(mutation!, selection!, rebin!);

# set up ensemble
E₀ = Dirac_to_Ensemble(x₀, n_particles);
rebin!(E₀, B₀, 0);

# run
E = deepcopy(E₀);
B = deepcopy(B₀);
Random.seed!(200)
run_we!(E, B, we_sampler, n_we_steps);
p_est = f.(E.ξ) ⋅ E.ω
p_est ≈ 0.04342107610697742;
