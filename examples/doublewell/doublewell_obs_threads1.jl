using Statistics
using HypothesisTests
using Printf
using WeightedEnsemble

include("doublewell_setup.jl");

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
B₀, bin_id, rebin! = setup_Voronoi_bins(voronoi_pts);

# define the mutation mapping
opts = MDOptions(n_iters=nΔt_coarse, n_save_iters = nΔt_coarse)
mutation! = x-> sample_trajectory!(x, sampler, options=opts);

selection! = (E, B, t)-> WeightedEnsemble.uniform_allocation!(E, B);

f = x-> Float64(-0.1 < x[1] < 0.1); # define observable
f1 = x-> Float64(-0.1 < x[1] < 0.1); # define observables
f2 = x-> Float64( x[1] > 1.1);
obs = (f1, f2);

# set up ensemble
E₀ = WeightedEnsemble.Dirac_to_EnsembleWithBins(x₀, n_particles);
rebin!(E₀, B₀, 0);

Random.seed!(200)
f_trajectory = WeightedEnsemble.trun_we_observables(E₀, B₀, mutation!, selection!, rebin!, n_we_steps, (f,));
