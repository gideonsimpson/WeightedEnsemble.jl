#=
Direct estimation of the probability of for a diffusion with X(0) = x₀ satisfying
X(T) ∈ B for the Muller potential.
=#

using Distributed
# nw = 4; # number of workers
# addprocs(nw);

using Statistics
using HypothesisTests
using Printf
@everywhere include("muller_setup.jl")
opts = MDOptions(n_iters=nΔt)

Random.seed!(100)
@everywhere n_trials = 10^4;
n_in_set = @sync @distributed (+) for j in 1:n_trials
    X = copy(x₀);
    sample_trajectory!(X, sampler, options=opts);
    f(X)
end
CI = confint(BinomialTest(n_in_set, n_trials));
@printf("95%% CI = (%g,%g)\n", CI[1], CI[2])
