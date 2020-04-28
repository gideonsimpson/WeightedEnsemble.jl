#=
Direct estimation of the probability of for a diffusion with X(0) = x₀ satisfying
X(T) ∈ B for the Muller potential.
=#

using StatsBase
using HypothesisTests
using Printf

include("muller_setup.jl");
opts = MDOptions(n_iters=nΔt)

Random.seed!(100);
n_trials = 10^4;
n_in_set = 0;
for j in 1:n_trials
    X = copy(x₀);
    sample_trajectory!(X, sampler, options=opts);
    global n_in_set+= f(X);
end
CI = confint(BinomialTest(n_in_set, n_trials));
@printf("95%% CI = (%g,%g)\n", CI[1], CI[2])
