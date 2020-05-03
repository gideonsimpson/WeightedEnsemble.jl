
using Distributed
# nw = 4; # number of workers
# addprocs(nw);

using StatsBase
using HypothesisTests
using Printf
using ProgressMeter
@everywhere include("schlogl_setup.jl")

Random.seed!(100)
@everywhere n_trials = 10^6;
n_in_set = @sync @showprogress @distributed (+) for j in 1:n_trials
    X₀ = copy(x₀);
    f(integrate_schlogl(X₀, Tmax));
end
CI = confint(BinomialTest(n_in_set, n_trials));
@printf("P Estimate = %g\n",n_in_set/n_trials);
@printf("95%% CI = (%g,%g)\n", CI[1], CI[2])
