#=
Direct estimation of the probability for the Schlögl chemical reaction network
with X(0)= a to reach X(T) ∈ (c,d).
=#

using StatsBase
using HypothesisTests
using Printf
using ProgressMeter

include("schlogl_setup.jl");

Random.seed!(100);
n_trials = 10^6;
n_in_set = 0;
@showprogress for j in 1:n_trials
    X₀ = copy(x₀);
    global n_in_set+= f(integrate_schlogl(X₀, Tmax));
end
CI = confint(BinomialTest(n_in_set, n_trials));
@printf("P Estimate = %g\n",n_in_set/n_trials);
@printf("95%% CI = (%g,%g)\n", CI[1], CI[2]);
