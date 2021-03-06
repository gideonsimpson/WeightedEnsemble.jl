#=
Parallel direct estimation of the probability for the Schlögl chemical reaction network
with X(0)= x₀ to reach X(T) ∈ A, the target set.

This reaction network corresponds to:
A + 2S → 3S
3S → A + 2S
B → S
S → B

Parameters were chosen so that the system is bistable, as in:
https://doi.org/10.1063/1.5017955
=#
using Distributed
# nw = 4; # number of workers
# addprocs(nw);

using Statistics
using HypothesisTests
using Printf

@everywhere include("schlogl_setup.jl")

Random.seed!(100)
@everywhere n_trials = 10^6;
n_in_set = @sync @distributed (+) for j in 1:n_trials
    X₀ = copy(x₀);
    f(integrate_schlogl(X₀, Tmax));
end
CI = confint(BinomialTest(n_in_set, n_trials));
@printf("P Estimate = %g\n",n_in_set/n_trials);
@printf("95%% CI = (%g,%g)\n", CI[1], CI[2])
