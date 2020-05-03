#=
Setup for estimation of the probability for the Schlögl chemical reaction network
with X(0)= x₀ to reach X(T) ∈ A, the target set.

This reaction network corresponds to:
A + 2S → 3S
3S → A + 2S
B → S
S → B

Parameters were chosen so that the system is bistable, as in:
https://doi.org/10.1063/1.5017955
=#

using LinearAlgebra
using Random
using DiffEqBiological

V = 25
a = 1
b = 2
c1 = 3 * 2
c2 = 0.6 * 6
c3 = 0.25
c4 = 2.95

params = (c1, c2, c3, c4, a, b, V);
Tmax = 1.0;

schlogl = @reaction_network begin
    c1 * a/ V, 2*S --> 3*S
    c2 / V^2, 3*S --> 2*S
    c3 * b * V, ∅ --> S
    c4, S --> ∅
end c1 c2 c3 c4 a b V

# preallocate data structure
dprob = DiscreteProblem(schlogl, [1], (0.0, 1.0), params);

# define rxn integrator
function integrate_schlogl(x, t)
    jprob = JumpProblem(remake(dprob, u0=x, tspan=(0,t)), Direct(),
        schlogl, save_positions = (false,false));
    jsol = solve(jprob, SSAStepper());
    return jsol.u[end];
end


x₀ = [5];       # initial state
A = (80,100);   # target set, (80, 100)
f = x-> Int(A[1]<x[1]<A[2]); # indicator function for A
