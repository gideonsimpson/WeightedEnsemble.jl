#=
Common parameters for estimating the probability of for a diffusion with X(0) =
a satisfying X(T) ∈ (b, ∞) for the double well potential V(x) = (x²-1)².
=#

using LinearAlgebra
using Random
using BasicMD

function V(X)
    return (X⋅X - 1)^2;
end
function gradV!(gradV, X)
    @. gradV = 4.0 * (X^2 - 1) * X;
    gradV
end

a = -1.0;   # starting point
x₀ = [a];
b = 0.5;    # target set, (b, ∞)
Δt = 1e-4;  # time step
β = 5.0;    # inverse temperature
T = 1.0;    # terminal time
nΔt = ceil(Int, T/Δt);

f = x-> Int(x[1]>b); # define observable

# define sampler
sampler = MALA(V, gradV!, β, Δt);
