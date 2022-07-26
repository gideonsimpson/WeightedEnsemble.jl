#=
Common parameters for estimating the probability of for a diffusion with X(0) = x₀
satisfying X(T) ∈ B for the Muller potential.
=#

using LinearAlgebra
using Random
using BasicMD
using ForwardDiff
using Optim
using TestLandscapes: Muller

#  set up potential
V = x-> Muller(x)
# commented out due to threading issue with ForwardDiff
# cfg = ForwardDiff.GradientConfig(V, zeros(Float64,2));
∇V! = (gradV, X) -> ForwardDiff.gradient!(gradV, V, X);

# Identify minima for defining x₀ and target set B
mins=[optimize(V, x) for x in [[-0.5, 1.5],[0., 0.5],[0.5, 0.]]];
x₀ = copy(mins[1].minimizer);   # initial point
xt = copy(mins[3].minimizer);   # point defining target set B

T = 10;     # terminal time
Δt = 1e-4;  # time step
nΔt = Int(T/Δt);
β = 0.1;    # inverse temperature

# define observable as the indicator function on B = {x∣ |x-xt|≤r}
r = 0.5;
f = X-> Int(norm(X.-xt)≤r);

sampler = MALA(V,∇V!, β, Δt);
