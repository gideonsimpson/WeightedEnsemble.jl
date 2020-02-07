#=
Common parameters for estimating the probability of for a diffusion with X(0) = x₀
satisfying X(T) ∈ B for the Muller potential.
=#

using LinearAlgebra
using Random
# adjust this path as neccessary
push!(LOAD_PATH, "../../../JuBasicMD/src")
using JuBasicMD: MALA, MALA!
using JuBasicMD: EM, EM!
using ForwardDiff
using Optim

#  set up potentail
function V(x)

    aa = Float64[-1, -1, -6.5, 0.7];
    bb = Float64[0., 0., 11., 0.6];
    cc = Float64[-10., -10., -6.5, 0.7];
    AA = Float64[-200., -100., -170., 15.];
    XX = Float64[1., 0., -0.5, -1.];
    YY = Float64[0., 0.5, 1.5, 1.];

    return ( AA[1]*exp(aa[1]*(x[1]-XX[1])^2+bb[1]*(x[1]-XX[1])*(x[2]-YY[1])+cc[1]*(x[2]-YY[1])^2)
                 +AA[2]*exp(aa[2]*(x[1]-XX[2])^2+bb[2]*(x[1]-XX[2])*(x[2]-YY[2])+cc[2]*(x[2]-YY[2])^2)
                 +AA[3]*exp(aa[3]*(x[1]-XX[3])^2+bb[3]*(x[1]-XX[3])*(x[2]-YY[3])+cc[3]*(x[2]-YY[3])^2)
                 +AA[4]*exp(aa[4]*(x[1]-XX[4])^2+bb[4]*(x[1]-XX[4])*(x[2]-YY[4])+cc[4]*(x[2]-YY[4])^2));
end
cfg = ForwardDiff.GradientConfig(V, zeros(Float64,2));
∇V! = (gradV, X) -> ForwardDiff.gradient!(gradV, V, X, cfg);

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
