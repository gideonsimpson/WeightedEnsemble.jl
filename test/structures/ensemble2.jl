let
    x0 = [[x_] for x_ in LinRange(-1,1,101)];
    E = Ensemble(x0);
    f(x) = x[1]^2;

    E.ω ⋅ f.(E.ξ) ≈ .34
end