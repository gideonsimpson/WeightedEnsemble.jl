let
    x0 = [[x_] for x_ in LinRange(-1,1,101)];
    n_particles = length(x0);
    E = Ensemble(x0, 1.0/n_particles * ones(n_particles) );
    f(x) = x[1]^2;

    E.ω ⋅ f.(E.ξ) ≈ .34
end