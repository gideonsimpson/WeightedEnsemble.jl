let
    n_particles = 100;
    x0 = [1.0];
    E = Ensemble(x0, n_particles);
    f(x) = x[1];

    E.ω ⋅ f.(E.ξ) ≈ 1
end