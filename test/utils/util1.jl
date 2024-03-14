n_particles = 100;
x1 = [1.0];

E1 = Dirac_to_Ensemble(x1, n_particles);
f1(x) = x[1];
E1.ω ⋅ f1.(E1.ξ) ≈ 1