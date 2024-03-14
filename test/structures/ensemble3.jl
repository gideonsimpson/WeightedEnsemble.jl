x3 = [[x_] for x_ in LinRange(-1,1,101)];
n_particles = length(x3);
E3 = Ensemble(x3, 1.0/n_particles * ones(n_particles) );
f3(x) = x[1]^2;

E3.ω ⋅ f3.(E3.ξ) ≈ .34