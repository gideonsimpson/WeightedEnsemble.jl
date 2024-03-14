x2 = [[x_] for x_ in LinRange(-1,1,101)];
E2 = Ensemble(x2);
f2(x) = x[1]^2;

E2.ω ⋅ f2.(E2.ξ) ≈ .34