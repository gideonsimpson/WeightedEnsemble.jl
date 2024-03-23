voronoi_centers = [[x_] for x_ in LinRange(-1.0, 1.0, 11)];
B, bin_id, rebin! = setup_Voronoi_bins(voronoi_centers);

n_particles = 100;
x = [[x_] for x_ in LinRange(-1, 1, n_particles ÷ 2)];
append!(x, [[-1.05] for _ in 1:n_particles÷2])
E = Ensemble(x);
rebin!(E, B, 0);

G(p, E, B, t) = 2 - abs(B.Ω[p][1]);

Random.seed!(100);
targeted_selection!(E, B, G, 0);
count(E.b̂ .== 6) == 12