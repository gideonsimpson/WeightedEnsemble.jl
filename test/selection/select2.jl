voronoi_centers = [[x_] for x_ in LinRange(-1.0, 1.0, 11)];
B, bin_id, rebin! = setup_Voronoi_bins(voronoi_centers);

n_particles = 100;
x = [[x_] for x_ in LinRange(-1, 1, n_particles÷2)];
append!(x, [[0.05] for _ in 1:n_particles÷2])
E = Ensemble(x);

rebin!(E, B, 0);

n_static = 3 .+ (1:1:11)|>collect;
n_static[end] = 15;

Random.seed!(100);
static_selection!(E, B, n_static, 0);
count(E.b̂ .== 10) == 13