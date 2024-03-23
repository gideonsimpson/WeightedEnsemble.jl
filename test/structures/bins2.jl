voronoi_centers = [[x_] for x_ in LinRange(-1., 1., 11)];
B, bin_id, rebin! = setup_Voronoi_bins(voronoi_centers);

n_particles = 500;

x = [[x_] for x_ in LinRange(-1,1,n_particles)];
E = Ensemble(x);

rebin!(E, B, 0);
sum(B.ν)≈1
