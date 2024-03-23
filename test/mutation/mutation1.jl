
voronoi_centers = [[x_] for x_ in LinRange(-1.0, 1.0, 11)];
B, bin_id, rebin! = setup_Voronoi_bins(voronoi_centers);

n_particles = 100;
x = [[x_] for x_ in LinRange(-1, 1, n_particles)];
E = Ensemble(x);
rebin!(E, B, 0);

function mutation1!(x)
    x[1] +=1.
    x
end

mutation1!.(E.ξ);
rebin!(E, B, 1);
B.n[1] == 0
