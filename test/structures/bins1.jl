voronoi_centers = [[x_] for x_ in LinRange(-1., 1., 11)];
B, bin_id, rebin! = setup_Voronoi_bins(voronoi_centers);
bin_id([-0.78]) == 2