# utility functions

"""
`update_bin_weights!`: Update the bin weights from the ensemble

### Arguments
* `B` - bin data structure
* `E` - particle ensemble
"""
function update_bin_weights!(B::TB, E::TE) where {TB<:AbstractBins, TE<:EnsembleWithBins}

   n_particles = length(E);
   n_bins = length(B);

   # this loops over bins
   for i in 1:n_bins
      particle_ids = findall(isequal(i), E.b);

      B.ν[i] = sum(E.ω[particle_ids]);
      B.n[i] = length(particle_ids);
   end
   B
end

"""
`Voronoi_to_Bins`: Convenience function for constructing a bin structure using a
sequence of points as the Voronoi sites

### Arguments
`sites` - An array of points defining the Voronoi cells
"""
function Voronoi_to_Bins(sites)

   B = Bins{typeof(sites[1]), Float64, Int, Int}([],[],[],[]);
   for site in sites
      push!(B, site , 0, 0, 0);
   end
   return B
end

"""
`Voronoi_bin_id`: Convenience function for bin id in Voronoi based bins

### Arguments
`X` - Point thats bin is to be determined
`tree` - A nearest neighbors tree structure constructed with `KDTree`
"""
function Voronoi_bin_id(X, tree)
   return knn(tree, X, 1)[1][1]
end

"""
`Dirac_to_EnsembleWithoutBins`: Convenience function for construction an ensemble from a
single initial walker.  This hard codes the weights to be Float64

### Arguments
`X` - Starting state of all walkers
`n_particles` - Number of walkers in the ensemble
"""
function Dirac_to_EnsembleWithoutBins(X::TP, n_particles::TI) where {TP, TI<:Integer}
   ω = 1.0/n_particles
   E = EnsembleWithoutBins{TP, Float64, TI}([deepcopy(X) for i = 1:n_particles],
                                    [deepcopy(X) for i = 1:n_particles],
                                    ω * ones(n_particles),ω * ones(n_particles),
                                    zeros(Float64, n_particles),
                                    zeros(TI, n_particles));
   return E
end


"""
`Dirac_to_EnsembleWithBins`: Convenience function for construction an ensemble from a
single initial walker.  This hard codes the weights to be Float64.

### Arguments
`X` - Starting state of all walkers
`n_particles` - Number of walkers in the ensemble
"""
function Dirac_to_EnsembleWithBins(X::TP, n_particles::TI) where {TP, TI<:Integer}
   ω = 1.0/n_particles
   E = EnsembleWithBins{TP, Float64, TI}([deepcopy(X) for i = 1:n_particles],
                                    [deepcopy(X) for i = 1:n_particles],
                                    ω * ones(n_particles),ω * ones(n_particles),
                                    zeros(TI, n_particles),zeros(TI, n_particles),
                                    zeros(TI, n_particles));
   return E
end
