module JuWeightedEnsemble

using LinearAlgebra
using SparseArrays
using Printf
using NearestNeighbors
using Distributed
using SharedArrays
using StatsBase
using Distributions

# Load data structures
include("ensemble.jl")
include("bins.jl")

export Ensemble, Bins

# Load resampling algorithms
include("resampling.jl")
# Load selection algorithms
include("selection.jl")

"""
`build_coarse_transition_matrix`: Contruct a transition matrix amongst the bins (serial version).

### Arguments
* `mutation!` - an in place mutation function
* `bin_id` - bin identification function
* `x0_vals` - an array of starting values
* `bin0_vals` - an array of the bins corresponding to `x0_vals`
* `n_bins` - total number of bins
* `n_samples` - number of trials for each sample
"""
function build_coarse_transition_matrix(mutation!, bin_id, x0_vals, bin0_vals, n_bins, n_samples)

   X = similar(x0_vals[1]);
   T = zeros(Float64, n_bins, n_bins);

   for k in 1:n_samples, l in 1:length(x0_vals)
      X .= copy(x0_vals[l]);
      i = bin0_vals[l];
      mutation!(X);
      j = bin_id(X);
      T[i,j] +=1.0/n_samples;
   end

   return T
end

"""
`pbuild_coarse_transition_matrix`: Contruct a transition matrix amongst the
bins in parallel.  It assumed that a pool of workers has already been
contructed.

### Arguments
* `mutation!` - an in place mutation function
* `bin_id` - bin identification function
* `x0_vals` - an array of starting values
* `bin0_vals` - an array of the bins corresponding to `x0_vals`
* `n_bins` - total number of bins
* `n_samples` - number of trials for each sample
"""
function pbuild_coarse_transition_matrix(mutation!, bin_id, x0_vals, bin0_vals, n_bins, n_samples)

    T = SharedArray{Float64}(n_bins,n_bins);
    X = similar(x0_vals[1]);
    @. T  = 0.0;

    @sync @distributed for k in 1:n_samples
        for l in 1:length(x0_vals)
            X .= copy(x0_vals[l]);
            i = bin0_vals[l];
            mutation!(X);
            j = bin_id(X);
            T[i,j] +=1.0/n_samples;
        end
   end

   return T
end


"""
`update_bin_weights!`: Update the bin weights from the ensemble

### Arguments
* `B` - bin data structure
* `E` - particle ensemble
"""
function update_bin_weights!(B::Bins, E::Ensemble)

   n_particles = length(E);
   n_bins = length(B);

   # this loops over bins
   for i in 1:n_bins
      particle_ids = findall(isequal(i), E.bin);

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
`run_we!`: Run a serial WE simulation with value vectors for optimal allocation.

### Arguments
* `E` - particle ensemble
* `B` - bin data structure
* `mutation` - mutation function
* `selection!` - selection scheme
* `rebin!` - rebin and update particles and bins
* `n_we_steps` - number of steps in the WE run
"""
function run_we!(E, B, mutation, selection!, rebin!, n_we_steps)

   for t in 1:n_we_steps
      selection!(E, B, t);
      @. E.ω = E.ω̂;
      E.ξ .= map(mutation, E.ξ̂);
      rebin!(E, B, t);
   end
   E, B

end

"""
`prun_we!`: Run a parallel WE simulation with value vectors for optimal
allocation.  This performs the mutation steps in parallel, and assumes a worker
pool has already been created.

### Arguments
* `E` - particle ensemble
* `B` - bin data structure
* `mutation` - mutation function
* `selection!` - selection scheme
* `rebin!` - rebin and update particles and bins
* `n_we_steps` - number of steps in the WE run
"""
function prun_we!(E, B, mutation, selection!, rebin!, n_we_steps)

   for t in 1:n_we_steps
      selection!(E, B, t);
      @. E.ω = E.ω̂
      E.ξ .= pmap(mutation, E.ξ̂);
      rebin!(E, B, t);
   end
   E, B

end



end # end module
