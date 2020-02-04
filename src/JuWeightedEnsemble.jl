module JuWeightedEnsemble

using LinearAlgebra
using SparseArrays
using Printf
using NearestNeighbors
using Distributed
using SharedArrays
using StatsBase
using Distributions

"""
`Ensemble{TP, TW<:AbstractFloat, TB<:Integer}`: A particle ensemble structure designed for WE

### Fields

* `ξ̂` - particle positions after selection, before mutation
* `ξ` - particle positions after mutation
* `ω̂` - partice weights after selection, before mutation
* `ω` - partice weights after mutation
* `bin` - particle bin
* `offspring` - number of off spring of the particle
"""
struct Ensemble{TP, TW<:AbstractFloat, TB<:Integer}
   # positions of the walkers after resampling, before mutation
   ξ̂::Vector{TP}
   # positions of the walkers after mutation
   ξ::Vector{TP}
   # weights of the walkers after resampling, before mutation
   ω̂::Vector{TW}
   # weights of the walkers after mutation
   ω::Vector{TW}
   # category ("bin") type of each of the walkers
   bin::Vector{TB}
   # number of offspring of each particle
   offspring::Vector{TB}
end

function Base.eltype(E::Ensemble)
   return typeof((E.ξ̂[1], E.ξ[1], E.ω̂[1], E.ω[1], E.bin[1],E.offspring[1]))
end

function Base.push!(E::Ensemble, ξ̂, ξ, ω̂, ω, bin, offspring)
   push!(E.ξ̂, ξ̂);
   push!(E.ξ, ξ);
   push!(E.ω̂, ω̂);
   push!(E.ω, ω);
   push!(E.bin, bin);
   push!(E.offspring, offspring)
end

function Base.pop!(E::Ensemble)
   ξ̂ = pop!(E.ξ̂)
   ξ = pop!(E.ξ);
   ω̂ = pop!(E.ω̂);
   ω = pop!(E.ω);
   bin = pop!(E.bin)
   offspring = pop!(E.offspring)
   return ξ̂, ξ, ω̂, ω, bin, offspring
end

function Base.popfirst!(E::Ensemble)
   ξ̂ = popfirst!(E.ξ̂)
   ξ = popfirst!(E.ξ);
   ω̂ = popfirst!(E.ω̂);
   ω = popfirst!(E.ω);
   bin = popfirst!(E.bin)
   offspring = popfirst!(E.offspring)
   return ξ̂, ξ, ω̂, ω, bin, offspring
end

function Base.length(E::Ensemble)
   return length(E.ξ)
end

function Base.isempty(E::Ensemble)
   return isempty(E.ξ)
end

function Base.iterate(E::Ensemble, state = 1)

   if state > length(E)
      return nothing
   end
   return (E.ξ̂[state],E.ξ[state], E.ω̂[state],E.ω[state], E.bin[state], E.offspring[state]), state+1
end

"""
`Bins{TS, TW<:AbstractFloat, TB<:Integer, TT<:Real}`: A bin structure designed for WE

### Fields

* `Ω` - structure containing information for uniquely identifying each bin
* `n` - number of particles in each bin
* `target` - target number of particles in each bin
* `ν` - weight of each bin
"""
struct Bins{TS, TW<:AbstractFloat, TB<:Integer, TT<:Real}
   # structure which identifies the bins
   Ω::Vector{TS}
   # number of walkers in each bin
   n::Vector{TB}
   # target number of walkers in each bin - specify if real/int with TT
   target::Vector{TT}
   # weight associated with each bin
   ν::Vector{TW}
end

function Base.push!(B::Bins, Ω, n, target, ν)
   push!(B.Ω, Ω);
   push!(B.n, n);
   push!(B.target, target);
   push!(B.ν, ν);
   #push!(B.particles, particles)
end

function Base.pop!(B::Bins)
   Ω = pop!(B.Ω);
   n = pop!(B.n);
   target = pop!(B.target);
   ν = pop!(B.ν);
   #particles = pop!(B.particles)
   return Ω, n, target, ν
end

function Base.popfirst!(B::Bins)
   Ω = popfirst!(B.Ω);
   n = popfirst!(B.n);
   target = popfirst!(B.target);
   ν = popfirst!(B.ν);
   #particles = pop!(B.particles)
   return Ω, n, target, ν
end

function Base.length(B::Bins)
   return length(B.Ω)
end

function Base.eltype(B::Bins)
   return typeof((B.Ω, B.n, B.target, B.ν))
end

function Base.iterate(B::Bins, state = 1)

   if state > length(B)
      return nothing
   end
   return (B.Ω[state], B.n[state], B.target[state], B.ν[state]), state+1
end


"""
`coarse_transition_matrix`: Contruct a transition matrix amongst the bins (serial version).

### Arguments
* `mutation!` - an in place mutation function
* `bin_id` - bin identification function
* `x0_vals` - an array of starting values
* `bin0_vals` - an array of the bins corresponding to `x0_vals`
* `n_bins` - total number of bins
* `n_samples` - number of trials for each sample
"""
function coarse_transition_matrix(mutation!, bin_id, x0_vals, bin0_vals, n_bins, n_samples)

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
`coarse_transition_matrix_parallel`: Contruct a transition matrix amongst the
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
function coarse_transition_matrix_parallel(mutation!, bin_id, x0_vals, bin0_vals, n_bins, n_samples)

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
`value_vectors`: Assemble the value vectors associated with argument `u` for `n`
steps with bin transition matrix `T`

### Arguments
* `n` - number of steps
* `T` - bin transition matrix
* `u` - quantity of interest vector on the bin space
* `tol` - tolerance for negative values which may appear due to roundoff
"""
function value_vectors(n, T, u; tol=1e-15)
   n_bins = length(u);
   vvals = zeros(n_bins,n);

   Tu = copy(u);
   v1 = similar(u);
   v2 = similar(u);

   # l = n - p - 1, p = 0,...,n-1, l = n-1,...,0
   for l in n-1:-1:0
      # at the begining of the loop Tu = T^(n-p-1) u = T^l u
      v1 = copy(Tu);
      # Tu = T^(n-p) u = T^(l+1) u
      Tu .= T * Tu;
      v2 = copy(Tu);
      @. v1 = v1^2;
      v1 .= T * v1;
      @. v2 = v2^2;
      if(minimum(v1-v2)> -tol)
         @. vvals[:,l+1] .= sqrt(max(v1 - v2, 0));
      else
         throw(DomainError(l,"Nontrivial negative entry in value vector"));
      end

   end
   return vvals
end


"""
`we_optimal_resample!`: Resample particles according to the bins, using the value vector
for optimal allocation.

### Arguments
* `E` - particle ensemble
* `B` - bin data structure
* `v` - value vector for resampling
* `resample` - resampling scheme
"""
function we_optimal_resample!(E::Ensemble, B::Bins, v, resample)

   n_particles = length(E);
   n_bins = length(B);

   # find target number of offspring for the bins based on
   # bin weights
   R = sum(B.n .>0); # count number of bins which must have offspring
   B.target .= (B.n .>0) .+ resample(n_particles-R, (B.ν .* v)/(B.ν ⋅ v));

   # compute number of offspring of each particle bin by bin
   for i in 1:n_bins
      # get particle indices for bin i
      particle_ids = findall(isequal(i), E.bin);
      if !isempty(particle_ids)
         E.offspring[particle_ids] = resample(B.target[i], E.ω[particle_ids]/B.ν[i]);
      end
   end

   # resample the particles
   n_spawned = 0;
   children_allocated = 0;
   weight_allocated = 0.0;
   for i in 1:n_particles
      # identify the bin of the current particle
      bin = E.bin[i];
      for k in 1:E.offspring[i]
         E.ξ̂[k+n_spawned] = deepcopy(E.ξ[i]);
         E.ω̂[k+n_spawned] = B.ν[bin]/B.target[bin];
      end
      n_spawned += E.offspring[i];
   end
   E, B
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
* `bin_id` - bin identification function
* `resample` - resampling scheme
* `value_vecs` - value vectors for optimal allocation
* `n_we_iters` - number of iterations
"""
function run_we!(E, B, mutation, bin_id, resample, value_vecs, n_we_iters)

   v = similar(value_vecs[:,1]);

   for j in 1:n_we_iters
      @. v = value_vecs[:,j];
      if(B.ν ⋅ v>0)
         we_optimal_resample!(E, B, v, resample);
      else
         throw(DomainError(j,"Bin Weights ⟂ Value Vector"));
      end
      @. E.ω = E.ω̂
      E.ξ .= map(mutation, E.ξ̂);
      @. E.bin = bin_id(E.ξ);
      update_bin_weights!(B, E);
   end
   E, B

end

"""
`run_we_parallel!`: Run a parallel WE simulation with value vectors for optimal
allocation.  This performs the mutation steps in parallel, and assumes a worker
pool has already been created.

### Arguments
* `E` - particle ensemble
* `B` - bin data structure
* `mutation` - mutation function
* `bin_id` - bin identification function
* `resample` - resampling scheme
* `value_vecs` - value vectors for optimal allocation
* `n_we_iters` - number of iterations
"""
function run_we_parallel!(E, B, mutation, bin_id, resample, value_vecs, n_we_iters)

   v = similar(value_vecs[:,1]);

   for j in 1:n_we_iters
      @. v = value_vecs[:,j];
      if(B.ν ⋅ v>0)
         we_optimal_resample!(E, B, v, resample);
      else
         throw(DomainError(j,"Bin Weights ⟂ Value Vector"));
      end
      @. E.ω = E.ω̂
      E.ξ .= pmap(mutation, E.ξ̂);
      #E.bin .= pmap(bin_id, E.ξ);
      @. E.bin = bin_id(E.ξ);
      update_bin_weights!(B, E);
   end
   E, B

end

"""
`Residual`: perform residual sampling

### Arguments
`n` - number of trials
`ω` - probabilities
"""
function Residual(n,ω)

    m = length(ω);

    R = sum(floor.(Int, n * ω));

    if (R < n)
        ω̄ = (n * ω .- floor.(n * ω)) / (n - R);
        N̄vals = rand(Multinomial(n-R, ω̄));
        Nvals = floor.(Int, n * ω) .+ N̄vals;
    else
        Nvals = floor.(Int, n * ω);
    end

    return Nvals

end

"""
`Stratified`: perform stratified sampling

### Arguments
`n` - number of trials
`ω` - probabilities
"""
function Stratified(n,ω)
    U = range(0,stop=n-1)/n .+ rand(n)/n;

    Nvals = counts(quantile.(Categorical(ω), U), 1:length(ω));

    return Nvals
end

"""
`Systematic`: perform stratified sampling

### Arguments
`n` - number of trials
`ω` - probabilities
"""
function Systematic(n,ω)

    U = range(0,stop=n-1)/n .+ rand()/n;

    Nvals = counts(quantile.(Categorical(ω), U), 1:length(ω));

    return Nvals
end


export Ensemble, Bins
end # end module
