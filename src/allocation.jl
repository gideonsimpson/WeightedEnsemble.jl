"""
`minimal_bin_allocation!`: Allocates a single particle to be spawned within each
nonempty bin

### Arguments
* `B` - bin data structure
"""
function minimal_bin_allocation!(B::TB) where {TB<:Bins}

   # increment the number of offspring by one for each nonempty bin
   @. B.target += (B.n > 0)
   B
end

"""
`targeted_allocation!`: Targeted allocation of particles using a specified
function, `G:(p, E, B, t) → [0,∞)` for bin `p`.  Empty bins are allocated zero
children and nonempty bins are allocated ≥ 1 child.

### Arguments
* `E` - particle ensemble
* `B` - bin data structure
* `G` - target function
* `t` - t-th seletion step
* `n_allocate` - number of particles to allocate
### Optional Arguments
* `allocation_resampler=systematic` - resampling scheme amongst bins
"""
function targeted_allocation!(E::TE, B::TB, G::F, t::Int, n_allocate::Int; allocation_resampler = systematic) where {TE<:Ensemble,TB<:Bins,F<:Function}

   # identify nonempty bins
   non_empty_bins = findall(n -> n > 0, B.n)
   Ñ = zeros(n_bins)
   ρ = zeros(n_bins)

   for p in non_empty_bins
      Ñ[p] = G(p, E, B, t)
   end

   if (sum(Ñ) > 0)
      # compute probabilities when this can be normalized
      ρ .= Ñ / sum(Ñ)
   else
      # uniformly allocate amongst the bins 
      ρ .= non_empty_bins / sum(non_empty_bins)
   end
   B.target .+= allocation_resampler(n_allocate, ρ)

   E, B
end

"""
`optimal_allocation!`: Optimally particles according to the bins,
using a value function to approximate mutation variance.

### Arguments
* `E` - particle ensemble
* `B` - bin data structure
* `v²` - v² variance function estimator
* `t` - t-th seletion step
* `n_allocate` - number of particles to allocate 
### Optional Arguments
* `allocation_resampler=systematic` - resampling scheme amongst bins
"""
function optimal_allocation!(E::TE, B::TB, v²::F, t::Int, n_allocate::Int; allocation_resampler = systematic) where {TE<:Ensemble,TB<:Bins,F<:Function}

   function G(p, E, B, t)
      particle_ids = findall(isequal(p), E.b)
      bin_target = sqrt(B.ν[p] * (E.ω[particle_ids] ⋅ v².(E.ξ[particle_ids], t)))
      return bin_target
   end

   targeted_allocation!(E, B, G, t, n_allocate, allocation_resampler = allocation_resampler)

   E, B
end

"""
`uniform_allocation!`: Uniformly allocate particles, ensuring each bin with
positive bin weight has at least one offspring.

### Arguments
* `E` - particle ensemble
* `B` - bin data structure
* `n_allocate` - number of particles to allocate
### Optional Arguments
* `allocation_resampler=systematic` - resampling scheme amongst bins
"""
function uniform_allocation!(E::TE, B::TB, n_allocate::Int; allocation_resampler = systematic) where {TE<:Ensemble,TB<:Bins}
   G = (p, E, B, t) -> 1.0
   targeted_allocation!(E, B, G, 0, n_allocate, allocation_resampler = allocation_resampler)
   E, B
end

