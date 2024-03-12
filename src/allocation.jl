
"""
    trivial_allocation!(E, B)

Trivially allocate each particle to have one offspring.
### Arguments
* `E` - particle ensemble
* `B` - bin data structure
"""
function trivial_allocation!(E::TE, B::TB) where {TE<:Ensemble,TB<:Bins}
   @. E.o = 1
   @. B.target = B.n
   E, B
end

"""
    minimal_bin_allocation!(B; νmin=νmin)

Allocates a single particle to be spawned within each nonempty bin and the
current number of particles in any bin with less than νmin total mass.
### Arguments
* `B` - bin data structure
### Optional Arguments
* `νmin=νmin` - minimal bin weight to avoid underflow
"""
function minimal_bin_allocation!(B::TB; νmin=νmin) where {TB<:Bins}

   n_bins = length(B)

   # increment the number of offspring by one for each nontrivial bin
   nontrivial_bins = findall(ν -> ν ≥ νmin, B.ν)
   @. B.target[nontrivial_bins] = 1
   # set the number of offspring to be the same for any trivial (small mass) bins
   trivial_bins = setdiff(1:n_bins, nontrivial_bins)
   @. B.target[trivial_bins] = B.n[trivial_bins]

   B
end

"""
    static_bin_allocation!(B, n_static; ωmin=ωmin)

Allocates a predeterimined number of particles to be spawned within each
nonempty bin and the current number of particles in any bin, provided no
particle would be given less than a specified minimal mass.

### Arguments
* `B` - bin data structure
* `n_static` - array of predetermined bin allocation numbers
### Optional Arguments
* `ωmin=ωmin` - minimal particle weight to avoid underflow
"""
function static_bin_allocation!(B::TB, n_static::Vector{Int}; ωmin=ωmin) where {TB<:Bins}

   # set the number of offspring according to the static allocation, subject to
   # the ωmin constraint
   @. B.target = Int(min(n_static, floor(B.ν / ωmin)))
   B
end

"""
    targeted_bin_allocation!(B::TB, E::TE, G::F, t::Int, n_allocate::Int; allocation_resampler=systematic, νmin=νmin) where {TE<:Ensemble,TB<:Bins,F<:Function}

Targeted allocation of particles amongst bins using a specified function, `G:(p,
E, B, t) → [0,∞)` for bin `p`.
### Arguments
* `B` - bin data structure
* `E` - particle ensemble
* `G` - target function
* `t` - t-th seletion step
* `n_allocate` - number of particles to allocate
  ### Optional Arguments
* `allocation_resampler=systematic` - resampling scheme amongst bins
* `νmin=νmin` - minimal bin weight to avoid underflow
"""
function targeted_bin_allocation!(B::TB, E::TE, G::F, t::Int, n_allocate::Int; allocation_resampler=systematic, νmin=νmin) where {TE<:Ensemble,TB<:Bins,F<:Function}

   n_bins = length(B)

   # identify bins with nontrivial amount of mass
   nontrivial_bins = findall(ν -> ν ≥ νmin, B.ν)

   Ñ = zeros(n_bins)
   ρ = zeros(n_bins)

   for p in nontrivial_bins
      Ñ[p] = G(p, E, B, t)
   end

   # compute probabilities when this can be normalized
   ρ .= Ñ / sum(Ñ)
   B.target .+= allocation_resampler(n_allocate, ρ)

   B
end

"""
    optimal_bin_allocation!(B, E, v², t, n_allocate; allocation_resampler=systematic, νmin=νmin)

Optimally particles according to the bins, using a value function, `v²` to approximate
mutation variance.
### Arguments
* `B` - bin data structure
* `E` - particle ensemble
* `v²` - v² variance function estimator
* `t` - t-th seletion step
* `n_allocate` - number of particles to allocate 
### Optional Arguments
* `allocation_resampler=systematic` - resampling scheme amongst bins
* `νmin=νmin` - minimal bin weight to avoid underflow
"""
function optimal_bin_allocation!(B::TB, E::TE, v²::F, t::Int, n_allocate::Int; allocation_resampler=systematic, νmin=νmin) where {TE<:Ensemble,TB<:Bins,F<:Function}

   function G(p, E, B, t)
      particle_ids = findall(isequal(p), E.b)
      bin_target = sqrt(B.ν[p] * (E.ω[particle_ids] ⋅ v².(E.ξ[particle_ids], t)))
      return bin_target
   end

   targeted_bin_allocation!(B, E, G, t, n_allocate, allocation_resampler=allocation_resampler, νmin=νmin)

   B
end

"""
    uniform_bin_allocation!(B, E, n_allocate; allocation_resampler=systematic, νmin=νmin)

Uniformly allocate particles amongst bins.
### Arguments
* `B` - bin data structure
* `E` - particle ensemble
* `n_allocate` - number of particles to allocate
### Optional Arguments
* `allocation_resampler=systematic` - resampling scheme amongst bins
* `νmin=νmin` - minimal bin weight to avoid underflow
"""
function uniform_bin_allocation!(B::TB, E::TE, n_allocate::Int; allocation_resampler=systematic, νmin=νmin) where {TE<:Ensemble,TB<:Bins}
   G = (p, E, B, t) -> 1.0
   targeted_bin_allocation!(B, E, G, 0, n_allocate, allocation_resampler=allocation_resampler, νmin=νmin)
   B
end

"""
    within_bin_allocation!(E, B; within_bin_resampler=multinomial)

Once the number of offspring within each bin are set, allocate them amongst the
particles within the bin.  This assumes that the bin allocations of the bins
have completed.
### Arguments
* `E` - particle ensemble
* `B` - bin data structure
### Optional Arguments
* `within_bin_resampler=multinomial` - resampling scheme within bins
"""
function within_bin_allocation!(E::TE, B::TB; within_bin_resampler=multinomial) where {TE<:Ensemble,TB<:Bins}

   # identify bins that have a nonzero number of offspring
   non_empty_bins = findall(j -> j > 0, B.target)

   # compute number of offspring of each particle bin by bin for bins that will have offspring
   for p in non_empty_bins
      # get particle indices for bin p
      particle_ids = findall(isequal(p), E.b)
      @inbounds E.o[particle_ids] .+= within_bin_resampler(B.target[p], E.ω[particle_ids] / B.ν[p])
   end
   E
end