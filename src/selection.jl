"""
`trivial_selection!`: Trivial selection, copying over particles

### Arguments
* `E` - particle ensemble
"""
function trivial_selection!(E::TE) where {TE<:Ensemble}
    @. E.o = 1
    @. E.ω̂ = E.ω
    @. E.b̂ = E.b
    @. E.ξ̂ = deepcopy(E.ξ)
    @. E.d̂ = deepcopy(E.d)

    E
end

"""
`repopulate!`: After allocating the number of offspring of each particle, copy
the particles over.

### Arguments
* `E` - particle ensemble
* `B` - bin data structure
### Optional Arguments
* `within_bin_resampler=multinomial` - resampling scheme within bins
"""
function repopulate!(E::TE, B::TB; within_bin_resampler = multinomial) where {TE<:Ensemble,TB<:Bins}
    n_particles = length(E)
    non_empty_bins = findall(n -> n > 0, B.n)

    # compute number of offspring of each particle bin by bin
    for p in non_empty_bins
        # get particle indices for bin p
        particle_ids = findall(isequal(p), E.b)
        @inbounds E.o[particle_ids] .+= within_bin_resampler(B.target[p], E.ω[particle_ids] / B.ν[p])
    end

    n_spawned = 0
    # copy over the particles
    for i = 1:n_particles
        # identify the bin of the current particle
        @inbounds bin = E.b[i]
        for k = 1:E.o[i]
            @inbounds E.ξ̂[k+n_spawned] = deepcopy(E.ξ[i])
            @inbounds E.ω̂[k+n_spawned] = B.ν[bin] / B.target[bin]
            @inbounds E.b̂[k+n_spawned] = bin
            @inbounds E.d̂[k+n_spawned] = deepcopy(E.d[i])
        end
        @inbounds n_spawned += E.o[i]
    end
    E, B
end

"""
`uniform_selection!`: Uniformly select particles, ensuring each bin with
positive bin weight has at least one offspring.

### Arguments
* `E` - particle ensemble
* `B` - bin data structure
### Optional Arguments
* `allocation_resampler=systematic` - resampling scheme amongst bins
* `within_bin_resampler=multinomial` - resampling scheme within bins
"""
function uniform_selection!(E::TE, B::TB; allocation_resampler = systematic, within_bin_resampler = multinomial) where {TE<:Ensemble,TB<:Bins}

   # zero out offspring counts
   @. E.o = 0
   @. B.target = 0
   # ensure each nonempty bin has at least one particle
   minimal_bin_allocation!(B)
   n_particles = length(E)
   # number of remaining particles to allocate
   n_allocate = n_particles - sum(B.target)
   # allocate remaining particles
   uniform_allocation!(B, E, n_allocate, allocation_resampler = allocation_resampler)
   # populate the particles
   repopulate!(E, B, within_bin_resampler = within_bin_resampler)

   E, B
end

"""
`optimal_selection!`: Perform optimal selection of the particles, ensuring each
non empty bin has at least one particle.

### Arguments
* `E` - particle ensemble
* `B` - bin data structure
* `v²` - v² variance function estimator
* `t` - t-th seletion step
### Optional Arguments
* `allocation_resampler=systematic` - resampling scheme amongst bins
* `within_bin_resampler=multinomial` - resampling scheme within bins
"""
function optimal_selection!(E::TE, B::TB, v²::F, t::Int; allocation_resampler = systematic, within_bin_resampler = multinomial) where {TE<:Ensemble,TB<:Bins,F<:Function}

   # zero out offspring counts
   @. E.o = 0
   @. B.target = 0
   # ensure each nonempty bin has at least one particle
   minimal_bin_allocation!(B)
   n_particles = length(E)
   # number of remaining particles to allocate
   n_allocate = n_particles - sum(B.target)
   # allocate remaining particles
   optimal_allocation!(B, E, v², t, n_allocate, allocation_resampler = allocation_resampler)
   # populate the particles
   repopulate!(E, B, within_bin_resampler = within_bin_resampler)

   E, B
end

"""
`targeted_selection!`: Perform targeted selection of the particles, ensuring each
non empty bin has at least one particle.

### Arguments
* `E` - particle ensemble
* `B` - bin data structure
* `G` - target function
* `t` - t-th seletion step
### Optional Arguments
* `allocation_resampler=systematic` - resampling scheme amongst bins
* `within_bin_resampler=multinomial` - resampling scheme within bins
"""
function targeted_selection!(E::TE, B::TB, G::F, t::Int; allocation_resampler = systematic, within_bin_resampler = multinomial) where {TE<:Ensemble,TB<:Bins,F<:Function}

   # zero out offspring counts
   @. E.o = 0
   @. B.target = 0
   # ensure each nonempty bin has at least one particle
   minimal_bin_allocation!(B)
   n_particles = length(E)
   # number of remaining particles to allocate
   n_allocate = n_particles - sum(B.target)
   # allocate remaining particles
   targeted_allocation!(B, E, G, t, n_allocate, allocation_resampler = allocation_resampler)
   # populate the particles
   repopulate!(E, B, within_bin_resampler = within_bin_resampler)

   E, B
end