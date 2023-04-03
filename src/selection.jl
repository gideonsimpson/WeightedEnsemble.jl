"""
`trivial_selection!`: Trivial selection, copying over particles

### Arguments
* `E` - particle ensemble
"""
function trivial_selection!(E::TE) where {TE<:Ensemble}
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
"""
function repopulate!(E::TE, B::TB) where {TE<:Ensemble,TB<:Bins}
    
    # number of allocated particles <= number of particles
    n_particles = length(E)
    n_allocated = sum(B.target); 
    n_zero_mass = n_particles - n_allocated;
    n_spawned = 0
    # copy over the particles allocated by the bin allocation
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

    # zero mass particles are allocated to maintain a fixed total particle
    # count.  For simplicity, just copy over particles 1,...,n_zero_mass, but
    # given them zero weight.
    if(n_zero_mass>0)
        @printf(" ZERO MASS ALLOCATION, %d\n", n_zero_mass)
    end

    for i in 1:n_zero_mass
        @inbounds E.ξ̂[i+n_allocated] = deepcopy(E.ξ[i])
        @inbounds E.ω̂[i+n_allocated] = 0;
        @inbounds E.b̂[i+n_allocated] = E.b[i];
        @inbounds E.d̂[i+n_allocated] = deepcopy(E.d[i])
    end

    E, B
end

"""
`uniform_selection!`: Uniformly select particles, ensuring each bin with
positive bin weight has at least one offspring.

### Arguments
* `E` - particle ensemble
* `B` - bin data structure
* `t` - t-th seletion step
### Optional Arguments
* `allocation_resampler=systematic` - resampling scheme amongst bins
* `within_bin_resampler=multinomial` - resampling scheme within bins
"""
function uniform_selection!(E::TE, B::TB, t::Int; allocation_resampler = systematic, within_bin_resampler = multinomial, νmin=νmin) where {TE<:Ensemble,TB<:Bins}

    # zero out offspring counts
    @. E.o = 0
    @. B.target = 0
    # ensure each nonempty bin has at least one particle
    minimal_bin_allocation!(B, νmin=νmin)
    n_particles = length(E)
    # number of remaining particles to allocate
    n_allocate = n_particles - sum(B.target)

    try
        # allocate remaining particles
        uniform_bin_allocation!(B, E, n_allocate, allocation_resampler = allocation_resampler,νmin=νmin)
        # set number of offspring of each particle
        within_bin_allocation!(E, B, within_bin_resampler = within_bin_resampler)
    catch e
        # fall back to trivial allocation if uniform fails
        if e isa DomainError
            @printf("[%d]: TRIVIAL ALLOCATION\n", t)
            trivial_allocation!(E, B)
        else
            rethrow()
        end
    end
    # populate the particles
    repopulate!(E, B)

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
function optimal_selection!(E::TE, B::TB, v²::F, t::Int; allocation_resampler = systematic, within_bin_resampler = multinomial,νmin=νmin) where {TE<:Ensemble,TB<:Bins,F<:Function}

    # zero out offspring counts
    @. E.o = 0
    @. B.target = 0
    # ensure each nonempty bin has at least one particle
    minimal_bin_allocation!(B, νmin=νmin)
    n_particles = length(E)
    # number of remaining particles to allocate
    n_allocate = n_particles - sum(B.target)
    try
        # allocate remaining particles
        optimal_bin_allocation!(B, E, v², t, n_allocate, allocation_resampler = allocation_resampler,νmin=νmin);
        # set number of offspring of each particle
        within_bin_allocation!(E, B, within_bin_resampler = within_bin_resampler);
    catch e
        # fall back to trivial allocation if optimal fails
        if e isa DomainError
            @printf("[%d]: TRIVIAL ALLOCATION\n", t)
            trivial_allocation!(E, B)
        else
            rethrow()
        end
    end
    # populate the particles
    repopulate!(E, B);

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
function targeted_selection!(E::TE, B::TB, G::F, t::Int; allocation_resampler = systematic, within_bin_resampler = multinomial,νmin=νmin) where {TE<:Ensemble,TB<:Bins,F<:Function}

    # zero out offspring counts
    @. E.o = 0
    @. B.target = 0
    # ensure each nonempty bin has at least one particle
    minimal_bin_allocation!(B, νmin=νmin)
    n_particles = length(E)
    # number of remaining particles to allocate
    n_allocate = n_particles - sum(B.target)
    try
        # allocate remaining particles
        targeted_bin_allocation!(B, E, G, t, n_allocate, allocation_resampler = allocation_resampler,νmin=νmin)
        # set number of offspring of each particle
        within_bin_allocation!(E, B, within_bin_resampler=within_bin_resampler)
    catch e
        # fall back to trivial allocation if targeted fails        
        if e isa DomainError
            @printf("[%d]: TRIVIAL ALLOCATION\n",t);
            trivial_allocation!(E, B)
        else
            rethrow()
        end
    end
    # populate the particles
    repopulate!(E, B)

    E, B
end

"""
`static_selection!`: Select particles according to a static allocation rule.  If
the number of allocated particles is < N, the remaining particles are assigned
zero weight.

### Arguments
* `E` - particle ensemble
* `B` - bin data structure
* `static_allocate` - array of predetermined bin allocation numbers
### Optional Arguments
* `within_bin_resampler=multinomial` - resampling scheme within bins
"""
function static_selection!(E::TE, B::TB, static_allocate::Vector{Int}; within_bin_resampler=multinomial, ωmin=ωmin) where {TE<:Ensemble,TB<:Bins}

    # zero out offspring counts
    @. E.o = 0
    @. B.target = 0
    try
        static_bin_allocation!(B, static_allocate, ωmin=ωmin)
        within_bin_allocation!(E, B, within_bin_resampler=within_bin_resampler)
    catch e
        # fall back to trivial allocation if uniform fails
        if e isa DomainError
            @printf("[%d]: TRIVIAL ALLOCATION\n", t)
            trivial_allocation!(E, B)
        else
            rethrow()
        end
    end
    # populate the particles
    repopulate!(E, B)

    E, B
end