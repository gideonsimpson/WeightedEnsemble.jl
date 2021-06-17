"""
`targeted_allocation!`: Targeted allocation of particles using a
specified function, `g`, such that the party allocation is proportional to `g ⋅
ω` for each bin.

### Arguments
* `E` - particle ensemble
* `B` - bin data structure
* `g` - target function
* `t` - t-th seletion step
* `resample` - resampling scheme
"""
function targeted_allocation!(E::TE, B::TB, g::F, t::Int; resample=Systematic) where{TE<:Ensemble, TB<:Bins, F<:Function}

   n_particles = length(E);
   n_bins = length(B);
   # zero out offspring counts
   @. E.o = 0;
   @. B.target = 0;

   # identify nonempty bins
   non_empty_bins = findall(n->n>0, B.n);
   R = length(non_empty_bins);
   Ñ = zeros(n_bins);
   ρ = zeros(n_bins);

   for p in non_empty_bins
      particle_ids = findall(isequal(p), E.b);
      @inbounds Ñ[p] = E.ω[particle_ids] ⋅ g.(E.ξ[particle_ids],t);
   end

   if(sum(Ñ)>0)
      # compute probabilities
      ρ .= Ñ/sum(Ñ);
      B.target .= (B.n .>0) .+ resample(n_particles-R, ρ);

      # compute number of offspring of each particle bin by bin
      for p in non_empty_bins
         # get particle indices for bin p
         particle_ids = findall(isequal(p), E.b);
         @inbounds E.o[particle_ids] .= resample(B.target[p], E.ω[particle_ids]/B.ν[p]);
      end

   else
      # every particle copies itself
      @. B.target = B.n;
      @. E.o = 1;
   end

   # resample the particles
   n_spawned = 0;
   for i in 1:n_particles
      # identify the bin of the current particle
      @inbounds bin = E.b[i];
      for k in 1:E.o[i]
         @inbounds E.ξ̂[k+n_spawned] = deepcopy(E.ξ[i]);
         @inbounds E.ω̂[k+n_spawned] = B.ν[bin]/B.target[bin];
         @inbounds E.b̂[k+n_spawned] = bin;
         @inbounds E.d̂[k+n_spawned] = deepcopy(E.d[i]);
      end
      @inbounds n_spawned += E.o[i];
   end
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
* `resample` - resampling scheme
"""
function optimal_allocation!(E::TE, B::TB, v²::F, t::Int; resample=Systematic) where{TE<:Ensemble, TB<:Bins, F<:Function}

   n_particles = length(E);
   n_bins = length(B);
   # zero out offspring counts
   @. E.o = 0;
   @. B.target = 0;

   # identify nonempty bins
   non_empty_bins = findall(n->n>0, B.n);
   R = length(non_empty_bins);
   Ñ = zeros(n_bins);
   ρ = zeros(n_bins);

   for p in non_empty_bins
      particle_ids = findall(isequal(p), E.b);
      @inbounds Ñ[p] = sqrt(B.ν[p] * (E.ω[particle_ids] ⋅ v².(E.ξ[particle_ids],t)));
   end

   if(sum(Ñ)>0)
      # compute probabilities
      ρ .= Ñ/sum(Ñ);
      B.target .= (B.n .>0) .+ resample(n_particles-R, ρ);

      # compute number of offspring of each particle bin by bin
      for p in non_empty_bins
         # get particle indices for bin p
         particle_ids = findall(isequal(p), E.b);
         @inbounds E.o[particle_ids] .= resample(B.target[p], E.ω[particle_ids]/B.ν[p]);
      end

   else
      # every particle copies itself
      @. B.target = B.n;
      @. E.o = 1;
   end

   # resample the particles
   n_spawned = 0;
   for i in 1:n_particles
      # identify the bin of the current particle
      @inbounds bin = E.b[i];
      for k in 1:E.o[i]
         @inbounds E.ξ̂[k+n_spawned] = deepcopy(E.ξ[i]);
         @inbounds E.ω̂[k+n_spawned] = B.ν[bin]/B.target[bin];
         @inbounds E.b̂[k+n_spawned] = bin;
         @inbounds E.d̂[k+n_spawned] = deepcopy(E.d[i]);         
      end
      @inbounds n_spawned += E.o[i];
   end
   E, B
end

"""
`uniform_allocation!`: Uniformly select particles, ensuring each bin with
positive bin weight has at least one offspring.

### Arguments
* `E` - particle ensemble
* `B` - bin data structure
* `resample` - resampling scheme
"""
function uniform_allocation!(E::TE, B::TB; resample=Systematic) where{TE<:Ensemble, TB<:Bins}
   n_particles = length(E);
   n_bins = length(B);
   # zero out offspring counts
   @. E.o = 0;
   @. B.target = 0;

   # ensure each bin with walkers has at least one offspring
   non_empty_bins = findall(n->n>0, B.n);
   R = length(non_empty_bins);
   @inbounds B.target[non_empty_bins] .= 1 .+ resample(n_particles-R, [1.0/R for j in 1:R]);

   # compute number of offspring of each particle bin by bin
   for p in non_empty_bins
      # get particle indices for bin p
      particle_ids = findall(isequal(p), E.b);
      @inbounds E.o[particle_ids] .= resample(B.target[p], E.ω[particle_ids]/B.ν[p]);
   end

   # resample the particles
   n_spawned = 0;
   for i in 1:n_particles
      # identify the bin of the current particle
      @inbounds bin = E.b[i];
      for k in 1:E.o[i]
         @inbounds E.ξ̂[k+n_spawned] = deepcopy(E.ξ[i]);
         @inbounds E.ω̂[k+n_spawned] = B.ν[bin]/B.target[bin];
         @inbounds E.b̂[k+n_spawned] = bin;
         @inbounds E.d̂[k+n_spawned] = deepcopy(E.d[i]);
      end
      @inbounds n_spawned += E.o[i];
   end
   E, B
end

function trivial_allocation!(E::TE) where{TE<:Ensemble}
   @. E.o = 1;
   @. E.ω̂ = E.ω;
   @. E.b̂ = E.b;
   @. E.ξ̂ = deepcopy(E.ξ);
   @. E.d̂ = deepcopy(E.d);

   E
end
