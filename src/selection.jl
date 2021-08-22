"""
`targeted_allocation!`: Targeted allocation of particles using a specified
function, `G:(p, E, B, t) → [0,∞)` for bin `p`.  Empty bins are allocated zero
children and nonempty bins are allocated ≥ 1 child.

### Arguments
* `E` - particle ensemble
* `B` - bin data structure
* `G` - target function
* `t` - t-th seletion step
  ### Optional Arguments
* `allocation_resampler=systematic` - resampling scheme amongst bins
* `within_bin_resampler=multinomial` - resampling scheme within bins
"""
function targeted_allocation!(E::TE, B::TB, G::F, t::Int; allocation_resampler=systematic, within_bin_resampler=multinomial) where{TE<:Ensemble, TB<:Bins, F<:Function}

   n_particles = length(E);
   n_bins = length(B);
   # zero out offspring counts
   @. E.o = 0;
   @. B.target = 0;

   # identify nonempty bins
   non_empty_bins = findall(n->n>0, B.n);
   n_occupied = length(non_empty_bins);
   Ñ = zeros(n_bins);
   ρ = zeros(n_bins);

   for p in non_empty_bins
      Ñ[p] = G(p, E, B, t);
   end

   if(sum(Ñ)>0)
      # compute probabilities
      ρ .= Ñ/sum(Ñ);
      B.target .= (B.n .>0) .+ allocation_resampler(n_particles-n_occupied, ρ);

      # compute number of offspring of each particle bin by bin
      for p in non_empty_bins
         # get particle indices for bin p
         particle_ids = findall(isequal(p), E.b);
         @inbounds E.o[particle_ids] .= within_bin_resampler(B.target[p], E.ω[particle_ids]/B.ν[p]);
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
### Optional Arguments
* `allocation_resampler=systematic` - resampling scheme amongst bins
* `within_bin_resampler=multinomial` - resampling scheme within bins
"""
function optimal_allocation!(E::TE, B::TB, v²::F, t::Int; allocation_resampler=systematic, within_bin_resampler=multinomial) where{TE<:Ensemble, TB<:Bins, F<:Function}

   function G(p, E, B, t)
      particle_ids = findall(isequal(p), E.b);
      bin_target = sqrt(B.ν[p] * (E.ω[particle_ids] ⋅ v².(E.ξ[particle_ids],t)));
      return bin_target
   end

   targeted_allocation!(E, B, G, t, allocation_resampler=allocation_resampler, within_bin_resampler=within_bin_resampler);
   
   E, B
end

"""
`uniform_allocation!`: Uniformly select particles, ensuring each bin with
positive bin weight has at least one offspring.

### Arguments
* `E` - particle ensemble
* `B` - bin data structure
### Optional Arguments
* `allocation_resampler=systematic` - resampling scheme amongst bins
* `within_bin_resampler=multinomial` - resampling scheme within bins
"""
function uniform_allocation!(E::TE, B::TB; allocation_resampler=systematic, within_bin_resampler=multinomial) where{TE<:Ensemble, TB<:Bins}
   G = (p, E, B, t)-> 1.;
   targeted_allocation!(E, B, G, 0, allocation_resampler=allocation_resampler, within_bin_resampler=within_bin_resampler);
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
