"""
`optimal_allocation_selection!`: Optimally particles according to the bins,
using a value function to approximate mutation variance.

### Arguments
* `E` - particle ensemble
* `B` - bin data structure
* `v²` - v² variance function estimator
* `t` - t-th seletion step
* `resample` - resampling scheme
"""
function optimal_allocation_selection!(E::TE, B::TB, v²::F, t::Int; resample=Systematic) where{TE<:EnsembleWithBinsType, TB<:AbstractBinsType, F<:Function}

   n_particles = length(E);
   n_bins = length(B);
   # zero out offspring counts
   @. E.o = 0;
   @. B.target = 0;

   # identify nonempty bins
   non_empty_bins = findall(n->n>0, B.n);
   R = length(non_empty_bins);
   Ñ = zeros(n_bins);

   for p in non_empty_bins
      particle_ids = findall(isequal(p), E.b);
      Ñ[p] = sqrt(B.ν[p] * sum(E.ω[particle_ids] .* v².(E.ξ[particle_ids],t)));
   end

   if(sum(Ñ)>0)
      # normalize
      Ñ .= n_particles * Ñ/sum(Ñ);
      B.target .= (B.n .>0) .+ resample(n_particles-R, Ñ./n_particles);

      # compute number of offspring of each particle bin by bin
      for p in non_empty_bins
         # get particle indices for bin p
         particle_ids = findall(isequal(p), E.b);
         E.o[particle_ids] .= resample(B.target[p], E.ω[particle_ids]/B.ν[p]);
      end

   else
      # every particle copies itself
      B.target .= B.n;
      @. E.o = 1;
   end

   # resample the particles
   n_spawned = 0;
   for i in 1:n_particles
      # identify the bin of the current particle
      bin = E.b[i];
      for k in 1:E.o[i]
         E.ξ̂[k+n_spawned] = deepcopy(E.ξ[i]);
         E.ω̂[k+n_spawned] = B.ν[bin]/B.target[bin];
         E.b̂[k+n_spawned] = bin;
      end
      n_spawned += E.o[i];
   end
   E, B
end

"""
`uniform_selection!`: Uniformly select particles, ensuring each bin with
positive bin weight has at least one offspring.

### Arguments
* `E` - particle ensemble
* `B` - bin data structure
* `resample` - resampling scheme
"""
function uniform_selection!(E::TE, B::TB; resample=Systematic) where{TE<:EnsembleWithBinsType, TB<:AbstractBinsType}
   n_particles = length(E);
   n_bins = length(B);
   # zero out offspring counts
   @. E.o = 0;
   @. B.target = 0;

   # ensure each bin with walkers has at least one offspring
   non_empty_bins = findall(n->n>0, B.n);
   R = length(non_empty_bins);
   B.target[non_empty_bins] .= 1 .+ resample(n_particles-R, [1.0/R for j in 1:R]);

   # compute number of offspring of each particle bin by bin
   for p in non_empty_bins
      # get particle indices for bin p
      particle_ids = findall(isequal(p), E.b);
      E.o[particle_ids] .= resample(B.target[p], E.ω[particle_ids]/B.ν[p]);
   end

   # resample the particles
   n_spawned = 0;
   for i in 1:n_particles
      # identify the bin of the current particle
      bin = E.b[i];
      for k in 1:E.o[i]
         E.ξ̂[k+n_spawned] = deepcopy(E.ξ[i]);
         E.ω̂[k+n_spawned] = B.ν[bin]/B.target[bin];
         E.b̂[k+n_spawned] = bin;
      end
      n_spawned += E.o[i];
   end
   E, B
end
