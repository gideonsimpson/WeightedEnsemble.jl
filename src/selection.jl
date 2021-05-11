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
function optimal_allocation_selection!(E::TE, B::TB, v²::F, t::Int; resample=Systematic) where{TE<:EnsembleWithBins, TB<:AbstractBins, F<:Function}

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
         @inbounds copy!(E.ξ̂[k+n_spawned], E.ξ[i]);
         @inbounds E.ω̂[k+n_spawned] = B.ν[bin]/B.target[bin];
         @inbounds E.b̂[k+n_spawned] = bin;
      end
      @inbounds n_spawned += E.o[i];
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
function uniform_selection!(E::TE, B::TB; resample=Systematic) where{TE<:EnsembleWithBins, TB<:AbstractBins}
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
         @inbounds copy!(E.ξ̂[k+n_spawned], E.ξ[i]);
         @inbounds E.ω̂[k+n_spawned] = B.ν[bin]/B.target[bin];
         @inbounds E.b̂[k+n_spawned] = bin;
      end
      @inbounds n_spawned += E.o[i];
   end
   E, B
end


"""
`trivial_selection!`: Each parent has exactly one offspring

### Arguments
* `E` - particle ensemble
"""
function trivial_selection!(E::TE) where{TE<:EnsembleWithoutBins}
   
   @. E.o = 1;
   @. E.ω̂ = E.ω;
   copy!.(E.ξ̂, E.ξ);
   E
end

function trivial_selection!(E::TE) where{TE<:EnsembleWithBins}
   @. E.o = 1;
   @. E.ω̂ = E.ω;
   @. E.b̂ = E.b;
   copy!.(E.ξ̂, E.ξ);

   E
end
