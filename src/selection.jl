"""
`build_value_vectors`: Assemble the value vectors associated with argument `u` for `n`
steps with bin transition matrix `T`.  Used for optimal allocation.

### Arguments
* `n` - number of steps
* `T` - bin transition matrix
* `u` - quantity of interest vector on the bin space
* `tol` - tolerance for negative values which may appear due to roundoff
"""
function build_value_vectors(n, T, u; tol=1e-15)
   n_bins = length(u);
   vvals = [zeros(n_bins) for j in 1:n];

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
         @. vvals[l+1] .= sqrt(max(v1 - v2, 0));
      else
         throw(DomainError(l,"Nontrivial negative entry in value vector"));
      end

   end
   return vvals
end


"""
`optimal_allocation_selection!`: Select particles according to the bins, using
the value vectors for optimal allocation.

### Arguments
* `E` - particle ensemble
* `B` - bin data structure
* `value_vectors` - value vectors for resampling
* `j` - j-th seletion step
* `resample` - resampling scheme
"""
function optimal_allocation_selection!(E::Ensemble, B::Bins, value_vectors,j; resample=Systematic)

   n_particles = length(E);
   n_bins = length(B);

   if(B.ν ⋅ value_vectors[j] <= 0)
      throw(DomainError(j,"Bin Weights ⟂ Value Vector"));
   end

   # find target number of offspring for the bins based on
   # bin weights
   R = count(B.n .>0); # count number of bins which must have offspring
   B.target .= (B.n .>0) .+ resample(n_particles-R, (B.ν .* value_vectors[j])/(B.ν ⋅ value_vectors[j]));

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
`uniform_selection!`: Uniformly select particles, ensuring each bin with
positive bin weight has at least one offspring.

### Arguments
* `E` - particle ensemble
* `B` - bin data structure
* `resample` - resampling scheme
"""
function uniform_selection!(E::Ensemble, B::Bins)
   n_particles = length(E);
   n_bins = length(B);
   # zero out offspring counts
   @. E.offspring = 0;

   # ensure each bin with walkers has at lesat one offspring
   R = count(B.n .>0); 
   for i in 1:n_bins
      # get particle indices for bin i
      particle_ids = findall(isequal(i), E.bin);
      if !isempty(particle_ids)
         E.offspring[particle_ids] .= rand(Multinomial(1,length(particle_ids)));
      end
   end
   # allocate the other walkers n_partilces - R walkers
   E.offspring .+= rand(Multinomial(n_particles-R, n_particles));
   # resample the particles
   n_spawned = 0;
   for i in 1:n_particles
      for k in 1:E.offspring[i]
         E.ξ̂[k+n_spawned] = deepcopy(E.ξ[i]);
         E.ω̂[k+n_spawned] = E.ω[i]/E.offspring[i];
      end
      n_spawned += E.offspring[i];
   end
   E.ω̂ .= E.ω̂/sum(E.ω̂);
   E, B
end
