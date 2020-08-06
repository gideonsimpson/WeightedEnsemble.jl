# WE routines

""" `run_we`: Run a serial WE simulation, optionally returning the ensemble at
each, step with

### Arguments
* `E₀` - initial particle ensemble
* `B₀` - initial bin data structure
* `mutation` - mutation function
* `selection!` - selection scheme
* `rebin!` - rebin and update particles and bins
* `n_we_steps` - number of steps in the WE run
* `save_trajectory=true` - save the ensemble and bins at each iteration.  if false, only returns the final state
"""
function run_we(E₀::TE, B₀::TB, mutation::FM, selection!::FS, rebin!::FR, n_we_steps::Int, save_trajectory=true) where
   {TE<:EnsembleWithBins, TB<:AbstractBins, FM<:Function, FS<:Function, FR<:Function}

   E = deepcopy(E₀);
   B = deepcopy(B₀);
   E_trajectory = TE[];
   B_trajectory = TB[];

   if(save_trajectory)
      push!(E_trajectory, deepcopy(E))
      push!(B_trajectory, deepcopy(B))
   end

   for t in 0:n_we_steps-1
      # first selection is at t = 0
      selection!(E, B, t);
      @. E.ω = E.ω̂;
      @. E.ξ = mutation(E.ξ̂);
      # after mutation, time is t ↦ t+1
      rebin!(E, B, t+1);

      if(save_trajectory || t==n_we_steps-1)
         push!(E_trajectory, deepcopy(E))
         push!(B_trajectory, deepcopy(B))
      end
   end

   return E_trajectory, B_trajectory

end

""" `run_we`: Run a serial WE simulation, optionally returning the ensemble at
each, step with

### Arguments
* `E₀` - initial particle ensemble
* `mutation` - mutation function
* `selection!` - selection scheme
* `analysis!` - perform any post mutation updates
* `n_we_steps` - number of steps in the WE run
* `save_trajectory=true` - save the ensemble and bins at each iteration.  if false, only returns the final state
"""
function run_we(E₀::TE, mutation::FM, selection!::FS, analysis!::FA, n_we_steps::Int, save_trajectory=true) where
   {TE<:AbstractEnsemble, FM<:Function, FS<:Function, FA<:Function}

   E = deepcopy(E₀);
   E_trajectory = TE[];

   if(save_trajectory)
      push!(E_trajectory, deepcopy(E))
   end

   for t in 0:n_we_steps-1
      # first selection is at t = 0
      selection!(E, B, t);
      @. E.ω = E.ω̂;
      @. E.ξ = mutation(E.ξ̂);
      # after mutation, time is t ↦ t+1
      analysis!(E, t+1);

      if(save_trajectory || t==n_we_steps-1)
         push!(E_trajectory, deepcopy(E))
      end
   end

   return E_trajectory

end


"""
`run_we!`: Run an in place serial WE simulation with

### Arguments
* `E` - particle ensemble
* `B` - bin data structure
* `mutation` - mutation function
* `selection!` - selection scheme
* `rebin!` - rebin and update particles and bins
* `n_we_steps` - number of steps in the WE run
"""
function run_we!(E::TE, B::TB, mutation::FM, selection!::FS, rebin!::FR, n_we_steps::Int) where
   {TE<:EnsembleWithBins, TB<:AbstractBins, FM<:Function, FS<:Function, FR<:Function}

   for t in 0:n_we_steps-1
      # first selection is at t = 0
      selection!(E, B, t);
      @. E.ω = E.ω̂;
      @. E.ξ = mutation(E.ξ̂);
      # after mutation, time is t ↦ t+1
      rebin!(E, B, t+1);
   end
   E, B

end

"""
`run_we!`: Run an in place serial WE simulation with

### Arguments
* `E` - particle ensemble
* `mutation` - mutation function
* `selection!` - selection scheme
* `analysis!` - perform any post mutation updates
* `n_we_steps` - number of steps in the WE run
"""
function run_we!(E::TE, mutation::FM, selection!::FS, analysis!::FA, n_we_steps::Int) where
   {TE<:AbstractEnsemble, FM<:Function, FS<:Function, FA<:Function}

   for t in 0:n_we_steps-1
      # first selection is at t = 0
      selection!(E, t);
      @. E.ω = E.ω̂;
      @. E.ξ = mutation(E.ξ̂);
      # after mutation, time is t ↦ t+1
      analysis!(E, t+1);
   end
   E

end


"""
`prun_we`: Run a parallel WE simulation, optionally returning the ensemble at
each step. This performs the mutation steps in parallel, and assumes a worker
pool has already been created.

### Arguments
* `E₀` - initial particle ensemble
* `B₀` - initial bin data structure
* `mutation` - mutation function
* `selection!` - selection scheme
* `rebin!` - rebin and update particles and bins
* `n_we_steps` - number of steps in the WE run
* `save_trajectory=true` - save the ensemble and bins at each iteration.  if false, only returns the final state
"""
function prun_we(E₀::TE, B₀::TB, mutation::FM, selection!::FS, rebin!::FR, n_we_steps::Int, save_trajectory=true) where
   {TE<:EnsembleWithBins, TB<:AbstractBins, FM<:Function, FS<:Function, FR<:Function}

   E = deepcopy(E₀);
   B = deepcopy(B₀);
   E_trajectory = TE[];
   B_trajectory = TB[];

   if(save_trajectory)
      push!(E_vals, deepcopy(E))
      push!(B_vals, deepcopy(B))
   end

   for t in 0:n_we_steps-1
      # first selection is at t = 0
      selection!(E, B, t);
      @. E.ω = E.ω̂;
      E.ξ .= pmap(mutation, E.ξ̂);
      # after mutation, time is t ↦ t+1
      rebin!(E, B, t+1);

      if(save_trajectory || t==n_we_steps-1)
         push!(E_vals, deepcopy(E))
         push!(B_vals, deepcopy(B))
      end
   end

   return E_trajectory, B_trajectory

end

"""
`prun_we`: Run a parallel WE simulation, optionally returning the ensemble at
each step. This performs the mutation steps in parallel, and assumes a worker
pool has already been created.

### Arguments
* `E₀` - initial particle ensemble
* `mutation` - mutation function
* `selection!` - selection scheme
* `analysis!` - perform any post mutation updates
* `n_we_steps` - number of steps in the WE run
* `save_trajectory=true` - save the ensemble and bins at each iteration.  if false, only returns the final state
"""
function prun_we(E₀::TE, mutation::FM, selection!::FS, analysis!::FA, n_we_steps::Int, save_trajectory=true) where
   {TE<:AbstractEnsemble, FM<:Function, FS<:Function, FA<:Function}

   E = deepcopy(E₀);
   E_trajectory = TE[];

   if(save_trajectory)
      push!(E_vals, deepcopy(E))
   end

   for t in 0:n_we_steps-1
      # first selection is at t = 0
      selection!(E, B, t);
      @. E.ω = E.ω̂;
      E.ξ .= pmap(mutation, E.ξ̂);
      # after mutation, time is t ↦ t+1
      analysis!(E, t+1);

      if(save_trajectory || t==n_we_steps-1)
         push!(E_vals, deepcopy(E))
      end
   end

   return E_trajectory

end


"""
`prun_we!`: Run an in place parallel WE simulation.  This performs the mutation
steps in parallel, and assumes a worker pool has already been created.

### Arguments
* `E` - particle ensemble
* `B` - bin data structure
* `mutation` - mutation function
* `selection!` - selection scheme
* `rebin!` - rebin and update particles and bins
* `n_we_steps` - number of steps in the WE run
"""
function prun_we!(E::TE, B::TB, mutation::FM, selection!::FS, rebin!::FR, n_we_steps::Int) where
   {TE<:EnsembleWithBins, TB<:AbstractBins, FM<:Function, FS<:Function, FR<:Function}

   for t in 0:n_we_steps-1
      # first selection is at t = 0
      selection!(E, B, t);
      @. E.ω = E.ω̂
      E.ξ .= pmap(mutation, E.ξ̂);
      # after mutation, time is t ↦ t+1
      rebin!(E, B, t+1);
   end
   E, B

end


"""
`prun_we!`: Run a parallel in place WE simulation.  This performs the mutation
steps in parallel, and assumes a worker pool has already been created.

### Arguments
* `E` - particle ensemble
* `mutation` - mutation function
* `selection!` - selection scheme
* `analysis!` - perform any post mutation updates
* `n_we_steps` - number of steps in the WE run
"""
function prun_we!(E::TE, mutation::FM, selection!::FS, analysis!::FA, n_we_steps::Int) where
   {TE<:AbstractEnsemble, FM<:Function, FS<:Function, FA<:Function}

   for t in 0:n_we_steps-1
      # first selection is at t = 0
      selection!(E, t);
      @. E.ω = E.ω̂
      E.ξ .= pmap(mutation, E.ξ̂);
      # after mutation, time is t ↦ t+1
      analysis!(E, t+1);
   end
   E, B

end
