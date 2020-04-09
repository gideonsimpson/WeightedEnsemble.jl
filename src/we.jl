# WE routines

"""
`run_we!`: Run a serial WE simulation with

### Arguments
* `E` - particle ensemble
* `B` - bin data structure
* `mutation` - mutation function
* `selection!` - selection scheme
* `rebin!` - rebin and update particles and bins
* `n_we_steps` - number of steps in the WE run
"""
function run_we!(E::TE, B::TB, mutation::FM, selection!::FS, rebin!::FR, n_we_steps::Int) where
   {TE<:EnsembleWithBinsType, TB<:AbstractBinsType, FM<:Function, FS<:Function, FR<:Function}

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
`run_we!`: Run a serial WE simulation with

### Arguments
* `E` - particle ensemble
* `mutation` - mutation function
* `selection!` - selection scheme
* `analysis!` - perform any post mutation updates
* `n_we_steps` - number of steps in the WE run
"""
function run_we!(E::TE, mutation::FM, selection!::FS, analysis!::FA, n_we_steps::Int) where
   {TE<:AbstractEnsembleType, FM<:Function, FS<:Function, FA<:Function}

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
`prun_we!`: Run a parallel WE simulation.  This performs the mutation steps
in parallel, and assumes a worker pool has already been created.

### Arguments
* `E` - particle ensemble
* `B` - bin data structure
* `mutation` - mutation function
* `selection!` - selection scheme
* `rebin!` - rebin and update particles and bins
* `n_we_steps` - number of steps in the WE run
"""
function prun_we!(E::TE, B::TB, mutation::FM, selection!::FS, rebin!::FR, n_we_steps::Int) where
   {TE<:EnsembleWithBinsType, TB<:AbstractBinsType, FM<:Function, FS<:Function, FR<:Function}

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
`prun_we!`: Run a parallel WE simulation.  This performs the mutation steps
in parallel, and assumes a worker pool has already been created.

### Arguments
* `E` - particle ensemble
* `mutation` - mutation function
* `selection!` - selection scheme
* `rebin!` - rebin and update particles and bins
* `n_we_steps` - number of steps in the WE run
"""
function prun_we!(E::TE, mutation::FM, selection!::FS, analysis!::FA, n_we_steps::Int) where
   {TE<:AbstractEnsembleType, FM<:Function, FS<:Function, FA<:Function}

   for t in 0:n_we_steps-1
      # first selection is at t = 0
      selection!(E, B, t);
      @. E.ω = E.ω̂
      E.ξ .= pmap(mutation, E.ξ̂);
      # after mutation, time is t ↦ t+1
      analysis!(E, B, t+1);
   end
   E, B

end
