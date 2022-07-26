# WE routines

""" 
`trun_we`: Run a multithreaded WE simulation, optionally returning the
ensemble at each, step with

### Arguments
* `E₀` - initial particle ensemble
* `B₀` - initial bin data structure
* `sampler` - WE sampler functions data structure
* `n_we_steps` - number of steps in the WE run
### Optional Arguments
* `n_save_iters = 1` - save the ensemble and bins every `n_save_iters` iterations.  Set `n_save_iters=n_we_steps` to only save the final values.
"""
function trun_we(E₀::TE, B₀::TB, sampler::TWES, n_we_steps::Int; n_save_iters = 1) where
   {TE<:Ensemble, TB<:Bins, TWES<:WEsampler}

    E = deepcopy(E₀);
    B = deepcopy(B₀);
    E_trajectory = TE[];
    B_trajectory = TB[];

    n_particles = length(E);
   
    for t in 0:n_we_steps-1
        # first selection is at t = 0
        sampler.selection!(E, B, t);
        copy!(E.ω, E.ω̂);

        Threads.@threads for k in 1:n_particles
            E.ξ[k] = deepcopy(E.ξ̂[k]);
            E.d[k] = deepcopy(E.d̂[k]);            
            sampler.mutation!(E.ξ[k]);
        end
      
        # after mutation, time is t ↦ t+1
        sampler.rebin!(E, B, t+1);
        # analysis
        sampler.analysis!(E, B, t+1);

        if(mod(t+1, n_save_iters)==0)         
            push!(E_trajectory, deepcopy(E))
            push!(B_trajectory, deepcopy(B))
        end

    end
    
    return E_trajectory, B_trajectory

end


""" 
`trun_we_observables`: Run a multithreaded WE simulation, returning the values a
specified fucntions, `observables`, along the trajecotry.

### Arguments
* `E₀` - initial particle ensemble
* `B₀` - initial bin data structure
* `sampler` - WE sampler functions data structure
* `n_we_steps` - number of steps in the WE run
* `observables` - Tuple of scalar observable functions for the ergodic average
"""
function trun_we_observables(E₀::TE, B₀::TB, sampler::TWES, n_we_steps::Int, observables::Tuple{Vararg{<:Function,NO}}) where
   {TE<:Ensemble, TB<:Bins, TWES<:WEsampler, NO}
   
    E = deepcopy(E₀);
    B = deepcopy(B₀);
    observables_trajectory = zeros(NO, n_we_steps);

    n_particles = length(E);

    for t in 0:n_we_steps-1
        # first selection is at t = 0
        sampler.selection!(E, B, t);
        copy!(E.ω, E.ω̂);
        Threads.@threads for k in 1:n_particles
            E.ξ[k] = deepcopy(E.ξ̂[k]);
            E.d[k] = deepcopy(E.d̂[k]);
            sampler.mutation!(E.ξ[k]);
        end
        # after mutation, time is t ↦ t+1
        sampler.rebin!(E, B, t+1);
        # analysis
        sampler.analysis!(E, B, t+1);

        ntuple(k-> observables_trajectory[k,t+1] =(observables[k]).(E.ξ) ⋅ E.ω, NO)
    end

    return observables_trajectory

end

"""
`trun_we!`: Run an in place multithreaded WE simulation with

### Arguments
* `E` - particle ensemble
* `B` - bin data structure
* `sampler` - WE sampler functions data structure
* `n_we_steps` - number of steps in the WE run
"""
function trun_we!(E::TE, B::TB, sampler::TWES, n_we_steps::Int) where
   {TE<:Ensemble, TB<:Bins, TWES<:WEsampler}

    n_particles = length(E);

    for t in 0:n_we_steps-1
        # first selection is at t = 0
        sampler.selection!(E, B, t);
        copy!(E.ω, E.ω̂);
        Threads.@threads for k in 1:n_particles
            E.ξ[k] = deepcopy(E.ξ̂[k]);
            E.d[k] = deepcopy(E.d̂[k]);            
            sampler.mutation!(E.ξ[k]);
        end
        # after mutation, time is t ↦ t+1
        sampler.rebin!(E, B, t+1);
        # analysis
        sampler.analysis!(E, B, t+1);

   end
   E, B

end

