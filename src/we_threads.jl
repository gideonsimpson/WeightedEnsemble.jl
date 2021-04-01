# WE routines

""" 
`trun_we`: Run a multithreaded WE simulation, optionally returning the
ensemble at each, step with

### Arguments
* `E₀` - initial particle ensemble
* `B₀` - initial bin data structure
* `mutation!` - mutation function
* `selection!` - selection scheme
* `rebin!` - rebin and update particles and bins
* `n_we_steps` - number of steps in the WE run
* `save_trajectory=true` - save the ensemble and bins at each iteration.  if false, only returns the final state
"""
function trun_we(E₀::TE, B₀::TB, mutation!::FM, selection!::FS, rebin!::FR, n_we_steps::Int, save_trajectory=true) where
   {TE<:EnsembleWithBins, TB<:AbstractBins, FM<:Function, FS<:Function, FR<:Function}

    E = deepcopy(E₀);
    B = deepcopy(B₀);
    E_trajectory = TE[];
    B_trajectory = TB[];

    if(save_trajectory)
        push!(E_trajectory, deepcopy(E))
        push!(B_trajectory, deepcopy(B))
    end

    n_particles = length(E);
   
    for t in 0:n_we_steps-1
        # first selection is at t = 0
        selection!(E, B, t);
        copy!(E.ω, E.ω̂);

        Threads.@threads for k in 1:n_particles
            copy!(E.ξ[k],E.ξ̂[k]);
            mutation!(E.ξ[k]);
        end
      
        # after mutation, time is t ↦ t+1
        rebin!(E, B, t+1);
        if(save_trajectory || t==n_we_steps-1)
            push!(E_trajectory, deepcopy(E))
            push!(B_trajectory, deepcopy(B))
        end
    end
    
    return E_trajectory, B_trajectory

end

""" 
`trun_we`: Run a multithreaded WE simulation, optionally returning the ensemble
at each, step with

### Arguments
* `E₀` - initial particle ensemble
* `mutation!` - mutation function
* `selection!` - selection scheme
* `analysis!` - perform any post mutation updates
* `n_we_steps` - number of steps in the WE run
* `save_trajectory=true` - save the ensemble and bins at each iteration.  if false, only returns the final state
"""
function trun_we(E₀::TE, mutation!::FM, selection!::FS, analysis!::FA, n_we_steps::Int, save_trajectory=true) where
   {TE<:AbstractEnsemble, FM<:Function, FS<:Function, FA<:Function}

   E = deepcopy(E₀);
   E_trajectory = TE[];

   if(save_trajectory)
      push!(E_trajectory, deepcopy(E))
   end

    n_particles = length(E);


    for t in 0:n_we_steps-1
        # first selection is at t = 0
        selection!(E, B, t);
        copy!(E.ω, E.ω̂);
        Threads.@threads for k in 1:n_particles
            copy!(E.ξ[k],E.ξ̂[k]);
            mutation!(E.ξ[k]);
        end
        # after mutation, time is t ↦ t+1
        analysis!(E, t+1);

        if(save_trajectory || t==n_we_steps-1)
            push!(E_trajectory, deepcopy(E))
        end
    end

   return E_trajectory

end

""" 
`trun_we_observables`: Run a multithreaded WE simulation, returning the values a
specified fucntions, `observables`, along the trajecotry.

### Arguments
* `E₀` - initial particle ensemble
* `B₀` - initial bin data structure
* `mutation!` - mutation function
* `selection!` - selection scheme
* `rebin!` - rebin and update particles and bins
* `n_we_steps` - number of steps in the WE run
* `observables` - Tuple of scalar observable functions for the ergodic average
"""
function trun_we_observables(E₀::TE, B₀::TB, mutation!::FM, selection!::FS, rebin!::FR, n_we_steps::Int, observables::Tuple{Vararg{<:Function,NO}}) where
   {TE<:EnsembleWithBins, TB<:AbstractBins, FM<:Function, FS<:Function, FR<:Function, NO}
   
    quote
        E = deepcopy(E₀);
        B = deepcopy(B₀);
        observables_trajectory = zeros($NO, n_we_steps);

        n_particles = length(E);

        for t in 0:n_we_steps-1
            # first selection is at t = 0
            selection!(E, B, t);
            copy!(E.ω, E.ω̂);
            Threads.@threads for k in 1:n_particles
                copy!(E.ξ[k],E.ξ̂[k]);
                mutation!(E.ξ[k]);
            end
            # after mutation, time is t ↦ t+1
            rebin!(E, B, t+1);
            Base.Cartesian.@nexprs $NO k -> observables_trajectory[k,t+1] = (observables[k]).(E.ξ) ⋅ E.ω;
        end

        return observables_trajectory
    end

end

""" 
`trun_we_observables`: Run a multithreaded WE simulation, returning the
values a specified fucntion, `f`, along the trajecotry.


### Arguments
* `E₀` - initial particle ensemble
* `mutation` - mutation function
* `selection!` - selection scheme
* `analysis!` - perform any post mutation updates
* `n_we_steps` - number of steps in the WE run
* `observables` - Tuple of scalar observable functions for the ergodic average
"""
function trun_we_observables(E₀::TE, mutation!::FM, selection!::FS, analysis!::FA, n_we_steps::Int, observables::Tuple{Vararg{<:Function,NO}}) where
   {TE<:AbstractEnsemble, FM<:Function, FS<:Function, FA<:Function, NO}

    quote
        E = deepcopy(E₀);
        observables_trajectory = zeros($NO, n_we_steps);

        n_particles = length(E);

        for t in 0:n_we_steps-1
            # first selection is at t = 0
            selection!(E, B, t);
            copy!(E.ω, E.ω̂);
            Threads.@threads for k in 1:n_particles
                copy!(E.ξ[k],E.ξ̂[k]);
                mutation!(E.ξ[k]);
            end
            # after mutation, time is t ↦ t+1
            analysis!(E, t+1);
            Base.Cartesian.@nexprs $NO k -> observables_trajectory[k,t+1] = (observables[k]).(E.ξ) ⋅ E.ω;
        end

        return observables_trajectory
    end

end


"""
`trun_we!`: Run an in place multithreaded WE simulation with

### Arguments
* `E` - particle ensemble
* `B` - bin data structure
* `mutation` - mutation function
* `selection!` - selection scheme
* `rebin!` - rebin and update particles and bins
* `n_we_steps` - number of steps in the WE run
"""
function trun_we!(E::TE, B::TB, mutation!::FM, selection!::FS, rebin!::FR, n_we_steps::Int) where
   {TE<:EnsembleWithBins, TB<:AbstractBins, FM<:Function, FS<:Function, FR<:Function}

    n_particles = length(E);

    for t in 0:n_we_steps-1
        # first selection is at t = 0
        selection!(E, B, t);
        copy!(E.ω, E.ω̂);
        Threads.@threads for k in 1:n_particles
            copy!(E.ξ[k],E.ξ̂[k]);
            mutation!(E.ξ[k]);
        end
        # after mutation, time is t ↦ t+1
        rebin!(E, B, t+1);
   end
   E, B

end

"""
`trun_we!`: Run an in place multithreaded WE simulation with

### Arguments
* `E` - particle ensemble
* `mutation` - mutation function
* `selection!` - selection scheme
* `analysis!` - perform any post mutation updates
* `n_we_steps` - number of steps in the WE run
"""
function trun_we!(E::TE, mutation!::FM, selection!::FS, analysis!::FA, n_we_steps::Int) where
   {TE<:AbstractEnsemble, FM<:Function, FS<:Function, FA<:Function}

    n_particles = length(E);

    for t in 0:n_we_steps-1
        # first selection is at t = 0
        selection!(E, t);
        copy!(E.ω, E.ω̂);
        Threads.@threads for k in 1:n_particles
            copy!(E.ξ[k],E.ξ̂[k]);
            mutation!(E.ξ[k]);
        end
        # after mutation, time is t ↦ t+1
        analysis!(E, t+1);
    end 
    E
end

