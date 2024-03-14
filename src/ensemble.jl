

"""
    Ensemble{TP, TF<:AbstractFloat, TI<:Integer, TD}

A particle ensemble structure designed for WE with bins.
### Fields

* `ξ̂` - particle positions after selection, before mutation
* `ξ` - particle positions after mutation
* `ω̂` - partice weights after selection, before mutation
* `ω` - partice weights after mutation
* `b̂` - particle bin after selection, before mutation
* `b` - particle bin after mutation
* `o` - number of offspring of the particle
* `d̂` - auxiliary data for each particle after selection, before mutation
* `d` - auxiliary data for each particle, after mutation
"""
struct Ensemble{TP, TF<:AbstractFloat, TI<:Integer, TD} <: AbstractEnsemble
   # positions of the walkers after selection, before mutation
   ξ̂::Vector{TP}
   # positions of the walkers after mutation
   ξ::Vector{TP}
   # weights of the walkers after selection, before mutation
   ω̂::Vector{TF}
   # weights of the walkers after mutation
   ω::Vector{TF}
   # category ("bin") type of each of the walkers after selection, before mutation
   b̂::Vector{TI}
   # category ("bin") type of each of the walkers after mutation
   b::Vector{TI}
   # number of offspring of each particle
   o::Vector{TI}
   # auxiliary data for each particle after selection, before mutation
   d̂::Vector{TD}
   # auxiliary data for each particle after mutation
   d::Vector{TD}
end

function Base.eltype(E::TE) where {TE<:Ensemble}
   return typeof((E.ξ̂[1], E.ξ[1], E.ω̂[1], E.ω[1],E.b̂[1], E.b[1],E.o[1], E.d̂[1], E.d[1]))
end

function Base.push!(E::TE, ξ̂, ξ, ω̂, ω, b̂, b, o, d̂, d) where {TE<:Ensemble}
   push!(E.ξ̂, ξ̂);
   push!(E.ξ, ξ);
   push!(E.ω̂, ω̂);
   push!(E.ω, ω);
   push!(E.b̂, b̂);
   push!(E.b, b);
   push!(E.o, o);
   push!(E.d̂, d̂);
   push!(E.d, d);
end

function Base.pop!(E::TE) where {TE<:Ensemble}
   ξ̂ = pop!(E.ξ̂);
   ξ = pop!(E.ξ);
   ω̂ = pop!(E.ω̂);
   ω = pop!(E.ω);
   b̂ = pop!(E.b̂);
   b = pop!(E.b);
   o = pop!(E.o);
   d̂ = pop!(E.d̂);
   d = pop!(E.d);
   return ξ̂, ξ, ω̂, ω, b̂, b, o, d̂, d
end

function Base.popfirst!(E::TE) where {TE<:Ensemble}
   ξ̂ = popfirst!(E.ξ̂);
   ξ = popfirst!(E.ξ);
   ω̂ = popfirst!(E.ω̂);
   ω = popfirst!(E.ω);
   b̂ = popfirst!(E.b̂);
   b = popfirst!(E.b);
   o = popfirst!(E.o);
   d̂ = popfirst!(E.d̂);
   d = popfirst!(E.d);
   return ξ̂, ξ, ω̂, ω, b̂, b, o, d̂, d
end

function Base.iterate(E::TE, state = 1) where {TE<:Ensemble}

   if state > length(E)
      return nothing
   end
   return (E.ξ̂[state],E.ξ[state], E.ω̂[state], E.ω[state], E.b̂[state], E.b[state], E.o[state], E.d̂[state], E.d[state]), state+1
end


# common ensmble functions

function Base.length(E::TE) where {TE<:AbstractEnsemble}
   return length(E.ξ)
end

function Base.isempty(E::TE) where {TE<:AbstractEnsemble}
   return isempty(E.ξ)
end

"""
    Ensemble(X, n_particles)

Construct an ensemble of `n_particles` with all particles in the state `X` with
uniform weights.
### Arguments
* `X` - state to be assigned to all particles
* `n_particles` - size of ensemble
"""
function Ensemble(X::TP, n_particles::TI) where {TP, TI<:Integer}
   ω = 1.0 / n_particles
   E = Ensemble{TP,Float64,TI, Nothing}([deepcopy(X) for _ in 1:n_particles],
      [deepcopy(X) for _ in 1:n_particles],
      ω * ones(n_particles), ω * ones(n_particles),
      zeros(TI, n_particles), zeros(TI, n_particles),
      zeros(TI, n_particles),
      [Nothing() for _ in 1:n_particles],
      [Nothing() for _ in 1:n_particles])
   return E
end

"""
    Ensemble(X)

Construct an ensemble with particles in the states specified in vector `X` with
uniform weights.
### Arguments
* `X` - array of states at which to initialize the particles
"""
function Ensemble(X::Vector{TP}) where {TP}
   n_particles = length(X);
   ω = 1.0 / n_particles
   E = Ensemble{TP,Float64,Integer,Nothing}(deepcopy(X),
      deepcopy(X),
      ω * ones(n_particles), ω * ones(n_particles),
      zeros(Integer, n_particles), zeros(Integer, n_particles),
      zeros(Integer, n_particles),
      [Nothing() for _ in 1:n_particles],
      [Nothing() for _ in 1:n_particles])
   return E
end

"""
    Ensemble(X, ω)

Construct an ensemble with particles in the states specified in vector `X` with
 weights specified in vector `ω`
### Arguments
* `X` - array of states at which to initialize the particles
* `ω` - array of weights at which to initialize the particles
"""
function Ensemble(X::Vector{TP}, ω::Vector{TF}) where {TP, TF<:AbstractFloat}
   n_particles = length(X)
   E = Ensemble{TP,Float64,Integer,Nothing}(deepcopy(X),
      deepcopy(X), deepcopy(ω), deepcopy(ω),
      zeros(Integer, n_particles), zeros(Integer, n_particles),
      zeros(Integer, n_particles),
      [Nothing() for _ in 1:n_particles],
      [Nothing() for _ in 1:n_particles])
   return E
end


