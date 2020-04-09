

"""
`EnsembleWithloutBins{TP, TF, TI}`: A particle ensemble structure
designed for WE without bins.

### Fields

* `ξ̂` - particle positions after selection, before mutation
* `ξ` - particle positions after mutation
* `ω̂` - partice weights after selection, before mutation
* `ω` - partice weights after mutation
* `ρ` - particle resampling weight
* `o` - number of offspring of the particle
"""
struct EnsembleWithoutBins{TP, TF<:AbstractFloat, TI<:Integer} <: EnsembleWithoutBinsType
   # positions of the walkers after selection, before mutation
   ξ̂::Vector{TP}
   # positions of the walkers after mutation
   ξ::Vector{TP}
   # weights of the walkers after selection, before mutation
   ω̂::Vector{TF}
   # weights of the walkers after mutation
   ω::Vector{TF}
   # category ("bin") type of each of the walkers after selection, before mutation
   ρ::Vector{TF}
   # number of offspring of each particle
   o::Vector{TI}
end

function Base.eltype(E::TE) where {TE<:EnsembleWithoutBinsType}
   return typeof((E.ξ̂[1], E.ξ[1], E.ω̂[1], E.ω[1], E.ρ[1], E.o[1]))
end

function Base.push!(E::TE, ξ̂, ξ, ω̂, ω, ρ, o) where {TE<:EnsembleWithoutBinsType}
   push!(E.ξ̂, ξ̂);
   push!(E.ξ, ξ);
   push!(E.ω̂, ω̂);
   push!(E.ω, ω);
   push!(E.ρ, ρ);
   push!(E.o, o)
end

function Base.pop!(E::TE) where {TE<:EnsembleWithoutBinsType}
   ξ̂ = pop!(E.ξ̂);
   ξ = pop!(E.ξ);
   ω̂ = pop!(E.ω̂);
   ω = pop!(E.ω);
   ρ = pop!(E.ρ);
   o = pop!(E.o);
   return ξ̂, ξ, ω̂, ω, ρ, o
end

function Base.popfirst!(E::TE) where {TE<:EnsembleWithoutBinsType}
   ξ̂ = popfirst!(E.ξ̂);
   ξ = popfirst!(E.ξ);
   ω̂ = popfirst!(E.ω̂);
   ω = popfirst!(E.ω);
   ρ = popfirst!(E.ρ);
   o = popfirst!(E.o)
   return ξ̂, ξ, ω̂, ω, ρ, o
end

function Base.iterate(E::TE, state = 1) where {TE<:EnsembleWithoutBinsType}

   if state > length(E)
      return nothing
   end
   return (E.ξ̂[state],E.ξ[state], E.ω̂[state], E.ω[state], E.ρ[state], E.o[state]), state+1
end
