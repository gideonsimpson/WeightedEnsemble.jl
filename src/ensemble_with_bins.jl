

"""
`EnsembleWithBins{TP, TF<:AbstractFloat, TI<:Integer}`: A particle ensemble structure
designed for WE with bins.

### Fields

* `ξ̂` - particle positions after selection, before mutation
* `ξ` - particle positions after mutation
* `ω̂` - partice weights after selection, before mutation
* `ω` - partice weights after mutation
* `b̂` - particle bin after selection, before mutation
* `b` - particle bin after mutation
* `o` - number of offspring of the particle
"""
struct EnsembleWithBins{TP, TF<:AbstractFloat, TI<:Integer} <: AbstractEnsemble
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
end

function Base.eltype(E::TE) where {TE<:EnsembleWithBins}
   return typeof((E.ξ̂[1], E.ξ[1], E.ω̂[1], E.ω[1],E.b̂[1], E.b[1],E.o[1]))
end

function Base.push!(E::TE, ξ̂, ξ, ω̂, ω, b̂, b, o) where {TE<:EnsembleWithBins}
   push!(E.ξ̂, ξ̂);
   push!(E.ξ, ξ);
   push!(E.ω̂, ω̂);
   push!(E.ω, ω);
   push!(E.b̂, b̂);
   push!(E.b, b);
   push!(E.o, o)
end

function Base.pop!(E::TE) where {TE<:EnsembleWithBins}
   ξ̂ = pop!(E.ξ̂);
   ξ = pop!(E.ξ);
   ω̂ = pop!(E.ω̂);
   ω = pop!(E.ω);
   b̂ = pop!(E.b̂);
   b = pop!(E.b);
   o = pop!(E.o);
   return ξ̂, ξ, ω̂, ω, b̂, b, o
end

function Base.popfirst!(E::TE) where {TE<:EnsembleWithBins}
   ξ̂ = popfirst!(E.ξ̂);
   ξ = popfirst!(E.ξ);
   ω̂ = popfirst!(E.ω̂);
   ω = popfirst!(E.ω);
   b̂ = popfirst!(E.b̂)
   b = popfirst!(E.b);
   o = popfirst!(E.o)
   return ξ̂, ξ, ω̂, ω, b̂, b, o
end

function Base.iterate(E::TE, state = 1) where {TE<:EnsembleWithBins}

   if state > length(E)
      return nothing
   end
   return (E.ξ̂[state],E.ξ[state], E.ω̂[state], E.ω[state], E.b̂[state], E.b[state], E.o[state]), state+1
end