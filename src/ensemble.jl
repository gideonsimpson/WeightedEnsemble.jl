"""
`Ensemble{TP, TW<:AbstractFloat, TB<:Integer}`: A particle ensemble structure designed for WE

### Fields

* `ξ̂` - particle positions after selection, before mutation
* `ξ` - particle positions after mutation
* `ω̂` - partice weights after selection, before mutation
* `ω` - partice weights after mutation
* `bin` - particle bin
* `offspring` - number of off spring of the particle
"""
struct Ensemble{TP, TW<:AbstractFloat, TB<:Integer}
   # positions of the walkers after resampling, before mutation
   ξ̂::Vector{TP}
   # positions of the walkers after mutation
   ξ::Vector{TP}
   # weights of the walkers after resampling, before mutation
   ω̂::Vector{TW}
   # weights of the walkers after mutation
   ω::Vector{TW}
   # category ("bin") type of each of the walkers
   bin::Vector{TB}
   # number of offspring of each particle
   offspring::Vector{TB}
end

function Base.eltype(E::Ensemble)
   return typeof((E.ξ̂[1], E.ξ[1], E.ω̂[1], E.ω[1], E.bin[1],E.offspring[1]))
end

function Base.push!(E::Ensemble, ξ̂, ξ, ω̂, ω, bin, offspring)
   push!(E.ξ̂, ξ̂);
   push!(E.ξ, ξ);
   push!(E.ω̂, ω̂);
   push!(E.ω, ω);
   push!(E.bin, bin);
   push!(E.offspring, offspring)
end

function Base.pop!(E::Ensemble)
   ξ̂ = pop!(E.ξ̂)
   ξ = pop!(E.ξ);
   ω̂ = pop!(E.ω̂);
   ω = pop!(E.ω);
   bin = pop!(E.bin)
   offspring = pop!(E.offspring)
   return ξ̂, ξ, ω̂, ω, bin, offspring
end

function Base.popfirst!(E::Ensemble)
   ξ̂ = popfirst!(E.ξ̂)
   ξ = popfirst!(E.ξ);
   ω̂ = popfirst!(E.ω̂);
   ω = popfirst!(E.ω);
   bin = popfirst!(E.bin)
   offspring = popfirst!(E.offspring)
   return ξ̂, ξ, ω̂, ω, bin, offspring
end

function Base.length(E::Ensemble)
   return length(E.ξ)
end

function Base.isempty(E::Ensemble)
   return isempty(E.ξ)
end

function Base.iterate(E::Ensemble, state = 1)

   if state > length(E)
      return nothing
   end
   return (E.ξ̂[state],E.ξ[state], E.ω̂[state],E.ω[state], E.bin[state], E.offspring[state]), state+1
end
