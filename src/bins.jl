abstract type AbstractBins end

"""
`Bins{TS, TW<:AbstractFloat, TB<:Integer, TT<:Real}`: A bin structure designed for WE

### Fields

* `Ω` - structure containing information for uniquely identifying each bin
* `n` - number of particles in each bin
* `target` - target number of particles in each bin
* `ν` - weight of each bin
"""
struct Bins{TS, TW<:AbstractFloat, TB<:Integer, TT<:Real} <: AbstractBins
   # structure which identifies the bins
   Ω::Vector{TS}
   # number of walkers in each bin
   n::Vector{TB}
   # target number of walkers in each bin - specify if real/int with TT
   target::Vector{TT}
   # weight associated with each bin
   ν::Vector{TW}
end

function Base.push!(B::Bins, Ω, n, target, ν)
   push!(B.Ω, Ω);
   push!(B.n, n);
   push!(B.target, target);
   push!(B.ν, ν);
   #push!(B.particles, particles)
end

function Base.pop!(B::Bins)
   Ω = pop!(B.Ω);
   n = pop!(B.n);
   target = pop!(B.target);
   ν = pop!(B.ν);
   #particles = pop!(B.particles)
   return Ω, n, target, ν
end

function Base.popfirst!(B::Bins)
   Ω = popfirst!(B.Ω);
   n = popfirst!(B.n);
   target = popfirst!(B.target);
   ν = popfirst!(B.ν);
   #particles = pop!(B.particles)
   return Ω, n, target, ν
end

function Base.length(B::Bins)
   return length(B.Ω)
end

function Base.eltype(B::Bins)
   return typeof((B.Ω, B.n, B.target, B.ν))
end

function Base.iterate(B::Bins, state = 1)

   if state > length(B)
      return nothing
   end
   return (B.Ω[state], B.n[state], B.target[state], B.ν[state]), state+1
end
