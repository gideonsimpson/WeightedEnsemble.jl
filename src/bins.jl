
"""
`Bins{TS, TW, TBI, TT}`: A bin structure designed for WE

### Fields

* `Ω` - structure containing information for uniquely identifying each bin
* `n` - number of particles in each bin
* `target` - target number of particles in each bin
* `ν` - weight of each bin
* `d` - auxiliary bin data
"""
struct Bins{TS, TW<:AbstractFloat, TBI<:Integer, TT<:Real, TD} <: AbstractBins
   # structure which identifies the bins
   Ω::Vector{TS}
   # number of walkers in each bin
   n::Vector{TBI}
   # target number of walkers in each bin - specify if real/int with TT
   target::Vector{TT}
   # weight associated with each bin
   ν::Vector{TW}
   # auxiliary data structure for each bin
   d::Vector{TD}
end

function Base.push!(B::TB, Ω, n, target, ν, d) where {TB<:Bins}
   push!(B.Ω, Ω);
   push!(B.n, n);
   push!(B.target, target);
   push!(B.ν, ν);
   push!(B.d, d)
end

function Base.pop!(B::TB) where {TB<:Bins}
   Ω = pop!(B.Ω);
   n = pop!(B.n);
   target = pop!(B.target);
   ν = pop!(B.ν);
   d = pop!(B.d)
   #particles = pop!(B.particles)
   return Ω, n, target, ν, d
end

function Base.popfirst!(B::TB) where {TB<:Bins}
   Ω = popfirst!(B.Ω);
   n = popfirst!(B.n);
   target = popfirst!(B.target);
   ν = popfirst!(B.ν);
   d = popfirst!(B.d);
   #particles = pop!(B.particles)
   return Ω, n, target, ν, d
end

function Base.length(B::TB) where {TB<:Bins}
   return length(B.Ω)
end

function Base.eltype(B::TB) where {TB<:Bins}
   return typeof((B.Ω, B.n, B.target, B.ν, B.d))
end

function Base.iterate(B::TB, state = 1) where {TB<:Bins}

   if state > length(B)
      return nothing
   end
   return (B.Ω[state], B.n[state], B.target[state], B.ν[state], B.d[state]), state+1
end
