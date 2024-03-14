using Test
using WeightedEnsemble
using LinearAlgebra

@testset "Structures" begin
    @test include("structures/ensemble1.jl")
    @test include("structures/ensemble2.jl")
    @test include("structures/ensemble3.jl")
end

@testset "Utilities" begin
    @test include("utils/util1.jl")
end

@testset "Serial Doublewell" begin
    @test include("doublewell/doublewell_serial1.jl")
    @test include("doublewell/doublewell_serial2.jl")
end

@testset "Serial Muller" begin
    @test include("muller/muller_serial1.jl")
    @test include("muller/muller_serial2.jl")
end

