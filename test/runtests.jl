using Test
using WeightedEnsemble

@testset "Serial Doublewell" begin
    @test include("doublewell/doublewell_serial1.jl")
    @test include("doublewell/doublewell_serial2.jl")
end

@testset "Serial Muller" begin
    @test include("muller/muller_serial1.jl")
    @test include("muller/muller_serial2.jl")
end