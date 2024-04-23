using Test
using WeightedEnsemble
using LinearAlgebra
using Random
using BasicMD
using ForwardDiff
using Optim
using TestLandscapes: Muller

@testset "Ensembles" begin
    @test include("structures/ensemble1.jl")
    @test include("structures/ensemble2.jl")
    @test include("structures/ensemble3.jl")
end

@testset "Bins" begin
    @test include("structures/bins1.jl")
    @test include("structures/bins2.jl")
end

@testset "Utilities" begin
    @test include("utils/util1.jl")
end

# @testset "Selection" begin
#     @test include("selection/select1.jl")
#     @test include("selection/select2.jl")
#     @test include("selection/select3.jl")
# end

@testset "Mutation" begin
    @test include("mutation/mutation1.jl")
end

@testset "Serial Doublewell" begin
    @test include("doublewell/doublewell_serial1.jl")
    @test include("doublewell/doublewell_serial2.jl")
end

@testset "Serial Muller" begin
    @test include("muller/muller_serial1.jl")
    @test include("muller/muller_serial2.jl")
end

