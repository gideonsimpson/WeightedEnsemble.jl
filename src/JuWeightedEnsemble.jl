module JuWeightedEnsemble

using LinearAlgebra
using SparseArrays
using Printf
using NearestNeighbors
using Distributed
using SharedArrays
using StatsBase
using Distributions

# Load types
include("types.jl")
# Load data structures
include("ensemble.jl")
include("ensemble_without_bins.jl")
include("ensemble_with_bins.jl")
include("bins.jl")
# Load resampling algorithms
include("resampling.jl")
# Load selection algorithms
include("selection.jl")
# Load coarse model algorithms
include("coarse.jl")
# Load utility functions
include("utils.jl")
# Load WE methods
include("we.jl")

export EnsembleWithBins, Bins, EnsembleWithoutBins

end # end module
