module WeightedEnsemble

using LinearAlgebra
using SparseArrays
using Printf
using NearestNeighbors
using Distributed
using SharedArrays
using StatsBase
using Distributions
using Arpack

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
export setup_Voronoi_bins
# Load WE methods
include("we_serial.jl")
export run_we, run_we!, run_we_observable
include("we_parallel.jl")
export prun_we, prun_we!, prun_we_observable
include("we_threads.jl")
export trun_we, trun_we!, trun_we_observable

export EnsembleWithBins, Bins, EnsembleWithoutBins



end # end module
