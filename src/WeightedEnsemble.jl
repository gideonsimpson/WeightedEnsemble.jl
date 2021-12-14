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
using IterativeSolvers

# Load types
include("types.jl")
# Load data structures
include("ensemble.jl")
include("bins.jl")
# Load resampling algorithms
include("resampling.jl")
# Load allocation algorithms
include("allocation.jl")
export trivial_allocation!, minimal_bin_allocation!, 
    targeted_bin_allocation!, optimal_bin_allocation!, 
    uniform_bin_allocation!, within_bin_allocation!
# Load selection algorithms
include("selection.jl")
export trivial_selection!, uniform_selection!,
    optimal_selection!, targeted_selection!, repopulate!
# Load coarse model algorithms
include("coarse_serial.jl")
include("coarse_parallel.jl")
include("coarse_threads.jl")
# Load utility functions
include("utils.jl")
export setup_Voronoi_bins, Dirac_to_Ensemble
# Load sampler structures
include("sampler.jl")
# Load WE methods
include("we_serial.jl")
export run_we, run_we!, run_we_observables
include("we_parallel.jl")
export prun_we, prun_we!, prun_we_observables
include("we_threads.jl")
export trun_we, trun_we!, trun_we_observables

export Ensemble, Bins, WEsampler, DistributedWEsampler

end # end module
