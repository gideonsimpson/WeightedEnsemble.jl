push!(LOAD_PATH,"../src/")
using WeightedEnsemble
using Documenter
makedocs(
         sitename = "WeightedEnsemble.jl",
         modules  = [WeightedEnsemble],
         pages=[
                "Home" => "index.md"
                "Structures" => "structures.md"
                "Resampling" => "resampling.md"
                "Bin Allocation" => "allocation.md"
                "Particle Selection" => "selection.md"
                "Coarse Modeling" => "coarse_models.md"
                "Utility Functions" => "utilities.md"
                "Weighted Ensemble Methods" => "weighted_ensemble_algs.md"
               ])
deploydocs(;
    repo="github.com/liamfdoherty/WeightedEnsemble.jl",
)