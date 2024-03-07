push!(LOAD_PATH,"../src/")
using WeightedEnsemble
using Documenter
makedocs(checkdocs=:none,
         sitename = "WeightedEnsemble.jl",
         modules  = [WeightedEnsemble],
         pages=[
                "Home" => "index.md"
                "Structures" => "struct1.md"
                "Resampling" => "resampling1.md"
                "Bin Allocation" => "allocation1.md"
                "Particle Selection" => "selection1.md"
                "Coarse Modeling" => "coarse1.md"
                "Utility Functions" => "util1.md"
                "Weighted Ensemble Methods" => "we1.md"
               ])
deploydocs(;
    repo="github.com/gideonsimpson/WeightedEnsemble.jl",
)