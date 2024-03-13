push!(LOAD_PATH,"../src/")
using WeightedEnsemble
using Documenter
makedocs(checkdocs=:none,
         sitename = "WeightedEnsemble.jl",
         modules  = [WeightedEnsemble],
         pages=[
                "Home" => "index.md",
                "Structures" => "struct1.md",
                "Mutation" =>"mutation1.md",
                "Selection" => ["selection1.md", "allocation1.md", "resampling1.md"],
                "Coarse Modeling" => "coarse1.md",
                "Utility Functions" => "util1.md",
                "Weighted Ensemble Methods" => "we1.md"
               ])
deploydocs(;
    repo="github.com/gideonsimpson/WeightedEnsemble.jl",
)