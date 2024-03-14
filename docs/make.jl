# push!(LOAD_PATH,"../src/")
if abspath(PROGRAM_FILE) == @__FILE__
    # When running the `make.jl` file as a script, automatically activate the
    # `docs` environment and dev-install the main package into that environment
    import Pkg
    Pkg.activate(@__DIR__)
    Pkg.develop(path=joinpath(@__DIR__, ".."))
    Pkg.instantiate()
end
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
                "Weighted Ensemble Methods" => "we1.md",
                "Examples"=>["examples/equil1.md"]
               ])
deploydocs(;
    repo="github.com/gideonsimpson/WeightedEnsemble.jl",
)