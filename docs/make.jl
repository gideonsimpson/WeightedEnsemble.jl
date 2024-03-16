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
using DocumenterCitations

bib = CitationBibliography(joinpath(@__DIR__, "src", "refs.bib"); style=:numeric)

makedocs(checkdocs=:none,
         sitename = "WeightedEnsemble.jl",
         modules  = [WeightedEnsemble],
         format=Documenter.HTML(
            # ...
            assets=String["assets/citations.css"],
            ),
         plugins=[bib],
         pages=[
                "Home" => "index.md",
                "Structures" => "struct1.md",
                "Mutation" =>"mutation/mutation1.md",
                "Selection" => ["selection/selection1.md", 
                "selection/allocation1.md", "selection/resampling1.md"],
                "Coarse Modeling" => "coarse1.md",
                "Utility Functions" => "util1.md",
                "Weighted Ensemble Methods" => "we1.md",
                "Examples"=>["examples/equil1.md"],
                "References"=>"refs1.md",
                "API"=>"api.md"
               ])
deploydocs(;
    repo="github.com/gideonsimpson/WeightedEnsemble.jl",
)