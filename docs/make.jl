using Documenter
using AstroScalingRelations

# The `format` below makes it so that urls are set to "pretty" if you are pushing them to a hosting service, and basic if you are just using them locally to make browsing easier.

DocMeta.setdocmeta!(AstroScalingRelations, :DocTestSetup, :(using AstroScalingRelations; import Unitful; import UnitfulAstro); recursive=true)

makedocs(
    sitename="AstroScalingRelations.jl",
    modules = [AstroScalingRelations],
    format = Documenter.HTML(;prettyurls = get(ENV, "CI", nothing) == "true"),
    authors = "Chris Garling",
    pages = ["index.md", "public_methods.md"],
    doctest=true
)

deploydocs(;
    repo = "github.com/cgarling/AstroScalingRelations.jl.git",
    versions = ["stable" => "v^", "v#.#"],
    push_preview=true,
)
