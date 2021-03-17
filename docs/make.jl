using Documenter
using GeometricObjects

# Setup for doctests in docstrings
DocMeta.setdocmeta!(GeometricObjects, :DocTestSetup, recursive = true,
    quote
        using GeometricObjects
    end
)

makedocs(;
    format = Documenter.HTML(prettyurls = get(ENV, "CI", nothing) == "true"),
    modules = [GeometricObjects],
    sitename = "GeometricObjects.jl",
    pages=[
        "Home" => "index.md",
    ],
    doctest = true, # :fix
)

deploydocs(
    repo = "github.com/KeitaNakamura/GeometricObjects.jl.git",
    devbranch = "main",
)
