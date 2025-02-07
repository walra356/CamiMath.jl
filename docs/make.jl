using Documenter
using DocumenterInterLinks
using CamiMath

makedocs(;
    modules=[CamiMath],
    authors="<walra356@planet.nl> and contributors",
    sitename="CamiMath.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        assets=String[]
    ),
    pages=[
        "Home" => "index.md",
        "Library" => "man/library.md",
        "Polynomials" => "man/polynomials.md",
        "Vector Coupling" => "man/vectorcoupling.md",
        "Toolbox" => "man/toolbox.md",
        "Index" => "man/index.md"
    ]    
)

deploydocs(;
    repo = "github.com/walra356/CamiMath.jl.git",
    devbranch = "main"
)
