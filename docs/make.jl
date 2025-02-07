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
        "Library" => "pages/library.md",
        "Polynomials" => "pages/polynomials.md",
        "Vector Coupling" => "pages/vectorcoupling.md",
        "Toolbox" => "pages/toolbox.md",
        "Index" => "pages/index.md"
    ]    
)

deploydocs(;
    repo = "github.com/walra356/CamiMath.jl.git",
    devbranch = "main"
)
