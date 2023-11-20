using Documenter
using CamiMath

makedocs(;
    modules=[CamiMath],
    authors="<walra356@planet.nl> and contributors",
    repo = "github.com/walra356/CamiMath.jl.git",
    sitename="CamiMath.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        repolink="https://walra356.github.io/CamiMath.jl",
        assets=String[]
    ),
    pages=[
        "Home" => "index.md",
        "Library" => "man/library.md",
        "Index" => "man/index.md",
    ]
)

deploydocs(;
    repo = "github.com/walra356/CamiMath.jl.git",
    devbranch = "main"
)
