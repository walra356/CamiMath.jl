using Documenter
using CamMath

makedocs(;
    modules = [CamMath],
    authors = "<walra356@planet.nl> and contributors",
    #repo = "github.com/walra356/CamMath.jl.git",
    sitename = "CamMath.jl",
    format = Documenter.HTML(;
        prettyurls = get(ENV, "CI", "false") == "true",
        #canonical="https://walra356.github.io/CamMath.jl",
        assets=String[]
    ),
    pages=[
        "Home" => "index.md",
    ]
)

deploydocs(;
    repo = "github.com/walra356/CamMath.jl.git",
    devbranch = "main"
)
