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
        "Documentation" => "pages/documentation.md",
        "VectorCoupling" => "pages/vectorcoupling.md",
        "Toolbox" => "pages/toolbox.md",
        "Index" => "pages/index.md"
    ],
    Depth=1,
    doctest=true,
    doctestargs=["--project"],
    doctestcoverage=true,       # coverage report           (default: false)                                                
    doctestcoverageargs=["--project"],
    doctestcoverageformat=:lcov, # coverage report format    (default: :lcov)
    doctestcoveragefile="lcov.info", # coverage report file      (default: "lcov.info")
    doctestcoverageexclude=["src/CamiMath.jl"], # coverage report exclude   (default: [])               
)

deploydocs(;
    repo = "github.com/walra356/CamiMath.jl.git",
    devbranch = "main"
)
