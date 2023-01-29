using Documenter
using CamMath

makedocs(
    sitename = "CamMath",
    format = Documenter.HTML(),
    modules = [CamMath]
)

# Documenter can also automatically deploy documentation to gh-pages.
# See "Hosting Documentation" and deploydocs() in the Documenter manual
# for more information.
#=deploydocs(
    repo = "<repository url>"
)=#
