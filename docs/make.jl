using UnfoldMixedModels
using Documenter

DocMeta.setdocmeta!(
    UnfoldMixedModels,
    :DocTestSetup,
    :(using UnfoldMixedModels);
    recursive = true,
)

const page_rename = Dict("developer.md" => "Developer docs") # Without the numbers
const numbered_pages = [
    file for file in readdir(joinpath(@__DIR__, "src")) if
    file != "index.md" && splitext(file)[2] == ".md"
]

makedocs(;
    modules = [UnfoldMixedModels],
    authors = "Benedikt Ehinger <benedikt.ehinger@vis.uni-stuttgart.de>",
    repo = "https://github.com/unfoldtoolbox/UnfoldMixedModels.jl/blob/{commit}{path}#{line}",
    sitename = "UnfoldMixedModels.jl",
    format = Documenter.HTML(;
        canonical = "https://unfoldtoolbox.github.io/UnfoldMixedModels.jl",
    ),
    pages = ["index.md"; numbered_pages],
)

deploydocs(; repo = "github.com/unfoldtoolbox/UnfoldMixedModels.jl")
