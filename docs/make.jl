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
        prettyurls = get(ENV, "CI", "false") == "true",
        canonical = "https://unfoldtoolbox.github.io/UnfoldMixedModels.jl",
        sidebar_sitename = false,
        edit_link = "main",
        assets = String[],
    ),
    pages = [
        "index.md",
        "Tutorials" => [
            "lmmERP (mass univariate)" => "tutorials/lmm_mu.md",
            "lmmERP (overlap correction)" => "tutorials/lmm_overlap.md",
        ],
        "HowTo" => ["P-values for mixedModels" => "howto/lmm_pvalues.md"],
        "Explanations" => [],
        "Reference" => [
            "API: Types" => "references/types.md",
            "API: Functions" => "references/functions.md",
        ],
        "Contributing" => ["90-contributing.md"],
        "Developer Guide" => ["91-developer.md"],
    ],
)

deploydocs(; repo = "github.com/unfoldtoolbox/UnfoldMixedModels.jl")
