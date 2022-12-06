using Documenter, Accumulators

makedocs(
    modules = [Accumulators],
    sitename = "Accumulators.jl",
    authors = "Jeffrey A. Sarnoff <jeffrey.sarnoff@gmail.com>",
    pages = Any[
        "Overview" => "index.md",
        "Guide" => "guide.md",
        "Supported Operations" => "operations.md",
        "Refs" => "references.md"
    ]
)

deploydocs(
    repo = "github.com/JeffreySarnoff/Accumulators.jl.git",
    target = "build"
)
