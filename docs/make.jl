using Documenter, GliderADCP

makedocs(
    sitename = "GliderADCP.jl",
    modules = [GliderADCP],
    pages = [
        "Home" => "index.md",
        "API" => "api.md",
    ],
    checkdocs = :exports,
    warnonly = [:missing_docs],
)
