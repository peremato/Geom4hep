using Documenter, Literate, Geom4hep

const examplesdir = joinpath(@__DIR__, "src/examples")

#---Convert the Literate examples (.lit) to Markdown pages (executed), notebooks and scripts---------
function process_literate(names...)
    mds = String[]
    for name in names
        run(`$(Base.julia_cmd()) --project=$(@__DIR__) $(joinpath(@__DIR__, "literate.jl")) $name`)
        push!(mds, "examples/$name.md")
    end
    return mds
end

example_mds = process_literate("Shapes", "Volumes", "Navigation", "GDML", "RoundCube", "Benchmarks")

#---Build the documentation-----------------------------------------------------------------------
makedocs(;
    modules = [Geom4hep],
    format = Documenter.HTML(
        prettyurls = get(ENV, "CI", nothing) == "true",
        size_threshold = 8_000_000,
        size_threshold_warn = 4_000_000,
    ),
    pages = [
        "Introduction" => "index.md",
        "User Guide" => "guide.md",
        "Examples" => ["Running the examples" => "examples.md", example_mds...],
        "Public API" => "api.md",
        "Development" => "development.md",
        "Release Notes" => "releases.md",
    ],
    checkdocs = :exports,
    warnonly = [:missing_docs],
    repo = Remotes.GitHub("JuliaHEP", "Geom4hep.jl"),
    sitename = "Geom4hep.jl",
    authors = "Pere Mato",
)

#---Deploy the documentation-----------------------------------------------------------------------
deploydocs(;
    repo = "github.com/JuliaHEP/Geom4hep.jl",
    push_preview = true,
)
