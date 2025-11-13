using Documenter
using EDM4hep
using EDM4hep.RootIO
using Literate
using Preferences

project = @__DIR__
function process_literate(names...)
    examples_mds = []
    for name in names
        run(`julia --threads=auto --project=$project docs/literate.jl $name`)
        push!(examples_mds, "examples/$name.md")
    end
    return examples_mds
end

savmodel = EDM4hep.get_edmodel()

# Set preferences for the examples "hep" model
EDM4hep.set_edmodel("hep")
FCC_mds    = process_literate("analysis_mH_recoil", "analysis_mH_recoil-MT")

# Set preferences for the examples "eic" model
EDM4hep.set_edmodel("eic")
EIC_mds    = process_literate("ePIC_simulation")

# Restore saved model
EDM4hep.set_edmodel(savmodel)

makedocs(;
    modules=[EDM4hep, EDM4hep.RootIO],
    format = Documenter.HTML(
        prettyurls = Base.get(ENV, "CI", nothing) == "true",
        repolink="https://github.com/JuliaHEP/EDM4hep.jl",
        size_threshold=1_000_000,
        size_threshold_warn=400_000,
    ),
    pages=[
        "Introduction" => "index.md",
        "Public APIs" => "api.md",
        "Release Notes" => "release_notes.md",
        "Examples" => [ 
            "FCC" => FCC_mds,
            "EIC" => EIC_mds
            ],
    ],
    checkdocs=:exports,
    repo="https://github.com/JuliaHEP/EDM4hep.jl/blob/{commit}{path}#L{line}",
    sitename="EDM4hep.jl",
    authors="Pere Mato",
)

deploydocs(;
    repo="github.com/JuliaHEP/EDM4hep.jl",
    devbranch="v0-patches",
    push_preview = true
)
