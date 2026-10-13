using Documenter
using EDM4hep
using EDM4hep.RootIO
using Literate
using Preferences

project = @__DIR__
function process_literate(subdir, names...)
    mds = []
    for name in names
        run(`julia --threads=auto --project=$project docs/literate.jl $subdir $name`)
        push!(mds, "$subdir/$name.md")
    end
    return mds
end

savmodel = EDM4hep.get_edmodel()

# Set preferences for the examples "hep" model
EDM4hep.set_edmodel("hep")
# FCC examples are not regenerated: their input sample uses EDM4hep schema < 1.0, which this version cannot read.
# The pages are the ones rendered with EDM4hep.jl 0.4 (v0-patches branch).
FCC_mds    = ["examples/analysis_mH_recoil.md", "examples/analysis_mH_recoil-MT.md"]
tutorials_mds = process_literate("tutorials", "tutorial_edm", "tutorial_io")

# Set preferences for the examples "eic" model
EDM4hep.set_edmodel("eic")
EIC_mds    = process_literate("examples", "ePIC_simulation")

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
        "Tutorials" => tutorials_mds,
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
    devbranch="main",
    push_preview = true
)
