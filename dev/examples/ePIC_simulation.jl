using EDM4hep
using EDM4hep.RootIO
using EDM4hep.Histograms
using EDM4hep.Analysis
using Plots
using LinearAlgebra

fname = "pythia8NCDIS_18x275_minQ2=10_beamEffects_xAngle=-0.025_hiDiv_1.0000.eicrecon.edm4eic.root"

# eic_server = "root://dtn-eic.jlab.org/"
# fpath = "/volatile/eic/EPIC/RECO/25.08.0/epic_craterlake/DIS/NC/18x275/minQ2=10/"
# file = eic_server * fpath * fname

file  = joinpath(@__DIR__, "../../../examples/ePIC", fname)
reader = RootIO.Reader(file)
events = RootIO.get(reader, "events");

histogram(vcat(events.ReconstructedChargedParticles_energy...), bins=range(1,50,100))

positive = vcat(events.ReconstructedChargedParticles_charge...) .> 0
negative = vcat(events.ReconstructedChargedParticles_charge...) .< 0
histogram(vcat(events.ReconstructedChargedParticles_energy...)[positive], bins=range(0,50,100))

histogram(vcat(events.ReconstructedChargedParticles_energy...)[negative], label="negative", bins=range(0,50,100), xlabel="GeV")
histogram!(vcat(events.ReconstructedChargedParticles_energy...)[positive], label="positive", bins=range(0,50,100), xlabel="GeV")

get_charged = RootIO.create_getter(reader, "ReconstructedChargedParticles"; selection=[:energy, :momentum, :charge])
get_asso    = RootIO.create_getter(reader, "ReconstructedChargedParticleAssociations")
get_mcps    = RootIO.create_getter(reader, "MCParticles"; selection=[:PDG, :momentum, :charge, :mass, :parents, :daughters])

# Basic resolution histogram
hresolu = H1D("Energy resolution [%]", 100, -10., 10.)

# Loop over events and fill histogram
for evt in events
    recps  = get_charged(evt)     # Get the coll. of reconstructed charged particles
    assocs = get_asso(evt)        # Get the coll. of associations `rec <-> mcp`
    mcps   = get_mcps(evt)        # Get the coll. of MC particles. Used to get .energy

    for recp in recps             # Loop over reconstructed particles
        ind = findfirst(x -> x.rec == recp, assocs)  # Find the association to ReconstructedParticle
        isnothing(ind) && continue                   # If no association found, skip to next recp
        ΔE = recp.energy - assocs[ind].sim.energy    # Calculate energy resolution
        push!(hresolu, 100*ΔE/recp.energy)           # Fill histogram
    end
end

plot(hresolu)

mcps = get_mcps(events[1])                   # Get the collection of MC particles in the first event
for mcp in mcps[1:10]                        # Loop over the first N MC particles
    length(mcp.daughters) == 2 || continue   # only consider particles with 2 daughters
    println(mcp.name, "  Energy: ", mcp.energy, " Status: ", mcp.generatorStatus)
    for d in mcp.daughters
        println("   --->: ", d.name, "  Energy: ", d.energy)
    end
end

mcps[length.(mcps.daughters) .== 0].energy |> sum

hsumenergy = H1D("Summed energy final particles", 100, 200., 300.; unit=:GeV)
for evt in events
    mcps   = get_mcps(evt)        # Get the coll. of MC particles. Used to get .energy
    sumenergy = mcps[length.(mcps.daughters) .== 0].energy |> sum
    push!(hsumenergy, sumenergy)
end
plot(hsumenergy)
