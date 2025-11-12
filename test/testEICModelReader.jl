using EDM4hep.RootIO

@testset "EICModelReader" begin
    f = joinpath(@__DIR__, "../examples/ePIC/pythia8NCDIS_18x275_minQ2=10_beamEffects_xAngle=-0.025_hiDiv_1.0000.eicrecon.edm4eic.root")

    reader = RootIO.Reader(f)
    events = RootIO.get(reader, "events")
    
    @test reader.isRNTuple == false
    @test length(reader.collectionIDs) == 452

    @test length(events) == 356

    get_charged = RootIO.create_getter(reader, "ReconstructedChargedParticles"; selection=[:energy, :momentum, :charge])
    get_asso    = RootIO.create_getter(reader, "ReconstructedChargedParticleAssociations")
    get_mcps    = RootIO.create_getter(reader, "MCParticles"; selection=[:PDG, :momentum, :charge, :mass, :parents, :daughters])

    for evt in Iterators.take(events, 10)
        recps = get_charged(evt)
        asso  = get_asso(evt)
        mcps  = get_mcps(evt)
        @test eltype(recps) == ReconstructedParticle
        @test eltype(asso)  == MCRecoParticleAssociation
        @test eltype(mcps)  == MCParticle
    end
end
 