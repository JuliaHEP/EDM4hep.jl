using Test
using EDM4hep
using EDM4hep.RootIO

const frames = 3      # How many frames or events will be written
const vectorsize = 5  # For vector members, each vector member will have this size
mutable struct Counter
    count::Int
end
++(c::Counter) = c.count += 1
reset!(c::Counter) = c.count = 41

@testset "ReadEDM4hep$(T)file" for T in (:TTree, :RNTuple)
    f = joinpath(@__DIR__, "EDM4hep_example_$T-1.root")
    reader = RootIO.Reader(f)
    events = RootIO.get(reader, "events")

    @test length(events) == frames
    count = Counter(-1)
    reset!(count)
    cov3f = CovMatrix3f([++(count) for i in 1:6])
    reset!(count)
    cov4f = CovMatrix4f([++(count) for i in 1:10])
    reset!(count)
    cov6f = CovMatrix6f([++(count) for i in 1:21])

    for evt in events
        #---EventHeader----------------------------------------------------------------------------
        evth = RootIO.get(reader, evt, "EventHeader")[1]
        reset!(count)
        @test evth.eventNumber == ++(count)
        @test evth.runNumber == ++(count)
        @test evth.timeStamp == ++(count)
        @test evth.weight == 1
        for i in 1:vectorsize
            @test evth.weights[i] == ++(count)
        end
    
        #---MCParticleCollection-------------------------------------------------------------------
        mcp = RootIO.get(reader, evt, "MCParticleCollection")
        reset!(count)
        for p in mcp
            @test p.PDG == ++(count)
            @test p.generatorStatus == ++(count)
            @test p.simulatorStatus == ++(count)
            @test p.charge == ++(count)
            @test p.time == ++(count)
            @test p.mass == ++(count)
            @test p.vertex.x == ++(count)
            @test p.vertex.y == ++(count)
            @test p.vertex.z == ++(count)
            @test p.endpoint.x == ++(count)
            @test p.endpoint.y == ++(count)
            @test p.endpoint.z == ++(count)
            @test p.momentum.x == ++(count)
            @test p.momentum.y == ++(count)
            @test p.momentum.z == ++(count)
            @test p.momentumAtEndpoint.x == ++(count)
            @test p.momentumAtEndpoint.y == ++(count)
            @test p.momentumAtEndpoint.z == ++(count)
            @test p.spin.x == ++(count)
            @test p.spin.y == ++(count)
            @test p.spin.z == ++(count)   
        end 
        @test mcp[1].daughters[1] == mcp[2]
        @test mcp[1].parents[1] == mcp[3]

        #---SimTrackerHitCollection----------------------------------------------------------------
        sth = RootIO.get(reader, evt, "SimTrackerHitCollection")
        reset!(count)
        for h in sth
            @test h.cellID == ++(count)
            @test h.eDep == ++(count)
            @test h.time == ++(count)
            @test h.pathLength == ++(count)
            @test h.quality == ++(count)
            @test h.position.x == ++(count)
            @test h.position.y == ++(count)
            @test h.position.z == ++(count)
            @test h.momentum.x == ++(count)
            @test h.momentum.y == ++(count)
            @test h.momentum.z == ++(count)
            @test h.particle == mcp[1]
        end

        #---CaloHitContributionCollection----------------------------------------------------------
        chc = RootIO.get(reader, evt, "CaloHitContributionCollection")
        reset!(count)
        for c in chc
            @test c.PDG == ++(count)
            @test c.energy  == ++(count)
            @test c.time == ++(count)
            @test c.stepPosition.x == ++(count)
            @test c.stepPosition.y == ++(count)
            @test c.stepPosition.z == ++(count)
            @test c.particle == mcp[1]
        end 

        #---SimCalorimeterHitCollection------------------------------------------------------------
        sch = RootIO.get(reader, evt, "SimCalorimeterHitCollection")
        reset!(count)
        for h in sch
            @test h.cellID == ++(count)
            @test h.energy == ++(count)
            @test h.position.x == ++(count)
            @test h.position.y == ++(count)
            @test h.position.z == ++(count)
            @test h.contributions[1] == chc[1]
        end

        #---RawCalorimeterHitCollection------------------------------------------------------------
        rch = RootIO.get(reader, evt, "RawCalorimeterHitCollection")
        reset!(count)
        for h in rch
            @test h.cellID == ++(count)
            @test h.amplitude == ++(count)
            @test h.timeStamp == ++(count)
        end

        #---CalorimeterHitCollection----------------------------------------------------------------
        ch = RootIO.get(reader, evt, "CalorimeterHitCollection")
        reset!(count)
        for h in ch
            @test h.cellID == ++(count)
            @test h.energy == ++(count)
            @test h.energyError == ++(count)
            @test h.time == ++(count)
            @test h.position.x == ++(count)
            @test h.position.y == ++(count)
            @test h.position.z == ++(count)
            @test h.type == ++(count)
        end

        #---ParticleIDCollection-------------------------------------------------------------------
        recp = RootIO.get(reader, evt, "ReconstructedParticleCollection")
        pidc = RootIO.get(reader, evt, "ParticleIDCollection")
        reset!(count)
        for pid in pidc
            @test pid.type == ++(count)
            @test pid.PDG == ++(count)
            @test pid.algorithmType == ++(count)
            @test pid.likelihood == ++(count)
            @test pid.parameters ==  [++(count) for i in 1:vectorsize]
            @test pid.particle == recp[1]
        end   

        #---ClusterCollection-----------------------------------------------------------------------
        cc = RootIO.get(reader, evt, "ClusterCollection")
        reset!(count)
        for c in cc
            @test c.type == ++(count)
            @test c.energy == ++(count)
            @test c.energyError == ++(count)
            @test c.position.x == ++(count)
            @test c.position.y == ++(count)
            @test c.position.z == ++(count)
            @test c.positionError == cov3f
            @test c.iTheta == ++(count)
            @test c.phi == ++(count)
            @test c.directionError.x == ++(count)
            @test c.directionError.y == ++(count)
            @test c.directionError.z == ++(count)
            for i in 1:vectorsize
                @test c.shapeParameters[i] == ++(count)
                @test c.subdetectorEnergies[i] == ++(count)
            end
            @test c.clusters[1] == cc[1]
            @test c.hits[1] == ch[1]
        end

        #---TrackerHit3DCollection------------------------------------------------------------------
        th3d = RootIO.get(reader, evt, "TrackerHit3DCollection")
        reset!(count)
        for h in th3d
            @test h.cellID == ++(count)
            @test h.type == ++(count)
            @test h.quality == ++(count)
            @test h.time == ++(count)
            @test h.eDep == ++(count)
            @test h.eDepError == ++(count)
            @test th3d[1].position == Vector3d(++(count), ++(count), ++(count))
            @test h.covMatrix == cov3f
        end
       
        #---TrackerHitPlaneCollection----------------------------------------------------------------
        thp = RootIO.get(reader, evt, "TrackerHitPlaneCollection")
        reset!(count)
        for h in thp
            @test h.cellID == ++(count)
            @test h.type == ++(count)
            @test h.quality == ++(count)
            @test h.time == ++(count)
            @test h.eDep == ++(count)
            @test h.eDepError == ++(count)
            @test h.u == Vector2f(++(count), ++(count))
            @test h.v == Vector2f(++(count), ++(count))
            @test h.du == ++(count)
            @test h.dv == ++(count)
            @test thp[1].position == Vector3d(++(count), ++(count), ++(count))
            @test h.covMatrix == cov3f
        end

        #---RawTimeSeriesCollection-----------------------------------------------------------------
        rts = RootIO.get(reader, evt, "RawTimeSeriesCollection")
        reset!(count)
        for ts in rts
            @test ts.cellID == ++(count)
            @test ts.quality == ++(count)
            @test ts.time == ++(count)
            @test ts.charge == ++(count)
            @test ts.interval == ++(count)
            @test ts.adcCounts == [++(count) for i in 1:vectorsize]
        end

        #---TrackCollection------------------------------------------------------------------------
        tc = RootIO.get(reader, evt, "TrackCollection")
        reset!(count)
        for t in tc
            @test t.type == ++(count)
            @test t.chi2 == ++(count)
            @test t.ndf == ++(count)
            for (i, s) in t.trackStates |> enumerate
                @test t.subdetectorHitNumbers[i] == ++(count)
                @test t.subdetectorHoleNumbers[i] == ++(count)
                @test s.location == ++(count)
                @test s.D0 == ++(count)
                @test s.phi == ++(count)
                @test s.omega == ++(count)
                @test s.Z0 == ++(count)
                @test s.tanLambda == ++(count)
                @test s.time == ++(count)
                @test s.referencePoint == Vector3f(++(count), ++(count), ++(count))
                @test s.covMatrix == cov6f
            end
            @test t.Nholes == ++(count)
            @test t.trackerHits[1] == th3d[1]
            @test t.tracks[1] == tc[1]
        end

        #---VertexCollection------------------------------------------------------------------------
        vc = RootIO.get(reader, evt, "VertexCollection")
        reset!(count)
        for v in vc
            @test v.type == ++(count)
            @test v.chi2 == ++(count)
            @test v.ndf == ++(count)
            @test v.position == Vector3f(++(count), ++(count), ++(count))
            @test v.covMatrix == cov3f
            @test v.algorithmType == ++(count)
            @test v.parameters == [++(count) for i in 1:vectorsize]
            @test v.particles[1] == recp[1]
        end

        #---ReconstructedParticleCollection--------------------------------------------------------
        reset!(count)
        for p in recp
            @test p.PDG == ++(count)
            @test p.energy == ++(count)
            @test p.momentum == Vector3f(++(count), ++(count), ++(count))
            @test p.referencePoint == Vector3f(++(count), ++(count), ++(count))
            @test p.charge == ++(count)
            @test p.mass == ++(count)
            @test p.goodnessOfPID == ++(count) 
            @test p.covMatrix == cov4f
            @test p.clusters[1] == cc[1]
            @test p.tracks[1] == tc[1]
            @test p.particles[1] == recp[1]
            @test p.decayVertex == vc[1]
        end

        #---Links-----------------------------------------------------------------------------------
        function checkLink(name, from, to)
            links = RootIO.get(reader, evt, "$(name)Collection")
            reset!(count)
            for l in links
                @test l.weight == ++(count)
                @test l.from == from
                @test l.to == to
            end
        end

        checkLink("RecoMCParticleLink", recp[1], mcp[1])
        checkLink("CaloHitSimCaloHitLink", ch[1], sch[1])
        checkLink("TrackerHitSimTrackerHitLink", th3d[1], sth[1])
        checkLink("CaloHitMCParticleLink", ch[1], mcp[1])
        checkLink("ClusterMCParticleLink", cc[1], mcp[1])
        checkLink("TrackMCParticleLink", tc[1], mcp[1])
        checkLink("VertexRecoParticleLink", vc[1], recp[1])

        #---TimeSeriesCollection-------------------------------------------------------------------
        tsc = RootIO.get(reader, evt, "TimeSeriesCollection")
        reset!(count)
        for ts in tsc
            @test ts.cellID == ++(count)
            @test ts.time == ++(count)
            @test ts.interval == ++(count)
            @test ts.amplitude == [++(count) for i in 1:vectorsize]
        end

        #---RecDqdxCollection-----------------------------------------------------------------------
        rdq = RootIO.get(reader, evt, "RecDqdxCollection")
        reset!(count)
        for dq in rdq
            @test dq.dQdx == Quantity(++(count), ++(count), ++(count))
            @test dq.track == tc[1]
        end

        #---GeneratorEventParametersCollection-----------------------------------------------------
        gep = RootIO.get(reader, evt, "GeneratorEventParametersCollection")
        reset!(count)
        for p in gep
            @test p.eventScale == ++(count)
            @test p.alphaQED == ++(count)
            @test p.alphaQCD == ++(count)
            @test p.signalProcessId == ++(count)
            @test p.sqrts == ++(count)
            for i in 1:vectorsize
                @test p.crossSections[i] == ++(count)
                @test p.crossSectionErrors[i] == ++(count)
            end
            @test p.signalVertex[1] == mcp[1]
        end
        
        #---GeneratorPdfInfoCollection-------------------------------------------------------------
        gpi = RootIO.get(reader, evt, "GeneratorPdfInfoCollection")
        reset!(count)
        for p in gpi
            @test p.partonId == [++(count), ++(count)]
            @test p.lhapdfId == [++(count), ++(count)]
            @test p.x == [++(count), ++(count)]
            @test p.xf == [++(count), ++(count)]
            @test p.scale == ++(count)
        end

    end
end
