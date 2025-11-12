using Test
using EDM4hep

@testset "EDM4hep tests" verbose = true begin
    include("testComponents.jl")      # unit test basic components
    include("testMCParticle.jl")      # one-to-many relation
    include("testSimTrackerHit.jl")   # one-to-one relation
    include("testParticleID.jl")      # vector members
    #---ROOT I/O----------------------
    if EDM4hep.get_edmodel() == "hep"
        include("testCluster.jl")         # several one-to-many and Vector members
        include("testRootReader.jl")      # TTree and RNTuple reader
        include("testRootReaderLegacy.jl")# Testing podio version < 0.17
    #---Analysis----------------------
        include("testAnalysis.jl")        # Testing analysis interface (NThreads)
    elseif EDM4hep.get_edmodel() == "eic"
        include("testEICModelReader.jl")  # EICModel specific reader tests
    end
end
