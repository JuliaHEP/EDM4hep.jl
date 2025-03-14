using Test
using EDM4hep

@testset "EDM4hep tests" verbose = true begin
    include("testComponents.jl")      # unit test basic components
    include("testMCParticle.jl")      # one-to-many relation
    include("testSimTrackerHit.jl")   # one-to-one relation
    include("testParticleID.jl")      # vector members
    include("testCluster.jl")         # several one-to-many and Vector members
    include("testCovMatrix.jl")       # CovMatrix utility functions
    #---ROOT I/O----------------------
    include("testRootReader.jl")      # TTree and RNTuple reader
    include("testRootReaderLegacy.jl")# Testing podio verion < 0.17
    include("testEDM4hepFile.jl")     # EDM4hep file reader
    #---Analysis----------------------
    include("testAnalysis.jl")        # Testing analysis interface (NThreads)
end
