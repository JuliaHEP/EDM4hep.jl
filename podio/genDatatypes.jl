"""
Calibrated Detector Data
- Author: EDM4hep authors
# Fields
- `cellID::UInt64`:  cell id 
- `time::Float32`:  begin time [ns]
- `interval::Float32`:  interval of each sampling [ns]
- `amplitude::PVector{Float32}`:  calibrated detector data 
# Methods
- `setAmplitude(object::TimeSeries, v::AbstractVector{Float32})`: assign a set of values to the `amplitude` vector member
"""
struct TimeSeries <: POD
    index::ObjectID{TimeSeries}      # ObjectID of himself
    #---Data Members
    cellID::UInt64                   #  cell id 
    time::Float32                    #  begin time [ns]
    interval::Float32                #  interval of each sampling [ns]
    #---VectorMembers
    amplitude::PVector{TimeSeries,Float32,1}  #  calibrated detector data 
end

function TimeSeries(;cellID=0, time=0, interval=0, amplitude=PVector{TimeSeries,Float32,1}())
    TimeSeries(-1, cellID, time, interval, amplitude)
end

function setAmplitude(o::TimeSeries, v::AbstractVector{Float32})
    iszero(o.index) && (o = register(o))
    o = @set o.amplitude = v
    update(o)
end
"""
Calorimeter hit
- Author: EDM4hep authors
# Fields
- `cellID::UInt64`:  detector specific (geometrical) cell id 
- `energy::Float32`:  energy of the hit [GeV]
- `energyError::Float32`:  error of the hit energy [GeV]
- `time::Float32`:  time of the hit [ns]
- `position::Vector3f`:  position of the hit in world coordinates [mm]
- `type::Int32`:  type of hit 
"""
struct CalorimeterHit <: POD
    index::ObjectID{CalorimeterHit}  # ObjectID of himself
    #---Data Members
    cellID::UInt64                   #  detector specific (geometrical) cell id 
    energy::Float32                  #  energy of the hit [GeV]
    energyError::Float32             #  error of the hit energy [GeV]
    time::Float32                    #  time of the hit [ns]
    position::Vector3f               #  position of the hit in world coordinates [mm]
    type::Int32                      #  type of hit 
end

function CalorimeterHit(;cellID=0, energy=0, energyError=0, time=0, position=Vector3f(), type=0)
    CalorimeterHit(-1, cellID, energy, energyError, time, position, type)
end

"""
Calorimeter Hit Cluster
- Author: EDM4hep authors
# Fields
- `type::Int32`:  flagword that defines the type of cluster 
- `energy::Float32`:  energy of the cluster [GeV]
- `energyError::Float32`:  error on the energy [GeV]
- `position::Vector3f`:  position of the cluster [mm]
- `positionError::CovMatrix3f`:  covariance matrix of the position 
- `iTheta::Float32`:  Polar angle of the cluster's intrinsic direction (used e.g. for vertexing). Not to be confused with the cluster position seen from IP 
- `phi::Float32`:  Azimuthal angle of the cluster's intrinsic direction (used e.g. for vertexing). Not to be confused with the cluster position seen from IP 
- `directionError::Vector3f`:  covariance matrix of the direction [mm**2]
- `shapeParameters::PVector{Float32}`:  shape parameters. This should be accompanied by a descriptive list of names in the shapeParameterNames collection level metadata, as a vector of strings with the same ordering 
- `subdetectorEnergies::PVector{Float32}`:  energy observed in a particular subdetector 
# Relations
- `clusters::Cluster`:  clusters that have been combined to this cluster 
- `hits::CalorimeterHit`:  hits that have been combined to this cluster 
# Methods
- `setShapeParameters(object::Cluster, v::AbstractVector{Float32})`: assign a set of values to the `shapeParameters` vector member
- `setSubdetectorEnergies(object::Cluster, v::AbstractVector{Float32})`: assign a set of values to the `subdetectorEnergies` vector member
- `pushToClusters(obj::Cluster, robj::Cluster)`: push related object to the `clusters` relation
- `popFromClusters(obj::Cluster)`: pop last related object from `clusters` relation
- `pushToHits(obj::Cluster, robj::CalorimeterHit)`: push related object to the `hits` relation
- `popFromHits(obj::Cluster)`: pop last related object from `hits` relation
"""
struct Cluster <: POD
    index::ObjectID{Cluster}         # ObjectID of himself
    #---Data Members
    type::Int32                      #  flagword that defines the type of cluster 
    energy::Float32                  #  energy of the cluster [GeV]
    energyError::Float32             #  error on the energy [GeV]
    position::Vector3f               #  position of the cluster [mm]
    positionError::CovMatrix3f       #  covariance matrix of the position 
    iTheta::Float32                  #  Polar angle of the cluster's intrinsic direction (used e.g. for vertexing). Not to be confused with the cluster position seen from IP 
    phi::Float32                     #  Azimuthal angle of the cluster's intrinsic direction (used e.g. for vertexing). Not to be confused with the cluster position seen from IP 
    directionError::Vector3f         #  covariance matrix of the direction [mm**2]
    #---VectorMembers
    shapeParameters::PVector{Cluster,Float32,1}  #  shape parameters. This should be accompanied by a descriptive list of names in the shapeParameterNames collection level metadata, as a vector of strings with the same ordering 
    subdetectorEnergies::PVector{Cluster,Float32,2}  #  energy observed in a particular subdetector 
    #---OneToManyRelations
    clusters::Relation{Cluster,Cluster,1}  #  clusters that have been combined to this cluster 
    hits::Relation{Cluster,CalorimeterHit,2}  #  hits that have been combined to this cluster 
end

function Cluster(;type=0, energy=0, energyError=0, position=Vector3f(), positionError=CovMatrix3f(), iTheta=0, phi=0, directionError=Vector3f(), shapeParameters=PVector{Cluster,Float32,1}(), subdetectorEnergies=PVector{Cluster,Float32,2}(), clusters=Relation{Cluster,Cluster,1}(), hits=Relation{Cluster,CalorimeterHit,2}())
    Cluster(-1, type, energy, energyError, position, positionError, iTheta, phi, directionError, shapeParameters, subdetectorEnergies, clusters, hits)
end

function pushToClusters(c::Cluster, o::Cluster)
    iszero(c.index) && (c = register(c))
    c = @set c.clusters = push(c.clusters, o)
    update(c)
end
function popFromClusters(c::Cluster)
    iszero(c.index) && (c = register(c))
    c = @set c.clusters = pop(c.clusters)
    update(c)
end
function pushToHits(c::Cluster, o::CalorimeterHit)
    iszero(c.index) && (c = register(c))
    c = @set c.hits = push(c.hits, o)
    update(c)
end
function popFromHits(c::Cluster)
    iszero(c.index) && (c = register(c))
    c = @set c.hits = pop(c.hits)
    update(c)
end
function setShapeParameters(o::Cluster, v::AbstractVector{Float32})
    iszero(o.index) && (o = register(o))
    o = @set o.shapeParameters = v
    update(o)
end
function setSubdetectorEnergies(o::Cluster, v::AbstractVector{Float32})
    iszero(o.index) && (o = register(o))
    o = @set o.subdetectorEnergies = v
    update(o)
end
"""
The Monte Carlo particle - based on the lcio::MCParticle.
- Author: EDM4hep authors
# Fields
- `PDG::Int32`:  PDG code of the particle 
- `generatorStatus::Int32`:  status of the particle as defined by the generator 
- `simulatorStatus::Int32`:  status of the particle from the simulation program - use BIT constants below 
- `charge::Float32`:  particle charge 
- `time::Float32`:  creation time of the particle in wrt. the event, e.g. for preassigned decays or decays in flight from the simulator [ns]
- `mass::Float64`:  mass of the particle [GeV]
- `vertex::Vector3d`:  production vertex of the particle [mm]
- `endpoint::Vector3d`:  endpoint of the particle [mm]
- `momentum::Vector3d`:  particle 3-momentum at the production vertex [GeV]
- `momentumAtEndpoint::Vector3d`:  particle 3-momentum at the endpoint [GeV]
- `spin::Vector3f`:  spin (helicity) vector of the particle 
# Relations
- `parents::MCParticle`:  The parents of this particle 
- `daughters::MCParticle`:  The daughters this particle 
# Methods
- `pushToParents(obj::MCParticle, robj::MCParticle)`: push related object to the `parents` relation
- `popFromParents(obj::MCParticle)`: pop last related object from `parents` relation
- `pushToDaughters(obj::MCParticle, robj::MCParticle)`: push related object to the `daughters` relation
- `popFromDaughters(obj::MCParticle)`: pop last related object from `daughters` relation
"""
struct MCParticle <: POD
    index::ObjectID{MCParticle}      # ObjectID of himself
    #---Data Members
    PDG::Int32                       #  PDG code of the particle 
    generatorStatus::Int32           #  status of the particle as defined by the generator 
    simulatorStatus::Int32           #  status of the particle from the simulation program - use BIT constants below 
    charge::Float32                  #  particle charge 
    time::Float32                    #  creation time of the particle in wrt. the event, e.g. for preassigned decays or decays in flight from the simulator [ns]
    mass::Float64                    #  mass of the particle [GeV]
    vertex::Vector3d                 #  production vertex of the particle [mm]
    endpoint::Vector3d               #  endpoint of the particle [mm]
    momentum::Vector3d               #  particle 3-momentum at the production vertex [GeV]
    momentumAtEndpoint::Vector3d     #  particle 3-momentum at the endpoint [GeV]
    spin::Vector3f                   #  spin (helicity) vector of the particle 
    #---OneToManyRelations
    parents::Relation{MCParticle,MCParticle,1}  #  The parents of this particle 
    daughters::Relation{MCParticle,MCParticle,2}  #  The daughters this particle 
end

function MCParticle(;PDG=0, generatorStatus=0, simulatorStatus=0, charge=0, time=0, mass=0, vertex=Vector3d(), endpoint=Vector3d(), momentum=Vector3d(), momentumAtEndpoint=Vector3d(), spin=Vector3f(), parents=Relation{MCParticle,MCParticle,1}(), daughters=Relation{MCParticle,MCParticle,2}())
    MCParticle(-1, PDG, generatorStatus, simulatorStatus, charge, time, mass, vertex, endpoint, momentum, momentumAtEndpoint, spin, parents, daughters)
end

function pushToParents(c::MCParticle, o::MCParticle)
    iszero(c.index) && (c = register(c))
    c = @set c.parents = push(c.parents, o)
    update(c)
end
function popFromParents(c::MCParticle)
    iszero(c.index) && (c = register(c))
    c = @set c.parents = pop(c.parents)
    update(c)
end
function pushToDaughters(c::MCParticle, o::MCParticle)
    iszero(c.index) && (c = register(c))
    c = @set c.daughters = push(c.daughters, o)
    update(c)
end
function popFromDaughters(c::MCParticle)
    iszero(c.index) && (c = register(c))
    c = @set c.daughters = pop(c.daughters)
    update(c)
end
"""
Generator pdf information
- Author: EDM4hep authors
# Fields
- `partonId::SVector{2,Int32}`:  Parton PDG id 
- `lhapdfId::SVector{2,Int32}`:  LHAPDF PDF id (see https://lhapdf.hepforge.org/pdfsets.html) 
- `x::SVector{2,Float64}`:  Parton momentum fraction 
- `xf::SVector{2,Float64}`:  PDF value 
- `scale::Float64`:  Factorisation scale [GeV]
"""
struct GeneratorPdfInfo <: POD
    index::ObjectID{GeneratorPdfInfo}  # ObjectID of himself
    #---Data Members
    partonId::SVector{2,Int32}       #  Parton PDG id 
    lhapdfId::SVector{2,Int32}       #  LHAPDF PDF id (see https://lhapdf.hepforge.org/pdfsets.html) 
    x::SVector{2,Float64}            #  Parton momentum fraction 
    xf::SVector{2,Float64}           #  PDF value 
    scale::Float64                   #  Factorisation scale [GeV]
end

function GeneratorPdfInfo(;partonId=zero(SVector{2,Int32}), lhapdfId=zero(SVector{2,Int32}), x=zero(SVector{2,Float64}), xf=zero(SVector{2,Float64}), scale=0)
    GeneratorPdfInfo(-1, partonId, lhapdfId, x, xf, scale)
end

"""
Monte Carlo contribution to SimCalorimeterHit
- Author: EDM4hep authors
# Fields
- `PDG::Int32`:  PDG code of the shower particle that caused this contribution 
- `energy::Float32`:  energy of the this contribution [G]
- `time::Float32`:  time of this contribution [ns]
- `stepPosition::Vector3f`:  position of this energy deposition (step) [mm]
# Relations
- `particle::MCParticle`:  primary MCParticle that caused the shower responsible for this contribution to the hit 
"""
struct CaloHitContribution <: POD
    index::ObjectID{CaloHitContribution}  # ObjectID of himself
    #---Data Members
    PDG::Int32                       #  PDG code of the shower particle that caused this contribution 
    energy::Float32                  #  energy of the this contribution [G]
    time::Float32                    #  time of this contribution [ns]
    stepPosition::Vector3f           #  position of this energy deposition (step) [mm]
    #---OneToOneRelations
    particle_idx::ObjectID{MCParticle}  #  primary MCParticle that caused the shower responsible for this contribution to the hit 
end

function CaloHitContribution(;PDG=0, energy=0, time=0, stepPosition=Vector3f(), particle=-1)
    CaloHitContribution(-1, PDG, energy, time, stepPosition, particle)
end

function Base.getproperty(obj::CaloHitContribution, sym::Symbol)
    if sym == :particle
        idx = getfield(obj, :particle_idx)
        return iszero(idx) ? nothing : convert(MCParticle, idx)
    else # fallback to getfield
        return getfield(obj, sym)
    end
end
"""
Simulated calorimeter hit
- Author: EDM4hep authors
# Fields
- `cellID::UInt64`:  ID of the sensor that created this hit 
- `energy::Float32`:  energy of the hit [GeV]
- `position::Vector3f`:  position of the hit in world coordinates [mm]
# Relations
- `contributions::CaloHitContribution`:  Monte Carlo step contributions 
# Methods
- `pushToContributions(obj::SimCalorimeterHit, robj::CaloHitContribution)`: push related object to the `contributions` relation
- `popFromContributions(obj::SimCalorimeterHit)`: pop last related object from `contributions` relation
"""
struct SimCalorimeterHit <: POD
    index::ObjectID{SimCalorimeterHit}  # ObjectID of himself
    #---Data Members
    cellID::UInt64                   #  ID of the sensor that created this hit 
    energy::Float32                  #  energy of the hit [GeV]
    position::Vector3f               #  position of the hit in world coordinates [mm]
    #---OneToManyRelations
    contributions::Relation{SimCalorimeterHit,CaloHitContribution,1}  #  Monte Carlo step contributions 
end

function SimCalorimeterHit(;cellID=0, energy=0, position=Vector3f(), contributions=Relation{SimCalorimeterHit,CaloHitContribution,1}())
    SimCalorimeterHit(-1, cellID, energy, position, contributions)
end

function pushToContributions(c::SimCalorimeterHit, o::CaloHitContribution)
    iszero(c.index) && (c = register(c))
    c = @set c.contributions = push(c.contributions, o)
    update(c)
end
function popFromContributions(c::SimCalorimeterHit)
    iszero(c.index) && (c = register(c))
    c = @set c.contributions = pop(c.contributions)
    update(c)
end
"""
Raw data of a detector readout
- Author: EDM4hep authors
# Fields
- `cellID::UInt64`:  detector specific cell id 
- `quality::Int32`:  quality flag for the hit 
- `time::Float32`:  time of the hit [ns]
- `charge::Float32`:  integrated charge of the hit [fC]
- `interval::Float32`:  interval of each sampling [ns]
- `adcCounts::PVector{Int32}`:  raw data (32-bit) word at i 
# Methods
- `setAdcCounts(object::RawTimeSeries, v::AbstractVector{Int32})`: assign a set of values to the `adcCounts` vector member
"""
struct RawTimeSeries <: POD
    index::ObjectID{RawTimeSeries}   # ObjectID of himself
    #---Data Members
    cellID::UInt64                   #  detector specific cell id 
    quality::Int32                   #  quality flag for the hit 
    time::Float32                    #  time of the hit [ns]
    charge::Float32                  #  integrated charge of the hit [fC]
    interval::Float32                #  interval of each sampling [ns]
    #---VectorMembers
    adcCounts::PVector{RawTimeSeries,Int32,1}  #  raw data (32-bit) word at i 
end

function RawTimeSeries(;cellID=0, quality=0, time=0, charge=0, interval=0, adcCounts=PVector{RawTimeSeries,Int32,1}())
    RawTimeSeries(-1, cellID, quality, time, charge, interval, adcCounts)
end

function setAdcCounts(o::RawTimeSeries, v::AbstractVector{Int32})
    iszero(o.index) && (o = register(o))
    o = @set o.adcCounts = v
    update(o)
end
"""
Generator event parameters
- Author: EDM4hep authors
# Fields
- `eventScale::Float64`:  event scale 
- `alphaQED::Float64`:  alpha_QED 
- `alphaQCD::Float64`:  alpha_QCD 
- `signalProcessId::Int32`:  id of signal process 
- `sqrts::Float64`:  sqrt(s) [GeV]
- `crossSections::PVector{Float64}`:  list of cross sections [pb]
- `crossSectionErrors::PVector{Float64}`:  list of cross section errors [pb]
# Relations
- `signalVertex::MCParticle`:  List of initial state MCParticle that are the source of the hard interaction 
# Methods
- `setCrossSections(object::GeneratorEventParameters, v::AbstractVector{Float64})`: assign a set of values to the `crossSections` vector member
- `setCrossSectionErrors(object::GeneratorEventParameters, v::AbstractVector{Float64})`: assign a set of values to the `crossSectionErrors` vector member
- `pushToSignalVertex(obj::GeneratorEventParameters, robj::MCParticle)`: push related object to the `signalVertex` relation
- `popFromSignalVertex(obj::GeneratorEventParameters)`: pop last related object from `signalVertex` relation
"""
struct GeneratorEventParameters <: POD
    index::ObjectID{GeneratorEventParameters}  # ObjectID of himself
    #---Data Members
    eventScale::Float64              #  event scale 
    alphaQED::Float64                #  alpha_QED 
    alphaQCD::Float64                #  alpha_QCD 
    signalProcessId::Int32           #  id of signal process 
    sqrts::Float64                   #  sqrt(s) [GeV]
    #---VectorMembers
    crossSections::PVector{GeneratorEventParameters,Float64,1}  #  list of cross sections [pb]
    crossSectionErrors::PVector{GeneratorEventParameters,Float64,2}  #  list of cross section errors [pb]
    #---OneToManyRelations
    signalVertex::Relation{GeneratorEventParameters,MCParticle,1}  #  List of initial state MCParticle that are the source of the hard interaction 
end

function GeneratorEventParameters(;eventScale=0, alphaQED=0, alphaQCD=0, signalProcessId=0, sqrts=0, crossSections=PVector{GeneratorEventParameters,Float64,1}(), crossSectionErrors=PVector{GeneratorEventParameters,Float64,2}(), signalVertex=Relation{GeneratorEventParameters,MCParticle,1}())
    GeneratorEventParameters(-1, eventScale, alphaQED, alphaQCD, signalProcessId, sqrts, crossSections, crossSectionErrors, signalVertex)
end

function pushToSignalVertex(c::GeneratorEventParameters, o::MCParticle)
    iszero(c.index) && (c = register(c))
    c = @set c.signalVertex = push(c.signalVertex, o)
    update(c)
end
function popFromSignalVertex(c::GeneratorEventParameters)
    iszero(c.index) && (c = register(c))
    c = @set c.signalVertex = pop(c.signalVertex)
    update(c)
end
function setCrossSections(o::GeneratorEventParameters, v::AbstractVector{Float64})
    iszero(o.index) && (o = register(o))
    o = @set o.crossSections = v
    update(o)
end
function setCrossSectionErrors(o::GeneratorEventParameters, v::AbstractVector{Float64})
    iszero(o.index) && (o = register(o))
    o = @set o.crossSectionErrors = v
    update(o)
end
"""
Event Header. Additional parameters are assumed to go into the metadata tree.
- Author: EDM4hep authors
# Fields
- `eventNumber::UInt64`:  event number 
- `runNumber::UInt32`:  run number 
- `timeStamp::UInt64`:  time stamp 
- `weight::Float64`:  event weight 
- `weights::PVector{Float64}`:  event weights in case there are multiple. **NOTE that weights[0] might not be the same as weight!** Event weight names should be stored using the edm4hep::EventWeights name in the file level metadata 
# Methods
- `setWeights(object::EventHeader, v::AbstractVector{Float64})`: assign a set of values to the `weights` vector member
"""
struct EventHeader <: POD
    index::ObjectID{EventHeader}     # ObjectID of himself
    #---Data Members
    eventNumber::UInt64              #  event number 
    runNumber::UInt32                #  run number 
    timeStamp::UInt64                #  time stamp 
    weight::Float64                  #  event weight 
    #---VectorMembers
    weights::PVector{EventHeader,Float64,1}  #  event weights in case there are multiple. **NOTE that weights[0] might not be the same as weight!** Event weight names should be stored using the edm4hep::EventWeights name in the file level metadata 
end

function EventHeader(;eventNumber=0, runNumber=0, timeStamp=0, weight=0, weights=PVector{EventHeader,Float64,1}())
    EventHeader(-1, eventNumber, runNumber, timeStamp, weight, weights)
end

function setWeights(o::EventHeader, v::AbstractVector{Float64})
    iszero(o.index) && (o = register(o))
    o = @set o.weights = v
    update(o)
end
"""
Raw calorimeter hit
- Author: EDM4hep authors
# Fields
- `cellID::UInt64`:  detector specific (geometrical) cell id 
- `amplitude::Int32`:  amplitude of the hit in ADC counts 
- `timeStamp::Int32`:  time stamp for the hit 
"""
struct RawCalorimeterHit <: POD
    index::ObjectID{RawCalorimeterHit}  # ObjectID of himself
    #---Data Members
    cellID::UInt64                   #  detector specific (geometrical) cell id 
    amplitude::Int32                 #  amplitude of the hit in ADC counts 
    timeStamp::Int32                 #  time stamp for the hit 
end

function RawCalorimeterHit(;cellID=0, amplitude=0, timeStamp=0)
    RawCalorimeterHit(-1, cellID, amplitude, timeStamp)
end

"""
Tracker hit
- Author: EDM4hep authors
# Fields
- `cellID::UInt64`:  ID of the sensor that created this hit 
- `type::Int32`:  type of raw data hit 
- `quality::Int32`:  quality bit flag of the hit 
- `time::Float32`:  time of the hit [ns]
- `eDep::Float32`:  energy deposited on the hit [GeV]
- `eDepError::Float32`:  error measured on eDep [GeV]
- `position::Vector3d`:  hit position [mm]
- `covMatrix::CovMatrix3f`:  covariance matrix of the position (x,y,z) 
"""
struct TrackerHit3D <: TrackerHit
    index::ObjectID{TrackerHit3D}    # ObjectID of himself
    #---Data Members
    cellID::UInt64                   #  ID of the sensor that created this hit 
    type::Int32                      #  type of raw data hit 
    quality::Int32                   #  quality bit flag of the hit 
    time::Float32                    #  time of the hit [ns]
    eDep::Float32                    #  energy deposited on the hit [GeV]
    eDepError::Float32               #  error measured on eDep [GeV]
    position::Vector3d               #  hit position [mm]
    covMatrix::CovMatrix3f           #  covariance matrix of the position (x,y,z) 
end

function TrackerHit3D(;cellID=0, type=0, quality=0, time=0, eDep=0, eDepError=0, position=Vector3d(), covMatrix=CovMatrix3f())
    TrackerHit3D(-1, cellID, type, quality, time, eDep, eDepError, position, covMatrix)
end

"""
Vertex
- Author: EDM4hep authors
# Fields
- `type::UInt32`:  Flagword that defines the type of the vertex, see reserved bits for more information 
- `chi2::Float32`:  chi-squared of the vertex fit 
- `ndf::Int32`:  number of degrees of freedom of the vertex fit 
- `position::Vector3f`:  [mm] position of the vertex 
- `covMatrix::CovMatrix3f`:  covariance matrix of the position 
- `algorithmType::Int32`:  type code for the algorithm that has been used to create the vertex 
- `parameters::PVector{Float32}`:  additional parameters related to this vertex 
# Relations
- `particles::POD`:  particles that have been used to form this vertex, aka the decay particles emerging from this vertex 
# Methods
- `setParameters(object::Vertex, v::AbstractVector{Float32})`: assign a set of values to the `parameters` vector member
- `pushToParticles(obj::Vertex, robj::POD)`: push related object to the `particles` relation
- `popFromParticles(obj::Vertex)`: pop last related object from `particles` relation
"""
struct Vertex <: POD
    index::ObjectID{Vertex}          # ObjectID of himself
    #---Data Members
    type::UInt32                     #  Flagword that defines the type of the vertex, see reserved bits for more information 
    chi2::Float32                    #  chi-squared of the vertex fit 
    ndf::Int32                       #  number of degrees of freedom of the vertex fit 
    position::Vector3f               #  [mm] position of the vertex 
    covMatrix::CovMatrix3f           #  covariance matrix of the position 
    algorithmType::Int32             #  type code for the algorithm that has been used to create the vertex 
    #---VectorMembers
    parameters::PVector{Vertex,Float32,1}  #  additional parameters related to this vertex 
    #---OneToManyRelations
    particles::Relation{Vertex,POD,1}  #  particles that have been used to form this vertex, aka the decay particles emerging from this vertex 
end

function Vertex(;type=0, chi2=0, ndf=0, position=Vector3f(), covMatrix=CovMatrix3f(), algorithmType=0, parameters=PVector{Vertex,Float32,1}(), particles=Relation{Vertex,POD,1}())
    Vertex(-1, type, chi2, ndf, position, covMatrix, algorithmType, parameters, particles)
end

function pushToParticles(c::Vertex, o::POD)
    iszero(c.index) && (c = register(c))
    c = @set c.particles = push(c.particles, o)
    update(c)
end
function popFromParticles(c::Vertex)
    iszero(c.index) && (c = register(c))
    c = @set c.particles = pop(c.particles)
    update(c)
end
function setParameters(o::Vertex, v::AbstractVector{Float32})
    iszero(o.index) && (o = register(o))
    o = @set o.parameters = v
    update(o)
end
"""
Reconstructed track
- Author: EDM4hep authors
# Fields
- `type::Int32`:  flagword that defines the type of track 
- `chi2::Float32`:  Chi^2 of the track fit 
- `ndf::Int32`:  number of degrees of freedom of the track fit 
- `Nholes::Int32`:  number of holes on track 
- `subdetectorHitNumbers::PVector{Int32}`:  number of hits in particular subdetectors 
- `subdetectorHoleNumbers::PVector{Int32}`:  number of holes in particular subdetectors 
- `trackStates::PVector{TrackState}`:  track states 
# Relations
- `trackerHits::TrackerHit`:  hits that have been used to create this track 
- `tracks::Track`:  tracks (segments) that have been combined to create this track 
# Methods
- `setSubdetectorHitNumbers(object::Track, v::AbstractVector{Int32})`: assign a set of values to the `subdetectorHitNumbers` vector member
- `setSubdetectorHoleNumbers(object::Track, v::AbstractVector{Int32})`: assign a set of values to the `subdetectorHoleNumbers` vector member
- `setTrackStates(object::Track, v::AbstractVector{TrackState})`: assign a set of values to the `trackStates` vector member
- `pushToTrackerHits(obj::Track, robj::TrackerHit)`: push related object to the `trackerHits` relation
- `popFromTrackerHits(obj::Track)`: pop last related object from `trackerHits` relation
- `pushToTracks(obj::Track, robj::Track)`: push related object to the `tracks` relation
- `popFromTracks(obj::Track)`: pop last related object from `tracks` relation
"""
struct Track <: POD
    index::ObjectID{Track}           # ObjectID of himself
    #---Data Members
    type::Int32                      #  flagword that defines the type of track 
    chi2::Float32                    #  Chi^2 of the track fit 
    ndf::Int32                       #  number of degrees of freedom of the track fit 
    Nholes::Int32                    #  number of holes on track 
    #---VectorMembers
    subdetectorHitNumbers::PVector{Track,Int32,1}  #  number of hits in particular subdetectors 
    subdetectorHoleNumbers::PVector{Track,Int32,2}  #  number of holes in particular subdetectors 
    trackStates::PVector{Track,TrackState,3}  #  track states 
    #---OneToManyRelations
    trackerHits::Relation{Track,TrackerHit,1}  #  hits that have been used to create this track 
    tracks::Relation{Track,Track,2}  #  tracks (segments) that have been combined to create this track 
end

function Track(;type=0, chi2=0, ndf=0, Nholes=0, subdetectorHitNumbers=PVector{Track,Int32,1}(), subdetectorHoleNumbers=PVector{Track,Int32,2}(), trackStates=PVector{Track,TrackState,3}(), trackerHits=Relation{Track,TrackerHit,1}(), tracks=Relation{Track,Track,2}())
    Track(-1, type, chi2, ndf, Nholes, subdetectorHitNumbers, subdetectorHoleNumbers, trackStates, trackerHits, tracks)
end

function pushToTrackerHits(c::Track, o::TrackerHit)
    iszero(c.index) && (c = register(c))
    c = @set c.trackerHits = push(c.trackerHits, o)
    update(c)
end
function popFromTrackerHits(c::Track)
    iszero(c.index) && (c = register(c))
    c = @set c.trackerHits = pop(c.trackerHits)
    update(c)
end
function pushToTracks(c::Track, o::Track)
    iszero(c.index) && (c = register(c))
    c = @set c.tracks = push(c.tracks, o)
    update(c)
end
function popFromTracks(c::Track)
    iszero(c.index) && (c = register(c))
    c = @set c.tracks = pop(c.tracks)
    update(c)
end
function setSubdetectorHitNumbers(o::Track, v::AbstractVector{Int32})
    iszero(o.index) && (o = register(o))
    o = @set o.subdetectorHitNumbers = v
    update(o)
end
function setSubdetectorHoleNumbers(o::Track, v::AbstractVector{Int32})
    iszero(o.index) && (o = register(o))
    o = @set o.subdetectorHoleNumbers = v
    update(o)
end
function setTrackStates(o::Track, v::AbstractVector{TrackState})
    iszero(o.index) && (o = register(o))
    o = @set o.trackStates = v
    update(o)
end
"""
Reconstructed Particle
- Author: EDM4hep authors
# Fields
- `PDG::Int32`:  PDG of the reconstructed particle. 
- `energy::Float32`:  energy of the reconstructed particle. Four momentum state is not kept consistent internally [GeV]
- `momentum::Vector3f`:   particle momentum. Four momentum state is not kept consistent internally [GeV]
- `referencePoint::Vector3f`:  reference, i.e. where the particle has been measured [mm]
- `charge::Float32`:  charge of the reconstructed particle 
- `mass::Float32`:   mass of the reconstructed particle, set independently from four vector. Four momentum state is not kept consistent internally [GeV]
- `goodnessOfPID::Float32`:  overall goodness of the PID on a scale of [0;1] 
- `covMatrix::CovMatrix4f`:  covariance matrix of the reconstructed particle 4vector 
# Relations
- `decayVertex::Vertex`:  decay vertex for the particle (if it is a composite particle) 
- `clusters::Cluster`:  clusters that have been used for this particle 
- `tracks::Track`:  tracks that have been used for this particle 
- `particles::ReconstructedParticle`:  reconstructed particles that have been combined to this particle 
# Methods
- `pushToClusters(obj::ReconstructedParticle, robj::Cluster)`: push related object to the `clusters` relation
- `popFromClusters(obj::ReconstructedParticle)`: pop last related object from `clusters` relation
- `pushToTracks(obj::ReconstructedParticle, robj::Track)`: push related object to the `tracks` relation
- `popFromTracks(obj::ReconstructedParticle)`: pop last related object from `tracks` relation
- `pushToParticles(obj::ReconstructedParticle, robj::ReconstructedParticle)`: push related object to the `particles` relation
- `popFromParticles(obj::ReconstructedParticle)`: pop last related object from `particles` relation
"""
struct ReconstructedParticle <: POD
    index::ObjectID{ReconstructedParticle}  # ObjectID of himself
    #---Data Members
    PDG::Int32                       #  PDG of the reconstructed particle. 
    energy::Float32                  #  energy of the reconstructed particle. Four momentum state is not kept consistent internally [GeV]
    momentum::Vector3f               #   particle momentum. Four momentum state is not kept consistent internally [GeV]
    referencePoint::Vector3f         #  reference, i.e. where the particle has been measured [mm]
    charge::Float32                  #  charge of the reconstructed particle 
    mass::Float32                    #   mass of the reconstructed particle, set independently from four vector. Four momentum state is not kept consistent internally [GeV]
    goodnessOfPID::Float32           #  overall goodness of the PID on a scale of [0;1] 
    covMatrix::CovMatrix4f           #  covariance matrix of the reconstructed particle 4vector 
    #---OneToManyRelations
    clusters::Relation{ReconstructedParticle,Cluster,1}  #  clusters that have been used for this particle 
    tracks::Relation{ReconstructedParticle,Track,2}  #  tracks that have been used for this particle 
    particles::Relation{ReconstructedParticle,ReconstructedParticle,3}  #  reconstructed particles that have been combined to this particle 
    #---OneToOneRelations
    decayVertex_idx::ObjectID{Vertex}  #  decay vertex for the particle (if it is a composite particle) 
end

function ReconstructedParticle(;PDG=0, energy=0, momentum=Vector3f(), referencePoint=Vector3f(), charge=0, mass=0, goodnessOfPID=0, covMatrix=CovMatrix4f(), clusters=Relation{ReconstructedParticle,Cluster,1}(), tracks=Relation{ReconstructedParticle,Track,2}(), particles=Relation{ReconstructedParticle,ReconstructedParticle,3}(), decayVertex=-1)
    ReconstructedParticle(-1, PDG, energy, momentum, referencePoint, charge, mass, goodnessOfPID, covMatrix, clusters, tracks, particles, decayVertex)
end

function Base.getproperty(obj::ReconstructedParticle, sym::Symbol)
    if sym == :decayVertex
        idx = getfield(obj, :decayVertex_idx)
        return iszero(idx) ? nothing : convert(Vertex, idx)
    else # fallback to getfield
        return getfield(obj, sym)
    end
end
function pushToClusters(c::ReconstructedParticle, o::Cluster)
    iszero(c.index) && (c = register(c))
    c = @set c.clusters = push(c.clusters, o)
    update(c)
end
function popFromClusters(c::ReconstructedParticle)
    iszero(c.index) && (c = register(c))
    c = @set c.clusters = pop(c.clusters)
    update(c)
end
function pushToTracks(c::ReconstructedParticle, o::Track)
    iszero(c.index) && (c = register(c))
    c = @set c.tracks = push(c.tracks, o)
    update(c)
end
function popFromTracks(c::ReconstructedParticle)
    iszero(c.index) && (c = register(c))
    c = @set c.tracks = pop(c.tracks)
    update(c)
end
function pushToParticles(c::ReconstructedParticle, o::ReconstructedParticle)
    iszero(c.index) && (c = register(c))
    c = @set c.particles = push(c.particles, o)
    update(c)
end
function popFromParticles(c::ReconstructedParticle)
    iszero(c.index) && (c = register(c))
    c = @set c.particles = pop(c.particles)
    update(c)
end
"""
ParticleID
- Author: EDM4hep authors
# Fields
- `type::Int32`:  userdefined type 
- `PDG::Int32`:  PDG code of this id - ( 999999 ) if unknown 
- `algorithmType::Int32`:  type of the algorithm/module that created this hypothesis 
- `likelihood::Float32`:  likelihood of this hypothesis - in a user defined normalization 
- `parameters::PVector{Float32}`:  parameters associated with this hypothesis 
# Relations
- `particle::ReconstructedParticle`:  the particle from which this PID has been computed 
# Methods
- `setParameters(object::ParticleID, v::AbstractVector{Float32})`: assign a set of values to the `parameters` vector member
"""
struct ParticleID <: POD
    index::ObjectID{ParticleID}      # ObjectID of himself
    #---Data Members
    type::Int32                      #  userdefined type 
    PDG::Int32                       #  PDG code of this id - ( 999999 ) if unknown 
    algorithmType::Int32             #  type of the algorithm/module that created this hypothesis 
    likelihood::Float32              #  likelihood of this hypothesis - in a user defined normalization 
    #---VectorMembers
    parameters::PVector{ParticleID,Float32,1}  #  parameters associated with this hypothesis 
    #---OneToOneRelations
    particle_idx::ObjectID{ReconstructedParticle}  #  the particle from which this PID has been computed 
end

function ParticleID(;type=0, PDG=0, algorithmType=0, likelihood=0, parameters=PVector{ParticleID,Float32,1}(), particle=-1)
    ParticleID(-1, type, PDG, algorithmType, likelihood, parameters, particle)
end

function Base.getproperty(obj::ParticleID, sym::Symbol)
    if sym == :particle
        idx = getfield(obj, :particle_idx)
        return iszero(idx) ? nothing : convert(ReconstructedParticle, idx)
    else # fallback to getfield
        return getfield(obj, sym)
    end
end
function setParameters(o::ParticleID, v::AbstractVector{Float32})
    iszero(o.index) && (o = register(o))
    o = @set o.parameters = v
    update(o)
end
"""
dN/dx or dE/dx info of a Track
- Author: EDM4hep authors
# Fields
- `dQdx::Quantity`:  the reconstructed dEdx or dNdx and its error 
# Relations
- `track::Track`:  the corresponding track 
"""
struct RecDqdx <: POD
    index::ObjectID{RecDqdx}         # ObjectID of himself
    #---Data Members
    dQdx::Quantity                   #  the reconstructed dEdx or dNdx and its error 
    #---OneToOneRelations
    track_idx::ObjectID{Track}       #  the corresponding track 
end

function RecDqdx(;dQdx=Quantity(), track=-1)
    RecDqdx(-1, dQdx, track)
end

function Base.getproperty(obj::RecDqdx, sym::Symbol)
    if sym == :track
        idx = getfield(obj, :track_idx)
        return iszero(idx) ? nothing : convert(Track, idx)
    else # fallback to getfield
        return getfield(obj, sym)
    end
end
"""
Tracker hit plane
- Author: EDM4hep authors
# Fields
- `cellID::UInt64`:  ID of the sensor that created this hit 
- `type::Int32`:  type of raw data hit 
- `quality::Int32`:  quality bit flag of the hit 
- `time::Float32`:  time of the hit [ns]
- `eDep::Float32`:  energy deposited on the hit [GeV]
- `eDepError::Float32`:  error measured on eDep [GeV]
- `u::Vector2f`:  measurement direction vector, u lies in the x-y plane 
- `v::Vector2f`:  measurement direction vector, v is along z 
- `du::Float32`:  measurement error along the direction 
- `dv::Float32`:  measurement error along the direction 
- `position::Vector3d`:  hit position [mm]
- `covMatrix::CovMatrix3f`:  covariance of the position (x,y,z) 
"""
struct TrackerHitPlane <: TrackerHit
    index::ObjectID{TrackerHitPlane} # ObjectID of himself
    #---Data Members
    cellID::UInt64                   #  ID of the sensor that created this hit 
    type::Int32                      #  type of raw data hit 
    quality::Int32                   #  quality bit flag of the hit 
    time::Float32                    #  time of the hit [ns]
    eDep::Float32                    #  energy deposited on the hit [GeV]
    eDepError::Float32               #  error measured on eDep [GeV]
    u::Vector2f                      #  measurement direction vector, u lies in the x-y plane 
    v::Vector2f                      #  measurement direction vector, v is along z 
    du::Float32                      #  measurement error along the direction 
    dv::Float32                      #  measurement error along the direction 
    position::Vector3d               #  hit position [mm]
    covMatrix::CovMatrix3f           #  covariance of the position (x,y,z) 
end

function TrackerHitPlane(;cellID=0, type=0, quality=0, time=0, eDep=0, eDepError=0, u=Vector2f(), v=Vector2f(), du=0, dv=0, position=Vector3d(), covMatrix=CovMatrix3f())
    TrackerHitPlane(-1, cellID, type, quality, time, eDep, eDepError, u, v, du, dv, position, covMatrix)
end

"""
Simulated tracker hit
- Author: EDM4hep authors
# Fields
- `cellID::UInt64`:  ID of the sensor that created this hit 
- `eDep::Float32`:  energy deposited in the hit [GeV]
- `time::Float32`:  proper time of the hit in the lab frame [ns]
- `pathLength::Float32`:  path length of the particle in the sensitive material that resulted in this hit 
- `quality::Int32`:  quality bit flag 
- `position::Vector3d`:  the hit position [mm]
- `momentum::Vector3f`:  the 3-momentum of the particle at the hits position [GeV]
# Relations
- `particle::MCParticle`:  MCParticle that caused the hit 
"""
struct SimTrackerHit <: POD
    index::ObjectID{SimTrackerHit}   # ObjectID of himself
    #---Data Members
    cellID::UInt64                   #  ID of the sensor that created this hit 
    eDep::Float32                    #  energy deposited in the hit [GeV]
    time::Float32                    #  proper time of the hit in the lab frame [ns]
    pathLength::Float32              #  path length of the particle in the sensitive material that resulted in this hit 
    quality::Int32                   #  quality bit flag 
    position::Vector3d               #  the hit position [mm]
    momentum::Vector3f               #  the 3-momentum of the particle at the hits position [GeV]
    #---OneToOneRelations
    particle_idx::ObjectID{MCParticle}  #  MCParticle that caused the hit 
end

function SimTrackerHit(;cellID=0, eDep=0, time=0, pathLength=0, quality=0, position=Vector3d(), momentum=Vector3f(), particle=-1)
    SimTrackerHit(-1, cellID, eDep, time, pathLength, quality, position, momentum, particle)
end

function Base.getproperty(obj::SimTrackerHit, sym::Symbol)
    if sym == :particle
        idx = getfield(obj, :particle_idx)
        return iszero(idx) ? nothing : convert(MCParticle, idx)
    else # fallback to getfield
        return getfield(obj, sym)
    end
end
export setAmplitude, TimeSeries, CalorimeterHit, pushToClusters, popFromClusters, pushToHits, popFromHits, setShapeParameters, setSubdetectorEnergies, Cluster, pushToParents, popFromParents, pushToDaughters, popFromDaughters, MCParticle, GeneratorPdfInfo, CaloHitContribution, pushToContributions, popFromContributions, SimCalorimeterHit, setAdcCounts, RawTimeSeries, pushToSignalVertex, popFromSignalVertex, setCrossSections, setCrossSectionErrors, GeneratorEventParameters, setWeights, EventHeader, RawCalorimeterHit, TrackerHit3D, pushToParticles, popFromParticles, setParameters, Vertex, pushToTrackerHits, popFromTrackerHits, pushToTracks, popFromTracks, setSubdetectorHitNumbers, setSubdetectorHoleNumbers, setTrackStates, Track, ReconstructedParticle, ParticleID, RecDqdx, TrackerHitPlane, SimTrackerHit
