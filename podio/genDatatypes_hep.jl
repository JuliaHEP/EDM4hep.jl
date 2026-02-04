"""
Calibrated Detector Data
- Author: EDM4hep authors
# Fields
- `cellID::UInt64`:  cell id
- `time::Float32`:  begin time
- `interval::Float32`:  interval of each sampling
- `amplitude::PVector{Float32}`:  calibrated detector data
# Methods
- `setAmplitude(object::edm4hep!TimeSeries, v::AbstractVector{Float32})`: assign a set of values to the `amplitude` vector member
"""
struct edm4hep!TimeSeries <: POD
    index::ObjectID{edm4hep!TimeSeries}  # ObjectID of himself
    #---Data Members
    cellID::UInt64                   #  cell id
    time::Float32                    #  begin time
    interval::Float32                #  interval of each sampling
    #---VectorMembers
    amplitude::PVector{edm4hep!TimeSeries,Float32,1}  #  calibrated detector data
end

function edm4hep!TimeSeries(;cellID=zero(UInt64), time=zero(Float32), interval=zero(Float32), amplitude=PVector{edm4hep!TimeSeries,Float32,1}())
    edm4hep!TimeSeries(-1, cellID, time, interval, amplitude)
end

function setAmplitude(o::edm4hep!TimeSeries, v::AbstractVector{Float32})
    iszero(o.index) && (o = register(o))
    o = @set o.amplitude = v
    update(o)
end
"""
Calorimeter hit
- Author: EDM4hep authors
# Fields
- `cellID::UInt64`:  detector specific (geometrical) cell id
- `energy::Float32`:  energy of the hit
- `energyError::Float32`:  error of the hit energy
- `time::Float32`:  time of the hit
- `position::edm4hep!Vector3f`:  position of the hit in world coordinates
- `type::Int32`:  type of hit
"""
struct edm4hep!CalorimeterHit <: POD
    index::ObjectID{edm4hep!CalorimeterHit}  # ObjectID of himself
    #---Data Members
    cellID::UInt64                   #  detector specific (geometrical) cell id
    energy::Float32                  #  energy of the hit
    energyError::Float32             #  error of the hit energy
    time::Float32                    #  time of the hit
    position::edm4hep!Vector3f       #  position of the hit in world coordinates
    type::Int32                      #  type of hit
end

function edm4hep!CalorimeterHit(;cellID=zero(UInt64), energy=zero(Float32), energyError=zero(Float32), time=zero(Float32), position=edm4hep!Vector3f(), type=zero(Int32))
    edm4hep!CalorimeterHit(-1, cellID, energy, energyError, time, position, type)
end

"""
Calorimeter Hit Cluster
- Author: EDM4hep authors
# Fields
- `type::Int32`:  flagword that defines the type of cluster
- `energy::Float32`:  energy of the cluster
- `energyError::Float32`:  error on the energy
- `position::edm4hep!Vector3f`:  position of the cluster
- `positionError::edm4hep!CovMatrix3f`:  covariance matrix of the position
- `iTheta::Float32`:  Polar angle of the cluster's intrinsic direction (used e.g. for vertexing). Not to be confused with the cluster position seen from IP
- `iPhi::Float32`:  Azimuthal angle of the cluster's intrinsic direction (used e.g. for vertexing). Not to be confused with the cluster position seen from IP
- `directionError::edm4hep!Vector3f`:  covariance matrix of the direction
- `shapeParameters::PVector{Float32}`:  shape parameters. The corresponding names of the shape parameters should be stored in the collection named by edm4hep::labels::ShapeParameterNames in the file-level metadata, as a vector of strings in the same order as the parameters.
- `subdetectorEnergies::PVector{Float32}`:  energy observed in a particular subdetector
# Relations
- `clusters::edm4hep!Cluster`:  clusters that have been combined to this cluster
- `hits::edm4hep!CalorimeterHit`:  hits that have been combined to this cluster
# Methods
- `setShapeParameters(object::edm4hep!Cluster, v::AbstractVector{Float32})`: assign a set of values to the `shapeParameters` vector member
- `setSubdetectorEnergies(object::edm4hep!Cluster, v::AbstractVector{Float32})`: assign a set of values to the `subdetectorEnergies` vector member
- `pushToClusters(obj::edm4hep!Cluster, robj::edm4hep!Cluster)`: push related object to the `clusters` relation
- `popFromClusters(obj::edm4hep!Cluster)`: pop last related object from `clusters` relation
- `pushToHits(obj::edm4hep!Cluster, robj::edm4hep!CalorimeterHit)`: push related object to the `hits` relation
- `popFromHits(obj::edm4hep!Cluster)`: pop last related object from `hits` relation
"""
struct edm4hep!Cluster <: POD
    index::ObjectID{edm4hep!Cluster} # ObjectID of himself
    #---Data Members
    type::Int32                      #  flagword that defines the type of cluster
    energy::Float32                  #  energy of the cluster
    energyError::Float32             #  error on the energy
    position::edm4hep!Vector3f       #  position of the cluster
    positionError::edm4hep!CovMatrix3f  #  covariance matrix of the position
    iTheta::Float32                  #  Polar angle of the cluster's intrinsic direction (used e.g. for vertexing). Not to be confused with the cluster position seen from IP
    iPhi::Float32                    #  Azimuthal angle of the cluster's intrinsic direction (used e.g. for vertexing). Not to be confused with the cluster position seen from IP
    directionError::edm4hep!Vector3f #  covariance matrix of the direction
    #---VectorMembers
    shapeParameters::PVector{edm4hep!Cluster,Float32,1}  #  shape parameters. The corresponding names of the shape parameters should be stored in the collection named by edm4hep::labels::ShapeParameterNames in the file-level metadata, as a vector of strings in the same order as the parameters.
    subdetectorEnergies::PVector{edm4hep!Cluster,Float32,2}  #  energy observed in a particular subdetector
    #---OneToManyRelations
    clusters::Relation{edm4hep!Cluster,edm4hep!Cluster,1}  #  clusters that have been combined to this cluster
    hits::Relation{edm4hep!Cluster,edm4hep!CalorimeterHit,2}  #  hits that have been combined to this cluster
end

function edm4hep!Cluster(;type=zero(Int32), energy=zero(Float32), energyError=zero(Float32), position=edm4hep!Vector3f(), positionError=edm4hep!CovMatrix3f(), iTheta=zero(Float32), iPhi=zero(Float32), directionError=edm4hep!Vector3f(), shapeParameters=PVector{edm4hep!Cluster,Float32,1}(), subdetectorEnergies=PVector{edm4hep!Cluster,Float32,2}(), clusters=Relation{edm4hep!Cluster,edm4hep!Cluster,1}(), hits=Relation{edm4hep!Cluster,edm4hep!CalorimeterHit,2}())
    edm4hep!Cluster(-1, type, energy, energyError, position, positionError, iTheta, iPhi, directionError, shapeParameters, subdetectorEnergies, clusters, hits)
end

function pushToClusters(c::edm4hep!Cluster, o::edm4hep!Cluster)
    iszero(c.index) && (c = register(c))
    c = @set c.clusters = push(c.clusters, o)
    update(c)
end
function popFromClusters(c::edm4hep!Cluster)
    iszero(c.index) && (c = register(c))
    c = @set c.clusters = pop(c.clusters)
    update(c)
end
function pushToHits(c::edm4hep!Cluster, o::edm4hep!CalorimeterHit)
    iszero(c.index) && (c = register(c))
    c = @set c.hits = push(c.hits, o)
    update(c)
end
function popFromHits(c::edm4hep!Cluster)
    iszero(c.index) && (c = register(c))
    c = @set c.hits = pop(c.hits)
    update(c)
end
function setShapeParameters(o::edm4hep!Cluster, v::AbstractVector{Float32})
    iszero(o.index) && (o = register(o))
    o = @set o.shapeParameters = v
    update(o)
end
function setSubdetectorEnergies(o::edm4hep!Cluster, v::AbstractVector{Float32})
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
- `time::Float32`:  creation time of the particle in wrt. the event, e.g. for preassigned decays or decays in flight from the simulator
- `mass::Float64`:  mass of the particle
- `vertex::edm4hep!Vector3d`:  production vertex of the particle
- `endpoint::edm4hep!Vector3d`:  endpoint of the particle
- `momentum::edm4hep!Vector3d`:  particle 3-momentum at the production vertex
- `momentumAtEndpoint::edm4hep!Vector3d`:  particle 3-momentum at the endpoint
- `helicity::Int32`:  particle helicity (9 if unset)
# Relations
- `parents::edm4hep!MCParticle`:  The parents of this particle
- `daughters::edm4hep!MCParticle`:  The daughters this particle
# Methods
- `pushToParents(obj::edm4hep!MCParticle, robj::edm4hep!MCParticle)`: push related object to the `parents` relation
- `popFromParents(obj::edm4hep!MCParticle)`: pop last related object from `parents` relation
- `pushToDaughters(obj::edm4hep!MCParticle, robj::edm4hep!MCParticle)`: push related object to the `daughters` relation
- `popFromDaughters(obj::edm4hep!MCParticle)`: pop last related object from `daughters` relation
"""
struct edm4hep!MCParticle <: POD
    index::ObjectID{edm4hep!MCParticle}  # ObjectID of himself
    #---Data Members
    PDG::Int32                       #  PDG code of the particle
    generatorStatus::Int32           #  status of the particle as defined by the generator
    simulatorStatus::Int32           #  status of the particle from the simulation program - use BIT constants below
    charge::Float32                  #  particle charge
    time::Float32                    #  creation time of the particle in wrt. the event, e.g. for preassigned decays or decays in flight from the simulator
    mass::Float64                    #  mass of the particle
    vertex::edm4hep!Vector3d         #  production vertex of the particle
    endpoint::edm4hep!Vector3d       #  endpoint of the particle
    momentum::edm4hep!Vector3d       #  particle 3-momentum at the production vertex
    momentumAtEndpoint::edm4hep!Vector3d  #  particle 3-momentum at the endpoint
    helicity::Int32                  #  particle helicity (9 if unset)
    #---OneToManyRelations
    parents::Relation{edm4hep!MCParticle,edm4hep!MCParticle,1}  #  The parents of this particle
    daughters::Relation{edm4hep!MCParticle,edm4hep!MCParticle,2}  #  The daughters this particle
end

function edm4hep!MCParticle(;PDG=zero(Int32), generatorStatus=zero(Int32), simulatorStatus=zero(Int32), charge=zero(Float32), time=zero(Float32), mass=zero(Float64), vertex=edm4hep!Vector3d(), endpoint=edm4hep!Vector3d(), momentum=edm4hep!Vector3d(), momentumAtEndpoint=edm4hep!Vector3d(), helicity=9, parents=Relation{edm4hep!MCParticle,edm4hep!MCParticle,1}(), daughters=Relation{edm4hep!MCParticle,edm4hep!MCParticle,2}())
    edm4hep!MCParticle(-1, PDG, generatorStatus, simulatorStatus, charge, time, mass, vertex, endpoint, momentum, momentumAtEndpoint, helicity, parents, daughters)
end

function pushToParents(c::edm4hep!MCParticle, o::edm4hep!MCParticle)
    iszero(c.index) && (c = register(c))
    c = @set c.parents = push(c.parents, o)
    update(c)
end
function popFromParents(c::edm4hep!MCParticle)
    iszero(c.index) && (c = register(c))
    c = @set c.parents = pop(c.parents)
    update(c)
end
function pushToDaughters(c::edm4hep!MCParticle, o::edm4hep!MCParticle)
    iszero(c.index) && (c = register(c))
    c = @set c.daughters = push(c.daughters, o)
    update(c)
end
function popFromDaughters(c::edm4hep!MCParticle)
    iszero(c.index) && (c = register(c))
    c = @set c.daughters = pop(c.daughters)
    update(c)
end
"""
Sense wire hit, knowing only the distance to the wire. The circle representing possible positions is parametrized with its center, radius and normal vector (given by the wire direction).
- Author: EDM4hep authors
# Fields
- `cellID::UInt64`:  ID of the sensor that created this hit
- `type::Int32`:  type of the raw data hit
- `quality::Int32`:  quality bit flag of the hit
- `time::Float32`:  time of the hit
- `eDep::Float32`:  energy deposited by the hit
- `eDepError::Float32`:  error on eDep
- `wireStereoAngle::Float32`:  angle between the sense wire axis and the drift chamber axis (usually the z-axis) - use it together with wireAzimuthalAngle to get the wire direction
- `wireAzimuthalAngle::Float32`:  azimuthal angle at the middle of the sense wire - use it together with wireStereoAngle to get the wire direction
- `position::edm4hep!Vector3d`:  point on the sense wire which is closest to the hit (center of the circle)
- `positionAlongWireError::Float64`:  error on the hit position along the wire direction
- `distanceToWire::Float32`:  distance between the hit and the wire (radius of the circle)
- `distanceToWireError::Float32`:  error on distanceToWire
- `nElectrons::PVector{UInt16}`:  number of electrons for each cluster (number of clusters = vector size)
# Methods
- `setNElectrons(object::edm4hep!SenseWireHit, v::AbstractVector{UInt16})`: assign a set of values to the `nElectrons` vector member
"""
struct edm4hep!SenseWireHit <: edm4hep!TrackerHit
    index::ObjectID{edm4hep!SenseWireHit}  # ObjectID of himself
    #---Data Members
    cellID::UInt64                   #  ID of the sensor that created this hit
    type::Int32                      #  type of the raw data hit
    quality::Int32                   #  quality bit flag of the hit
    time::Float32                    #  time of the hit
    eDep::Float32                    #  energy deposited by the hit
    eDepError::Float32               #  error on eDep
    wireStereoAngle::Float32         #  angle between the sense wire axis and the drift chamber axis (usually the z-axis) - use it together with wireAzimuthalAngle to get the wire direction
    wireAzimuthalAngle::Float32      #  azimuthal angle at the middle of the sense wire - use it together with wireStereoAngle to get the wire direction
    position::edm4hep!Vector3d       #  point on the sense wire which is closest to the hit (center of the circle)
    positionAlongWireError::Float64  #  error on the hit position along the wire direction
    distanceToWire::Float32          #  distance between the hit and the wire (radius of the circle)
    distanceToWireError::Float32     #  error on distanceToWire
    #---VectorMembers
    nElectrons::PVector{edm4hep!SenseWireHit,UInt16,1}  #  number of electrons for each cluster (number of clusters = vector size)
end

function edm4hep!SenseWireHit(;cellID=zero(UInt64), type=zero(Int32), quality=zero(Int32), time=zero(Float32), eDep=zero(Float32), eDepError=zero(Float32), wireStereoAngle=zero(Float32), wireAzimuthalAngle=zero(Float32), position=edm4hep!Vector3d(), positionAlongWireError=zero(Float64), distanceToWire=zero(Float32), distanceToWireError=zero(Float32), nElectrons=PVector{edm4hep!SenseWireHit,UInt16,1}())
    edm4hep!SenseWireHit(-1, cellID, type, quality, time, eDep, eDepError, wireStereoAngle, wireAzimuthalAngle, position, positionAlongWireError, distanceToWire, distanceToWireError, nElectrons)
end

function setNElectrons(o::edm4hep!SenseWireHit, v::AbstractVector{UInt16})
    iszero(o.index) && (o = register(o))
    o = @set o.nElectrons = v
    update(o)
end
"""
Monte Carlo contribution to SimCalorimeterHit
- Author: EDM4hep authors
# Fields
- `PDG::Int32`:  PDG code of the shower particle that caused this contribution
- `energy::Float32`:  energy of this contribution
- `time::Float32`:  time of this contribution
- `stepPosition::edm4hep!Vector3f`:  position of this energy deposition (step)
- `stepLength::Float32`:  Geant4 step length for this contribution
# Relations
- `particle::edm4hep!MCParticle`:  MCParticle responsible for this contribution to the hit. Only particles that are kept in the MCParticle record will appear here. Hence, this will point to the first mother appearing in the record.
"""
struct edm4hep!CaloHitContribution <: POD
    index::ObjectID{edm4hep!CaloHitContribution}  # ObjectID of himself
    #---Data Members
    PDG::Int32                       #  PDG code of the shower particle that caused this contribution
    energy::Float32                  #  energy of this contribution
    time::Float32                    #  time of this contribution
    stepPosition::edm4hep!Vector3f   #  position of this energy deposition (step)
    stepLength::Float32              #  Geant4 step length for this contribution
    #---OneToOneRelations
    particle_idx::ObjectID{edm4hep!MCParticle}  #  MCParticle responsible for this contribution to the hit. Only particles that are kept in the MCParticle record will appear here. Hence, this will point to the first mother appearing in the record.
end

function edm4hep!CaloHitContribution(;PDG=zero(Int32), energy=zero(Float32), time=zero(Float32), stepPosition=edm4hep!Vector3f(), stepLength=zero(Float32), particle=-1)
    edm4hep!CaloHitContribution(-1, PDG, energy, time, stepPosition, stepLength, particle)
end

function Base.getproperty(obj::edm4hep!CaloHitContribution, sym::Symbol)
    if sym == :particle
        idx = getfield(obj, :particle_idx)
        return iszero(idx) ? nothing : convert(edm4hep!MCParticle, idx)
    else # fallback to getfield
        return getfield(obj, sym)
    end
end
"""
Simulated calorimeter hit
- Author: EDM4hep authors
# Fields
- `cellID::UInt64`:  ID of the sensor that created this hit
- `energy::Float32`:  energy of the hit
- `position::edm4hep!Vector3f`:  position of the hit in world coordinates
# Relations
- `contributions::edm4hep!CaloHitContribution`:  Monte Carlo step contributions
# Methods
- `pushToContributions(obj::edm4hep!SimCalorimeterHit, robj::edm4hep!CaloHitContribution)`: push related object to the `contributions` relation
- `popFromContributions(obj::edm4hep!SimCalorimeterHit)`: pop last related object from `contributions` relation
"""
struct edm4hep!SimCalorimeterHit <: POD
    index::ObjectID{edm4hep!SimCalorimeterHit}  # ObjectID of himself
    #---Data Members
    cellID::UInt64                   #  ID of the sensor that created this hit
    energy::Float32                  #  energy of the hit
    position::edm4hep!Vector3f       #  position of the hit in world coordinates
    #---OneToManyRelations
    contributions::Relation{edm4hep!SimCalorimeterHit,edm4hep!CaloHitContribution,1}  #  Monte Carlo step contributions
end

function edm4hep!SimCalorimeterHit(;cellID=zero(UInt64), energy=zero(Float32), position=edm4hep!Vector3f(), contributions=Relation{edm4hep!SimCalorimeterHit,edm4hep!CaloHitContribution,1}())
    edm4hep!SimCalorimeterHit(-1, cellID, energy, position, contributions)
end

function pushToContributions(c::edm4hep!SimCalorimeterHit, o::edm4hep!CaloHitContribution)
    iszero(c.index) && (c = register(c))
    c = @set c.contributions = push(c.contributions, o)
    update(c)
end
function popFromContributions(c::edm4hep!SimCalorimeterHit)
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
- `time::Float32`:  time of the hit
- `charge::Float32`:  integrated charge of the hit
- `interval::Float32`:  interval of each sampling
- `adcCounts::PVector{Int32}`:  raw data (32-bit) word at i
# Methods
- `setAdcCounts(object::edm4hep!RawTimeSeries, v::AbstractVector{Int32})`: assign a set of values to the `adcCounts` vector member
"""
struct edm4hep!RawTimeSeries <: POD
    index::ObjectID{edm4hep!RawTimeSeries}  # ObjectID of himself
    #---Data Members
    cellID::UInt64                   #  detector specific cell id
    quality::Int32                   #  quality flag for the hit
    time::Float32                    #  time of the hit
    charge::Float32                  #  integrated charge of the hit
    interval::Float32                #  interval of each sampling
    #---VectorMembers
    adcCounts::PVector{edm4hep!RawTimeSeries,Int32,1}  #  raw data (32-bit) word at i
end

function edm4hep!RawTimeSeries(;cellID=zero(UInt64), quality=zero(Int32), time=zero(Float32), charge=zero(Float32), interval=zero(Float32), adcCounts=PVector{edm4hep!RawTimeSeries,Int32,1}())
    edm4hep!RawTimeSeries(-1, cellID, quality, time, charge, interval, adcCounts)
end

function setAdcCounts(o::edm4hep!RawTimeSeries, v::AbstractVector{Int32})
    iszero(o.index) && (o = register(o))
    o = @set o.adcCounts = v
    update(o)
end
"""
Generator Event Parameters and information
- Author: EDM4hep authors
# Fields
- `sqrts::Float64`:  sqrt(s) - The nominal beam center of mass energy
- `beamsPz::SVector{2,Float64}`:  nominal z-momentum of the two incoming particle (beams)
- `partonIds::SVector{2,Int32}`:  PDG id of the partons undergoing the hard scatter
- `beamPolarizations::SVector{2,Float32}`:  Polarization of the incoming beam particles
- `crossSections::PVector{Float64}`:  List of cross sections
- `crossSectionErrors::PVector{Float64}`:  List of cross section errors
- `weights::PVector{Float64}`:  event weights. The corresponding names are stored using the edm4hep::labels::GeneratorWeightNames in the file level metadata.
# Relations
- `signalVertexParticles::edm4hep!MCParticle`:  List of initial state MCParticles that are the source of the hard interaction
# Methods
- `setCrossSections(object::edm4hep!GeneratorEventParameters, v::AbstractVector{Float64})`: assign a set of values to the `crossSections` vector member
- `setCrossSectionErrors(object::edm4hep!GeneratorEventParameters, v::AbstractVector{Float64})`: assign a set of values to the `crossSectionErrors` vector member
- `setWeights(object::edm4hep!GeneratorEventParameters, v::AbstractVector{Float64})`: assign a set of values to the `weights` vector member
- `pushToSignalVertexParticles(obj::edm4hep!GeneratorEventParameters, robj::edm4hep!MCParticle)`: push related object to the `signalVertexParticles` relation
- `popFromSignalVertexParticles(obj::edm4hep!GeneratorEventParameters)`: pop last related object from `signalVertexParticles` relation
"""
struct edm4hep!GeneratorEventParameters <: POD
    index::ObjectID{edm4hep!GeneratorEventParameters}  # ObjectID of himself
    #---Data Members
    sqrts::Float64                   #  sqrt(s) - The nominal beam center of mass energy
    beamsPz::SVector{2,Float64}      #  nominal z-momentum of the two incoming particle (beams)
    partonIds::SVector{2,Int32}      #  PDG id of the partons undergoing the hard scatter
    beamPolarizations::SVector{2,Float32}  #  Polarization of the incoming beam particles
    #---VectorMembers
    crossSections::PVector{edm4hep!GeneratorEventParameters,Float64,1}  #  List of cross sections
    crossSectionErrors::PVector{edm4hep!GeneratorEventParameters,Float64,2}  #  List of cross section errors
    weights::PVector{edm4hep!GeneratorEventParameters,Float64,3}  #  event weights. The corresponding names are stored using the edm4hep::labels::GeneratorWeightNames in the file level metadata.
    #---OneToManyRelations
    signalVertexParticles::Relation{edm4hep!GeneratorEventParameters,edm4hep!MCParticle,1}  #  List of initial state MCParticles that are the source of the hard interaction
end

function edm4hep!GeneratorEventParameters(;sqrts=zero(Float64), beamsPz=zero(SVector{2,Float64}), partonIds=zero(SVector{2,Int32}), beamPolarizations=zero(SVector{2,Float32}), crossSections=PVector{edm4hep!GeneratorEventParameters,Float64,1}(), crossSectionErrors=PVector{edm4hep!GeneratorEventParameters,Float64,2}(), weights=PVector{edm4hep!GeneratorEventParameters,Float64,3}(), signalVertexParticles=Relation{edm4hep!GeneratorEventParameters,edm4hep!MCParticle,1}())
    edm4hep!GeneratorEventParameters(-1, sqrts, beamsPz, partonIds, beamPolarizations, crossSections, crossSectionErrors, weights, signalVertexParticles)
end

function pushToSignalVertexParticles(c::edm4hep!GeneratorEventParameters, o::edm4hep!MCParticle)
    iszero(c.index) && (c = register(c))
    c = @set c.signalVertexParticles = push(c.signalVertexParticles, o)
    update(c)
end
function popFromSignalVertexParticles(c::edm4hep!GeneratorEventParameters)
    iszero(c.index) && (c = register(c))
    c = @set c.signalVertexParticles = pop(c.signalVertexParticles)
    update(c)
end
function setCrossSections(o::edm4hep!GeneratorEventParameters, v::AbstractVector{Float64})
    iszero(o.index) && (o = register(o))
    o = @set o.crossSections = v
    update(o)
end
function setCrossSectionErrors(o::edm4hep!GeneratorEventParameters, v::AbstractVector{Float64})
    iszero(o.index) && (o = register(o))
    o = @set o.crossSectionErrors = v
    update(o)
end
function setWeights(o::edm4hep!GeneratorEventParameters, v::AbstractVector{Float64})
    iszero(o.index) && (o = register(o))
    o = @set o.weights = v
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
- `weights::PVector{Float64}`:  event weights in case there are multiple. **NOTE that weights[0] might not be the same as weight!** The corresponding names of the event weights should be stored in the collection named by edm4hep::labels::EventWeightsNames in the file-level metadata.
# Methods
- `setWeights(object::edm4hep!EventHeader, v::AbstractVector{Float64})`: assign a set of values to the `weights` vector member
"""
struct edm4hep!EventHeader <: POD
    index::ObjectID{edm4hep!EventHeader}  # ObjectID of himself
    #---Data Members
    eventNumber::UInt64              #  event number
    runNumber::UInt32                #  run number
    timeStamp::UInt64                #  time stamp
    weight::Float64                  #  event weight
    #---VectorMembers
    weights::PVector{edm4hep!EventHeader,Float64,1}  #  event weights in case there are multiple. **NOTE that weights[0] might not be the same as weight!** The corresponding names of the event weights should be stored in the collection named by edm4hep::labels::EventWeightsNames in the file-level metadata.
end

function edm4hep!EventHeader(;eventNumber=zero(UInt64), runNumber=zero(UInt32), timeStamp=zero(UInt64), weight=zero(Float64), weights=PVector{edm4hep!EventHeader,Float64,1}())
    edm4hep!EventHeader(-1, eventNumber, runNumber, timeStamp, weight, weights)
end

function setWeights(o::edm4hep!EventHeader, v::AbstractVector{Float64})
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
struct edm4hep!RawCalorimeterHit <: POD
    index::ObjectID{edm4hep!RawCalorimeterHit}  # ObjectID of himself
    #---Data Members
    cellID::UInt64                   #  detector specific (geometrical) cell id
    amplitude::Int32                 #  amplitude of the hit in ADC counts
    timeStamp::Int32                 #  time stamp for the hit
end

function edm4hep!RawCalorimeterHit(;cellID=zero(UInt64), amplitude=zero(Int32), timeStamp=zero(Int32))
    edm4hep!RawCalorimeterHit(-1, cellID, amplitude, timeStamp)
end

"""
Tracker hit
- Author: EDM4hep authors
# Fields
- `cellID::UInt64`:  ID of the sensor that created this hit
- `type::Int32`:  type of raw data hit
- `quality::Int32`:  quality bit flag of the hit
- `time::Float32`:  time of the hit
- `eDep::Float32`:  energy deposited on the hit
- `eDepError::Float32`:  error measured on eDep
- `position::edm4hep!Vector3d`:  hit position
- `covMatrix::edm4hep!CovMatrix3f`:  covariance matrix of the position (x,y,z)
"""
struct edm4hep!TrackerHit3D <: edm4hep!TrackerHit
    index::ObjectID{edm4hep!TrackerHit3D}  # ObjectID of himself
    #---Data Members
    cellID::UInt64                   #  ID of the sensor that created this hit
    type::Int32                      #  type of raw data hit
    quality::Int32                   #  quality bit flag of the hit
    time::Float32                    #  time of the hit
    eDep::Float32                    #  energy deposited on the hit
    eDepError::Float32               #  error measured on eDep
    position::edm4hep!Vector3d       #  hit position
    covMatrix::edm4hep!CovMatrix3f   #  covariance matrix of the position (x,y,z)
end

function edm4hep!TrackerHit3D(;cellID=zero(UInt64), type=zero(Int32), quality=zero(Int32), time=zero(Float32), eDep=zero(Float32), eDepError=zero(Float32), position=edm4hep!Vector3d(), covMatrix=edm4hep!CovMatrix3f())
    edm4hep!TrackerHit3D(-1, cellID, type, quality, time, eDep, eDepError, position, covMatrix)
end

"""
Vertex
- Author: EDM4hep authors
# Fields
- `type::UInt32`:  flagword that defines the type of the vertex, see reserved bits for more information
- `chi2::Float32`:  chi-squared of the vertex fit
- `ndf::Int32`:  number of degrees of freedom of the vertex fit
- `position::edm4hep!Vector3f`:   position of the vertex
- `covMatrix::edm4hep!CovMatrix3f`:  covariance matrix of the position
- `algorithmType::Int32`:  type code for the algorithm that has been used to create the vertex
- `parameters::PVector{Float32}`:  additional parameters related to this vertex
# Relations
- `particles::edm4hep!POD`:  particles that have been used to form this vertex, aka the decay particles emerging from this vertex
# Methods
- `setParameters(object::edm4hep!Vertex, v::AbstractVector{Float32})`: assign a set of values to the `parameters` vector member
- `pushToParticles(obj::edm4hep!Vertex, robj::edm4hep!POD)`: push related object to the `particles` relation
- `popFromParticles(obj::edm4hep!Vertex)`: pop last related object from `particles` relation
"""
struct edm4hep!Vertex <: POD
    index::ObjectID{edm4hep!Vertex}  # ObjectID of himself
    #---Data Members
    type::UInt32                     #  flagword that defines the type of the vertex, see reserved bits for more information
    chi2::Float32                    #  chi-squared of the vertex fit
    ndf::Int32                       #  number of degrees of freedom of the vertex fit
    position::edm4hep!Vector3f       #   position of the vertex
    covMatrix::edm4hep!CovMatrix3f   #  covariance matrix of the position
    algorithmType::Int32             #  type code for the algorithm that has been used to create the vertex
    #---VectorMembers
    parameters::PVector{edm4hep!Vertex,Float32,1}  #  additional parameters related to this vertex
    #---OneToManyRelations
    particles::Relation{edm4hep!Vertex,edm4hep!POD,1}  #  particles that have been used to form this vertex, aka the decay particles emerging from this vertex
end

function edm4hep!Vertex(;type=zero(UInt32), chi2=zero(Float32), ndf=zero(Int32), position=edm4hep!Vector3f(), covMatrix=edm4hep!CovMatrix3f(), algorithmType=zero(Int32), parameters=PVector{edm4hep!Vertex,Float32,1}(), particles=Relation{edm4hep!Vertex,edm4hep!POD,1}())
    edm4hep!Vertex(-1, type, chi2, ndf, position, covMatrix, algorithmType, parameters, particles)
end

function pushToParticles(c::edm4hep!Vertex, o::edm4hep!POD)
    iszero(c.index) && (c = register(c))
    c = @set c.particles = push(c.particles, o)
    update(c)
end
function popFromParticles(c::edm4hep!Vertex)
    iszero(c.index) && (c = register(c))
    c = @set c.particles = pop(c.particles)
    update(c)
end
function setParameters(o::edm4hep!Vertex, v::AbstractVector{Float32})
    iszero(o.index) && (o = register(o))
    o = @set o.parameters = v
    update(o)
end
"""
Reconstructed track
- Author: EDM4hep authors
# Fields
- `type::Int32`:  flagword that defines the type of track
- `chi2::Float32`:  chi-squared of the track fit
- `ndf::Int32`:  number of degrees of freedom of the track fit
- `Nholes::Int32`:  number of holes on track
- `subdetectorHitNumbers::PVector{Int32}`:  number of hits in particular subdetectors
- `subdetectorHoleNumbers::PVector{Int32}`:  number of holes in particular subdetectors
- `trackStates::PVector{edm4hep!TrackState}`:  track states
# Relations
- `trackerHits::edm4hep!TrackerHit`:  hits that have been used to create this track
- `tracks::edm4hep!Track`:  tracks (segments) that have been combined to create this track
# Methods
- `setSubdetectorHitNumbers(object::edm4hep!Track, v::AbstractVector{Int32})`: assign a set of values to the `subdetectorHitNumbers` vector member
- `setSubdetectorHoleNumbers(object::edm4hep!Track, v::AbstractVector{Int32})`: assign a set of values to the `subdetectorHoleNumbers` vector member
- `setTrackStates(object::edm4hep!Track, v::AbstractVector{edm4hep!TrackState})`: assign a set of values to the `trackStates` vector member
- `pushToTrackerHits(obj::edm4hep!Track, robj::edm4hep!TrackerHit)`: push related object to the `trackerHits` relation
- `popFromTrackerHits(obj::edm4hep!Track)`: pop last related object from `trackerHits` relation
- `pushToTracks(obj::edm4hep!Track, robj::edm4hep!Track)`: push related object to the `tracks` relation
- `popFromTracks(obj::edm4hep!Track)`: pop last related object from `tracks` relation
"""
struct edm4hep!Track <: POD
    index::ObjectID{edm4hep!Track}   # ObjectID of himself
    #---Data Members
    type::Int32                      #  flagword that defines the type of track
    chi2::Float32                    #  chi-squared of the track fit
    ndf::Int32                       #  number of degrees of freedom of the track fit
    Nholes::Int32                    #  number of holes on track
    #---VectorMembers
    subdetectorHitNumbers::PVector{edm4hep!Track,Int32,1}  #  number of hits in particular subdetectors
    subdetectorHoleNumbers::PVector{edm4hep!Track,Int32,2}  #  number of holes in particular subdetectors
    trackStates::PVector{edm4hep!Track,edm4hep!TrackState,3}  #  track states
    #---OneToManyRelations
    trackerHits::Relation{edm4hep!Track,edm4hep!TrackerHit,1}  #  hits that have been used to create this track
    tracks::Relation{edm4hep!Track,edm4hep!Track,2}  #  tracks (segments) that have been combined to create this track
end

function edm4hep!Track(;type=zero(Int32), chi2=zero(Float32), ndf=zero(Int32), Nholes=zero(Int32), subdetectorHitNumbers=PVector{edm4hep!Track,Int32,1}(), subdetectorHoleNumbers=PVector{edm4hep!Track,Int32,2}(), trackStates=PVector{edm4hep!Track,edm4hep!TrackState,3}(), trackerHits=Relation{edm4hep!Track,edm4hep!TrackerHit,1}(), tracks=Relation{edm4hep!Track,edm4hep!Track,2}())
    edm4hep!Track(-1, type, chi2, ndf, Nholes, subdetectorHitNumbers, subdetectorHoleNumbers, trackStates, trackerHits, tracks)
end

function pushToTrackerHits(c::edm4hep!Track, o::edm4hep!TrackerHit)
    iszero(c.index) && (c = register(c))
    c = @set c.trackerHits = push(c.trackerHits, o)
    update(c)
end
function popFromTrackerHits(c::edm4hep!Track)
    iszero(c.index) && (c = register(c))
    c = @set c.trackerHits = pop(c.trackerHits)
    update(c)
end
function pushToTracks(c::edm4hep!Track, o::edm4hep!Track)
    iszero(c.index) && (c = register(c))
    c = @set c.tracks = push(c.tracks, o)
    update(c)
end
function popFromTracks(c::edm4hep!Track)
    iszero(c.index) && (c = register(c))
    c = @set c.tracks = pop(c.tracks)
    update(c)
end
function setSubdetectorHitNumbers(o::edm4hep!Track, v::AbstractVector{Int32})
    iszero(o.index) && (o = register(o))
    o = @set o.subdetectorHitNumbers = v
    update(o)
end
function setSubdetectorHoleNumbers(o::edm4hep!Track, v::AbstractVector{Int32})
    iszero(o.index) && (o = register(o))
    o = @set o.subdetectorHoleNumbers = v
    update(o)
end
function setTrackStates(o::edm4hep!Track, v::AbstractVector{edm4hep!TrackState})
    iszero(o.index) && (o = register(o))
    o = @set o.trackStates = v
    update(o)
end
"""
Reconstructed Particle
- Author: EDM4hep authors
# Fields
- `PDG::Int32`:  PDG of the reconstructed particle.
- `energy::Float32`:  energy of the reconstructed particle. Four momentum state is not kept consistent internally
- `momentum::edm4hep!Vector3f`:   particle momentum. Four momentum state is not kept consistent internally
- `referencePoint::edm4hep!Vector3f`:  reference, i.e. where the particle has been measured
- `charge::Float32`:  charge of the reconstructed particle
- `mass::Float32`:   mass of the reconstructed particle, set independently from four vector. Four momentum state is not kept consistent internally
- `goodnessOfPID::Float32`:  overall goodness of the PID on a scale of [0;1]
- `covMatrix::edm4hep!CovMatrix4f`:  covariance matrix of the reconstructed particle 4vector
# Relations
- `decayVertex::edm4hep!Vertex`:  decay vertex for the particle (if it is a composite particle)
- `clusters::edm4hep!Cluster`:  clusters that have been used for this particle
- `tracks::edm4hep!Track`:  tracks that have been used for this particle
- `particles::edm4hep!ReconstructedParticle`:  reconstructed particles that have been combined to this particle
# Methods
- `pushToClusters(obj::edm4hep!ReconstructedParticle, robj::edm4hep!Cluster)`: push related object to the `clusters` relation
- `popFromClusters(obj::edm4hep!ReconstructedParticle)`: pop last related object from `clusters` relation
- `pushToTracks(obj::edm4hep!ReconstructedParticle, robj::edm4hep!Track)`: push related object to the `tracks` relation
- `popFromTracks(obj::edm4hep!ReconstructedParticle)`: pop last related object from `tracks` relation
- `pushToParticles(obj::edm4hep!ReconstructedParticle, robj::edm4hep!ReconstructedParticle)`: push related object to the `particles` relation
- `popFromParticles(obj::edm4hep!ReconstructedParticle)`: pop last related object from `particles` relation
"""
struct edm4hep!ReconstructedParticle <: POD
    index::ObjectID{edm4hep!ReconstructedParticle}  # ObjectID of himself
    #---Data Members
    PDG::Int32                       #  PDG of the reconstructed particle.
    energy::Float32                  #  energy of the reconstructed particle. Four momentum state is not kept consistent internally
    momentum::edm4hep!Vector3f       #   particle momentum. Four momentum state is not kept consistent internally
    referencePoint::edm4hep!Vector3f #  reference, i.e. where the particle has been measured
    charge::Float32                  #  charge of the reconstructed particle
    mass::Float32                    #   mass of the reconstructed particle, set independently from four vector. Four momentum state is not kept consistent internally
    goodnessOfPID::Float32           #  overall goodness of the PID on a scale of [0;1]
    covMatrix::edm4hep!CovMatrix4f   #  covariance matrix of the reconstructed particle 4vector
    #---OneToManyRelations
    clusters::Relation{edm4hep!ReconstructedParticle,edm4hep!Cluster,1}  #  clusters that have been used for this particle
    tracks::Relation{edm4hep!ReconstructedParticle,edm4hep!Track,2}  #  tracks that have been used for this particle
    particles::Relation{edm4hep!ReconstructedParticle,edm4hep!ReconstructedParticle,3}  #  reconstructed particles that have been combined to this particle
    #---OneToOneRelations
    decayVertex_idx::ObjectID{edm4hep!Vertex}  #  decay vertex for the particle (if it is a composite particle)
end

function edm4hep!ReconstructedParticle(;PDG=zero(Int32), energy=zero(Float32), momentum=edm4hep!Vector3f(), referencePoint=edm4hep!Vector3f(), charge=zero(Float32), mass=zero(Float32), goodnessOfPID=zero(Float32), covMatrix=edm4hep!CovMatrix4f(), clusters=Relation{edm4hep!ReconstructedParticle,edm4hep!Cluster,1}(), tracks=Relation{edm4hep!ReconstructedParticle,edm4hep!Track,2}(), particles=Relation{edm4hep!ReconstructedParticle,edm4hep!ReconstructedParticle,3}(), decayVertex=-1)
    edm4hep!ReconstructedParticle(-1, PDG, energy, momentum, referencePoint, charge, mass, goodnessOfPID, covMatrix, clusters, tracks, particles, decayVertex)
end

function Base.getproperty(obj::edm4hep!ReconstructedParticle, sym::Symbol)
    if sym == :decayVertex
        idx = getfield(obj, :decayVertex_idx)
        return iszero(idx) ? nothing : convert(edm4hep!Vertex, idx)
    else # fallback to getfield
        return getfield(obj, sym)
    end
end
function pushToClusters(c::edm4hep!ReconstructedParticle, o::edm4hep!Cluster)
    iszero(c.index) && (c = register(c))
    c = @set c.clusters = push(c.clusters, o)
    update(c)
end
function popFromClusters(c::edm4hep!ReconstructedParticle)
    iszero(c.index) && (c = register(c))
    c = @set c.clusters = pop(c.clusters)
    update(c)
end
function pushToTracks(c::edm4hep!ReconstructedParticle, o::edm4hep!Track)
    iszero(c.index) && (c = register(c))
    c = @set c.tracks = push(c.tracks, o)
    update(c)
end
function popFromTracks(c::edm4hep!ReconstructedParticle)
    iszero(c.index) && (c = register(c))
    c = @set c.tracks = pop(c.tracks)
    update(c)
end
function pushToParticles(c::edm4hep!ReconstructedParticle, o::edm4hep!ReconstructedParticle)
    iszero(c.index) && (c = register(c))
    c = @set c.particles = push(c.particles, o)
    update(c)
end
function popFromParticles(c::edm4hep!ReconstructedParticle)
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
- `particle::edm4hep!ReconstructedParticle`:  the particle from which this PID has been computed
# Methods
- `setParameters(object::edm4hep!ParticleID, v::AbstractVector{Float32})`: assign a set of values to the `parameters` vector member
"""
struct edm4hep!ParticleID <: POD
    index::ObjectID{edm4hep!ParticleID}  # ObjectID of himself
    #---Data Members
    type::Int32                      #  userdefined type
    PDG::Int32                       #  PDG code of this id - ( 999999 ) if unknown
    algorithmType::Int32             #  type of the algorithm/module that created this hypothesis
    likelihood::Float32              #  likelihood of this hypothesis - in a user defined normalization
    #---VectorMembers
    parameters::PVector{edm4hep!ParticleID,Float32,1}  #  parameters associated with this hypothesis
    #---OneToOneRelations
    particle_idx::ObjectID{edm4hep!ReconstructedParticle}  #  the particle from which this PID has been computed
end

function edm4hep!ParticleID(;type=zero(Int32), PDG=zero(Int32), algorithmType=zero(Int32), likelihood=zero(Float32), parameters=PVector{edm4hep!ParticleID,Float32,1}(), particle=-1)
    edm4hep!ParticleID(-1, type, PDG, algorithmType, likelihood, parameters, particle)
end

function Base.getproperty(obj::edm4hep!ParticleID, sym::Symbol)
    if sym == :particle
        idx = getfield(obj, :particle_idx)
        return iszero(idx) ? nothing : convert(edm4hep!ReconstructedParticle, idx)
    else # fallback to getfield
        return getfield(obj, sym)
    end
end
function setParameters(o::edm4hep!ParticleID, v::AbstractVector{Float32})
    iszero(o.index) && (o = register(o))
    o = @set o.parameters = v
    update(o)
end
"""
dN/dx or dE/dx info of a Track
- Author: EDM4hep authors
# Fields
- `dQdx::edm4hep!Quantity`:  the reconstructed dEdx or dNdx and its error
# Relations
- `track::edm4hep!Track`:  the corresponding track
"""
struct edm4hep!RecDqdx <: POD
    index::ObjectID{edm4hep!RecDqdx} # ObjectID of himself
    #---Data Members
    dQdx::edm4hep!Quantity           #  the reconstructed dEdx or dNdx and its error
    #---OneToOneRelations
    track_idx::ObjectID{edm4hep!Track}  #  the corresponding track
end

function edm4hep!RecDqdx(;dQdx=edm4hep!Quantity(), track=-1)
    edm4hep!RecDqdx(-1, dQdx, track)
end

function Base.getproperty(obj::edm4hep!RecDqdx, sym::Symbol)
    if sym == :track
        idx = getfield(obj, :track_idx)
        return iszero(idx) ? nothing : convert(edm4hep!Track, idx)
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
- `time::Float32`:  time of the hit
- `eDep::Float32`:  energy deposited on the hit
- `eDepError::Float32`:  error measured on eDep
- `u::edm4hep!Vector2f`:  direction of the first measurement given as (theta, phi) in spherical coordinates
- `v::edm4hep!Vector2f`:  direction of the second measurement given as (theta, phi) in spherical coordinates
- `du::Float32`:  measurement error along the direction
- `dv::Float32`:  measurement error along the direction
- `position::edm4hep!Vector3d`:  hit position
- `covMatrix::edm4hep!CovMatrix3f`:  covariance of the position (x,y,z)
"""
struct edm4hep!TrackerHitPlane <: edm4hep!TrackerHit
    index::ObjectID{edm4hep!TrackerHitPlane}  # ObjectID of himself
    #---Data Members
    cellID::UInt64                   #  ID of the sensor that created this hit
    type::Int32                      #  type of raw data hit
    quality::Int32                   #  quality bit flag of the hit
    time::Float32                    #  time of the hit
    eDep::Float32                    #  energy deposited on the hit
    eDepError::Float32               #  error measured on eDep
    u::edm4hep!Vector2f              #  direction of the first measurement given as (theta, phi) in spherical coordinates
    v::edm4hep!Vector2f              #  direction of the second measurement given as (theta, phi) in spherical coordinates
    du::Float32                      #  measurement error along the direction
    dv::Float32                      #  measurement error along the direction
    position::edm4hep!Vector3d       #  hit position
    covMatrix::edm4hep!CovMatrix3f   #  covariance of the position (x,y,z)
end

function edm4hep!TrackerHitPlane(;cellID=zero(UInt64), type=zero(Int32), quality=zero(Int32), time=zero(Float32), eDep=zero(Float32), eDepError=zero(Float32), u=edm4hep!Vector2f(), v=edm4hep!Vector2f(), du=zero(Float32), dv=zero(Float32), position=edm4hep!Vector3d(), covMatrix=edm4hep!CovMatrix3f())
    edm4hep!TrackerHitPlane(-1, cellID, type, quality, time, eDep, eDepError, u, v, du, dv, position, covMatrix)
end

"""
Simulated tracker hit
- Author: EDM4hep authors
# Fields
- `cellID::UInt64`:  ID of the sensor that created this hit
- `eDep::Float32`:  energy deposited in the hit
- `time::Float32`:  proper time of the hit in the lab frame
- `pathLength::Float32`:  path length of the particle in the sensitive material that resulted in this hit
- `quality::Int32`:  quality bit flag
- `position::edm4hep!Vector3d`:  the hit position
- `momentum::edm4hep!Vector3f`:  the 3-momentum of the particle at the hits position
# Relations
- `particle::edm4hep!MCParticle`:  MCParticle that caused the hit
"""
struct edm4hep!SimTrackerHit <: POD
    index::ObjectID{edm4hep!SimTrackerHit}  # ObjectID of himself
    #---Data Members
    cellID::UInt64                   #  ID of the sensor that created this hit
    eDep::Float32                    #  energy deposited in the hit
    time::Float32                    #  proper time of the hit in the lab frame
    pathLength::Float32              #  path length of the particle in the sensitive material that resulted in this hit
    quality::Int32                   #  quality bit flag
    position::edm4hep!Vector3d       #  the hit position
    momentum::edm4hep!Vector3f       #  the 3-momentum of the particle at the hits position
    #---OneToOneRelations
    particle_idx::ObjectID{edm4hep!MCParticle}  #  MCParticle that caused the hit
end

function edm4hep!SimTrackerHit(;cellID=zero(UInt64), eDep=zero(Float32), time=zero(Float32), pathLength=zero(Float32), quality=zero(Int32), position=edm4hep!Vector3d(), momentum=edm4hep!Vector3f(), particle=-1)
    edm4hep!SimTrackerHit(-1, cellID, eDep, time, pathLength, quality, position, momentum, particle)
end

function Base.getproperty(obj::edm4hep!SimTrackerHit, sym::Symbol)
    if sym == :particle
        idx = getfield(obj, :particle_idx)
        return iszero(idx) ? nothing : convert(edm4hep!MCParticle, idx)
    else # fallback to getfield
        return getfield(obj, sym)
    end
end

#---Aliases for easier access
const SimTrackerHit = edm4hep!SimTrackerHit
const TrackerHitPlane = edm4hep!TrackerHitPlane
const Track = edm4hep!Track
const Vertex = edm4hep!Vertex
const TrackerHit3D = edm4hep!TrackerHit3D
const RawCalorimeterHit = edm4hep!RawCalorimeterHit
const EventHeader = edm4hep!EventHeader
const GeneratorEventParameters = edm4hep!GeneratorEventParameters
const RawTimeSeries = edm4hep!RawTimeSeries
const CaloHitContribution = edm4hep!CaloHitContribution
const SenseWireHit = edm4hep!SenseWireHit
const MCParticle = edm4hep!MCParticle
const ReconstructedParticle = edm4hep!ReconstructedParticle
const SimCalorimeterHit = edm4hep!SimCalorimeterHit
const Cluster = edm4hep!Cluster
const RecDqdx = edm4hep!RecDqdx
const CalorimeterHit = edm4hep!CalorimeterHit
const TimeSeries = edm4hep!TimeSeries
const ParticleID = edm4hep!ParticleID


#---Exports
export setAmplitude, pushToClusters, popFromClusters, pushToHits, popFromHits, setShapeParameters, setSubdetectorEnergies, pushToParents, popFromParents, pushToDaughters, popFromDaughters, setNElectrons, pushToContributions, popFromContributions, setAdcCounts, pushToSignalVertexParticles, popFromSignalVertexParticles, setCrossSections, setCrossSectionErrors, setWeights, pushToParticles, popFromParticles, setParameters, pushToTrackerHits, popFromTrackerHits, pushToTracks, popFromTracks, setSubdetectorHitNumbers, setSubdetectorHoleNumbers, setTrackStates, SimTrackerHit, TrackerHitPlane, Track, Vertex, TrackerHit3D, RawCalorimeterHit, EventHeader, GeneratorEventParameters, RawTimeSeries, CaloHitContribution, SenseWireHit, MCParticle, ReconstructedParticle, SimCalorimeterHit, Cluster, RecDqdx, CalorimeterHit, TimeSeries, ParticleID
