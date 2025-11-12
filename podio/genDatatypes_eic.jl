"""
Calorimeter hit
- Author: F.Gaede, DESY
# Fields
- `cellID::UInt64`: detector specific (geometrical) cell id.
- `energy::Float32`: energy of the hit in [GeV].
- `energyError::Float32`: error of the hit energy in [GeV].
- `time::Float32`: time of the hit in [ns].
- `position::edm4hep!Vector3f`: position of the hit in world coordinates in [mm].
- `type::Int32`: type of hit. Mapping of integer types to names via collection parameters "CalorimeterHitTypeNames" and "CalorimeterHitTypeValues".
"""
struct edm4hep!CalorimeterHit <: POD
    index::ObjectID{edm4hep!CalorimeterHit}  # ObjectID of himself
    #---Data Members
    cellID::UInt64                   # detector specific (geometrical) cell id.
    energy::Float32                  # energy of the hit in [GeV].
    energyError::Float32             # error of the hit energy in [GeV].
    time::Float32                    # time of the hit in [ns].
    position::edm4hep!Vector3f       # position of the hit in world coordinates in [mm].
    type::Int32                      # type of hit. Mapping of integer types to names via collection parameters "CalorimeterHitTypeNames" and "CalorimeterHitTypeValues".
end

function edm4hep!CalorimeterHit(;cellID=0, energy=0, energyError=0, time=0, position=edm4hep!Vector3f(), type=0)
    edm4hep!CalorimeterHit(-1, cellID, energy, energyError, time, position, type)
end

"""
The Monte Carlo particle - based on the lcio::MCParticle.
- Author: F.Gaede, DESY
# Fields
- `PDG::Int32`: PDG code of the particle
- `generatorStatus::Int32`: status of the particle as defined by the generator
- `simulatorStatus::Int32`: status of the particle from the simulation program - use BIT constants below
- `charge::Float32`: particle charge
- `time::Float32`: creation time of the particle in [ns] wrt. the event, e.g. for preassigned decays or decays in flight from the simulator.
- `mass::Float64`: mass of the particle in [GeV]
- `vertex::edm4hep!Vector3d`: production vertex of the particle in [mm].
- `endpoint::edm4hep!Vector3d`: endpoint of the particle in [mm]
- `momentum::edm4hep!Vector3d`: particle 3-momentum at the production vertex in [GeV]
- `momentumAtEndpoint::edm4hep!Vector3d`: particle 3-momentum at the endpoint in [GeV]
- `spin::edm4hep!Vector3f`: spin (helicity) vector of the particle.
- `colorFlow::edm4hep!Vector2i`: color flow as defined by the generator
# Relations
- `parents::edm4hep!MCParticle`: The parents of this particle.
- `daughters::edm4hep!MCParticle`: The daughters this particle.
# Methods
- `pushToParents(obj::edm4hep!MCParticle, robj::edm4hep!MCParticle)`: push related object to the `parents` relation
- `popFromParents(obj::edm4hep!MCParticle)`: pop last related object from `parents` relation
- `pushToDaughters(obj::edm4hep!MCParticle, robj::edm4hep!MCParticle)`: push related object to the `daughters` relation
- `popFromDaughters(obj::edm4hep!MCParticle)`: pop last related object from `daughters` relation
"""
struct edm4hep!MCParticle <: POD
    index::ObjectID{edm4hep!MCParticle}  # ObjectID of himself
    #---Data Members
    PDG::Int32                       # PDG code of the particle
    generatorStatus::Int32           # status of the particle as defined by the generator
    simulatorStatus::Int32           # status of the particle from the simulation program - use BIT constants below
    charge::Float32                  # particle charge
    time::Float32                    # creation time of the particle in [ns] wrt. the event, e.g. for preassigned decays or decays in flight from the simulator.
    mass::Float64                    # mass of the particle in [GeV]
    vertex::edm4hep!Vector3d         # production vertex of the particle in [mm].
    endpoint::edm4hep!Vector3d       # endpoint of the particle in [mm]
    momentum::edm4hep!Vector3d       # particle 3-momentum at the production vertex in [GeV]
    momentumAtEndpoint::edm4hep!Vector3d  # particle 3-momentum at the endpoint in [GeV]
    spin::edm4hep!Vector3f           # spin (helicity) vector of the particle.
    colorFlow::edm4hep!Vector2i      # color flow as defined by the generator
    #---OneToManyRelations
    parents::Relation{edm4hep!MCParticle,edm4hep!MCParticle,1}  # The parents of this particle.
    daughters::Relation{edm4hep!MCParticle,edm4hep!MCParticle,2}  # The daughters this particle.
end

function edm4hep!MCParticle(;PDG=0, generatorStatus=0, simulatorStatus=0, charge=0, time=0, mass=0, vertex=edm4hep!Vector3d(), endpoint=edm4hep!Vector3d(), momentum=edm4hep!Vector3d(), momentumAtEndpoint=edm4hep!Vector3d(), spin=edm4hep!Vector3f(), colorFlow=edm4hep!Vector2i(), parents=Relation{edm4hep!MCParticle,edm4hep!MCParticle,1}(), daughters=Relation{edm4hep!MCParticle,edm4hep!MCParticle,2}())
    edm4hep!MCParticle(-1, PDG, generatorStatus, simulatorStatus, charge, time, mass, vertex, endpoint, momentum, momentumAtEndpoint, spin, colorFlow, parents, daughters)
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
Simulated Primary Ionization
- Author: Wenxing Fang, IHEP
# Fields
- `cellID::UInt64`: cell id.
- `time::Float32`: the primary ionization's time in the lab frame [ns].
- `position::edm4hep!Vector3d`: the primary ionization's position [mm].
- `type::Int16`: type.
- `electronCellID::PVector{UInt64}`: cell id.
- `electronTime::PVector{Float32}`: the time in the lab frame [ns].
- `electronPosition::PVector{edm4hep!Vector3d}`: the position in the lab frame [mm].
- `pulseTime::PVector{Float32}`: the pulse's time in the lab frame [ns].
- `pulseAmplitude::PVector{Float32}`: the pulse's amplitude [fC].
# Relations
- `mcparticle::edm4hep!MCParticle`: the particle that caused the ionizing collisions.
# Methods
- `setElectronCellID(object::edm4hep!SimPrimaryIonizationCluster, v::AbstractVector{UInt64})`: assign a set of values to the `electronCellID` vector member
- `setElectronTime(object::edm4hep!SimPrimaryIonizationCluster, v::AbstractVector{Float32})`: assign a set of values to the `electronTime` vector member
- `setElectronPosition(object::edm4hep!SimPrimaryIonizationCluster, v::AbstractVector{edm4hep!Vector3d})`: assign a set of values to the `electronPosition` vector member
- `setPulseTime(object::edm4hep!SimPrimaryIonizationCluster, v::AbstractVector{Float32})`: assign a set of values to the `pulseTime` vector member
- `setPulseAmplitude(object::edm4hep!SimPrimaryIonizationCluster, v::AbstractVector{Float32})`: assign a set of values to the `pulseAmplitude` vector member
"""
struct edm4hep!SimPrimaryIonizationCluster <: POD
    index::ObjectID{edm4hep!SimPrimaryIonizationCluster}  # ObjectID of himself
    #---Data Members
    cellID::UInt64                   # cell id.
    time::Float32                    # the primary ionization's time in the lab frame [ns].
    position::edm4hep!Vector3d       # the primary ionization's position [mm].
    type::Int16                      # type.
    #---VectorMembers
    electronCellID::PVector{edm4hep!SimPrimaryIonizationCluster,UInt64,1}  # cell id.
    electronTime::PVector{edm4hep!SimPrimaryIonizationCluster,Float32,2}  # the time in the lab frame [ns].
    electronPosition::PVector{edm4hep!SimPrimaryIonizationCluster,edm4hep!Vector3d,3}  # the position in the lab frame [mm].
    pulseTime::PVector{edm4hep!SimPrimaryIonizationCluster,Float32,4}  # the pulse's time in the lab frame [ns].
    pulseAmplitude::PVector{edm4hep!SimPrimaryIonizationCluster,Float32,5}  # the pulse's amplitude [fC].
    #---OneToOneRelations
    mcparticle_idx::ObjectID{edm4hep!MCParticle}  # the particle that caused the ionizing collisions.
end

function edm4hep!SimPrimaryIonizationCluster(;cellID=0, time=0, position=edm4hep!Vector3d(), type=0, electronCellID=PVector{edm4hep!SimPrimaryIonizationCluster,UInt64,1}(), electronTime=PVector{edm4hep!SimPrimaryIonizationCluster,Float32,2}(), electronPosition=PVector{edm4hep!SimPrimaryIonizationCluster,edm4hep!Vector3d,3}(), pulseTime=PVector{edm4hep!SimPrimaryIonizationCluster,Float32,4}(), pulseAmplitude=PVector{edm4hep!SimPrimaryIonizationCluster,Float32,5}(), mcparticle=-1)
    edm4hep!SimPrimaryIonizationCluster(-1, cellID, time, position, type, electronCellID, electronTime, electronPosition, pulseTime, pulseAmplitude, mcparticle)
end

function Base.getproperty(obj::edm4hep!SimPrimaryIonizationCluster, sym::Symbol)
    if sym == :mcparticle
        idx = getfield(obj, :mcparticle_idx)
        return iszero(idx) ? nothing : convert(edm4hep!MCParticle, idx)
    else # fallback to getfield
        return getfield(obj, sym)
    end
end
function setElectronCellID(o::edm4hep!SimPrimaryIonizationCluster, v::AbstractVector{UInt64})
    iszero(o.index) && (o = register(o))
    o = @set o.electronCellID = v
    update(o)
end
function setElectronTime(o::edm4hep!SimPrimaryIonizationCluster, v::AbstractVector{Float32})
    iszero(o.index) && (o = register(o))
    o = @set o.electronTime = v
    update(o)
end
function setElectronPosition(o::edm4hep!SimPrimaryIonizationCluster, v::AbstractVector{edm4hep!Vector3d})
    iszero(o.index) && (o = register(o))
    o = @set o.electronPosition = v
    update(o)
end
function setPulseTime(o::edm4hep!SimPrimaryIonizationCluster, v::AbstractVector{Float32})
    iszero(o.index) && (o = register(o))
    o = @set o.pulseTime = v
    update(o)
end
function setPulseAmplitude(o::edm4hep!SimPrimaryIonizationCluster, v::AbstractVector{Float32})
    iszero(o.index) && (o = register(o))
    o = @set o.pulseAmplitude = v
    update(o)
end
"""
Raw data of a detector readout
- Author: F.Gaede, DESY
# Fields
- `cellID::UInt64`: detector specific cell id.
- `quality::Int32`: quality flag for the hit.
- `time::Float32`: time of the hit [ns].
- `charge::Float32`: integrated charge of the hit [fC].
- `interval::Float32`: interval of each sampling [ns].
- `adcCounts::PVector{Int32}`: raw data (32-bit) word at i.
# Methods
- `setAdcCounts(object::edm4hep!RawTimeSeries, v::AbstractVector{Int32})`: assign a set of values to the `adcCounts` vector member
"""
struct edm4hep!RawTimeSeries <: POD
    index::ObjectID{edm4hep!RawTimeSeries}  # ObjectID of himself
    #---Data Members
    cellID::UInt64                   # detector specific cell id.
    quality::Int32                   # quality flag for the hit.
    time::Float32                    # time of the hit [ns].
    charge::Float32                  # integrated charge of the hit [fC].
    interval::Float32                # interval of each sampling [ns].
    #---VectorMembers
    adcCounts::PVector{edm4hep!RawTimeSeries,Int32,1}  # raw data (32-bit) word at i.
end

function edm4hep!RawTimeSeries(;cellID=0, quality=0, time=0, charge=0, interval=0, adcCounts=PVector{edm4hep!RawTimeSeries,Int32,1}())
    edm4hep!RawTimeSeries(-1, cellID, quality, time, charge, interval, adcCounts)
end

function setAdcCounts(o::edm4hep!RawTimeSeries, v::AbstractVector{Int32})
    iszero(o.index) && (o = register(o))
    o = @set o.adcCounts = v
    update(o)
end
"""
ACTS Bound Track parameters
- Author: W. Armstrong, S. Joosten, J. Osborn
# Fields
- `type::Int32`: Type of track parameters (-1/seed, 0/head, ...)
- `surface::UInt64`: Surface for bound parameters (geometryID)
- `loc::edm4hep!Vector2f`: 2D location on surface
- `theta::Float32`: Track polar angle [rad]
- `phi::Float32`: Track azimuthal angle [rad]
- `qOverP::Float32`: [e/GeV]
- `time::Float32`: Track time [ns]
- `pdg::Int32`: pdg pid for these parameters
- `covariance::edm4eic!Cov6f`: Full covariance in basis [l0,l1,theta,phi,q/p,t]
"""
struct edm4eic!TrackParameters <: POD
    index::ObjectID{edm4eic!TrackParameters}  # ObjectID of himself
    #---Data Members
    type::Int32                      # Type of track parameters (-1/seed, 0/head, ...)
    surface::UInt64                  # Surface for bound parameters (geometryID)
    loc::edm4hep!Vector2f            # 2D location on surface
    theta::Float32                   # Track polar angle [rad]
    phi::Float32                     # Track azimuthal angle [rad]
    qOverP::Float32                  # [e/GeV]
    time::Float32                    # Track time [ns]
    pdg::Int32                       # pdg pid for these parameters
    covariance::edm4eic!Cov6f        # Full covariance in basis [l0,l1,theta,phi,q/p,t]
end

function edm4eic!TrackParameters(;type=0, surface=0, loc=edm4hep!Vector2f(), theta=0, phi=0, qOverP=0, time=0, pdg=0, covariance=edm4eic!Cov6f())
    edm4eic!TrackParameters(-1, type, surface, loc, theta, phi, qOverP, time, pdg, covariance)
end

"""
Tracker hit
- Author: F.Gaede, DESY
# Fields
- `cellID::UInt64`: ID of the sensor that created this hit
- `type::Int32`: type of raw data hit, either one of edm4hep::RawTimeSeries, edm4hep::SIMTRACKERHIT - see collection parameters "TrackerHitTypeNames" and "TrackerHitTypeValues".
- `quality::Int32`: quality bit flag of the hit.
- `time::Float32`: time of the hit [ns].
- `eDep::Float32`: energy deposited on the hit [GeV].
- `eDepError::Float32`: error measured on EDep [GeV].
- `position::edm4hep!Vector3d`: hit position in [mm].
- `covMatrix::SVector{6,Float32}`: covariance of the position (x,y,z), stored as lower triangle matrix. i.e. cov(x,x) , cov(y,x) , cov(y,y) , cov(z,x) , cov(z,y) , cov(z,z)
- `rawHits::PVector{edm4hep!ObjectID}`: raw data hits. Check getType to get actual data type.
# Methods
- `setRawHits(object::edm4hep!TrackerHit, v::AbstractVector{edm4hep!ObjectID})`: assign a set of values to the `rawHits` vector member
"""
struct edm4hep!TrackerHit <: POD
    index::ObjectID{edm4hep!TrackerHit}  # ObjectID of himself
    #---Data Members
    cellID::UInt64                   # ID of the sensor that created this hit
    type::Int32                      # type of raw data hit, either one of edm4hep::RawTimeSeries, edm4hep::SIMTRACKERHIT - see collection parameters "TrackerHitTypeNames" and "TrackerHitTypeValues".
    quality::Int32                   # quality bit flag of the hit.
    time::Float32                    # time of the hit [ns].
    eDep::Float32                    # energy deposited on the hit [GeV].
    eDepError::Float32               # error measured on EDep [GeV].
    position::edm4hep!Vector3d       # hit position in [mm].
    covMatrix::SVector{6,Float32}    # covariance of the position (x,y,z), stored as lower triangle matrix. i.e. cov(x,x) , cov(y,x) , cov(y,y) , cov(z,x) , cov(z,y) , cov(z,z)
    #---VectorMembers
    rawHits::PVector{edm4hep!TrackerHit,edm4hep!ObjectID,1}  # raw data hits. Check getType to get actual data type.
end

function edm4hep!TrackerHit(;cellID=0, type=0, quality=0, time=0, eDep=0, eDepError=0, position=edm4hep!Vector3d(), covMatrix=zero(SVector{6,Float32}), rawHits=PVector{edm4hep!TrackerHit,edm4hep!ObjectID,1}())
    edm4hep!TrackerHit(-1, cellID, type, quality, time, eDep, eDepError, position, covMatrix, rawHits)
end

function setRawHits(o::edm4hep!TrackerHit, v::AbstractVector{edm4hep!ObjectID})
    iszero(o.index) && (o = register(o))
    o = @set o.rawHits = v
    update(o)
end
"""
Raw calorimeter hit
- Author: F.Gaede, DESY
# Fields
- `cellID::UInt64`: detector specific (geometrical) cell id.
- `amplitude::Int32`: amplitude of the hit in ADC counts.
- `timeStamp::Int32`: time stamp for the hit.
"""
struct edm4hep!RawCalorimeterHit <: POD
    index::ObjectID{edm4hep!RawCalorimeterHit}  # ObjectID of himself
    #---Data Members
    cellID::UInt64                   # detector specific (geometrical) cell id.
    amplitude::Int32                 # amplitude of the hit in ADC counts.
    timeStamp::Int32                 # time stamp for the hit.
end

function edm4hep!RawCalorimeterHit(;cellID=0, amplitude=0, timeStamp=0)
    edm4hep!RawCalorimeterHit(-1, cellID, amplitude, timeStamp)
end

"""
Tracker hit plane
- Author: Placido Fernandez Declara, CERN
# Fields
- `cellID::UInt64`: ID of the sensor that created this hit
- `type::Int32`: type of raw data hit, either one of edm4hep::RawTimeSeries, edm4hep::SIMTRACKERHIT - see collection parameters "TrackerHitTypeNames" and "TrackerHitTypeValues".
- `quality::Int32`: quality bit flag of the hit.
- `time::Float32`: time of the hit [ns].
- `eDep::Float32`: energy deposited on the hit [GeV].
- `eDepError::Float32`: error measured on EDep [GeV].
- `u::edm4hep!Vector2f`: measurement direction vector, u lies in the x-y plane
- `v::edm4hep!Vector2f`: measurement direction vector, v is along z
- `du::Float32`: measurement error along the direction
- `dv::Float32`: measurement error along the direction
- `position::edm4hep!Vector3d`: hit position in [mm].
- `covMatrix::SVector{6,Float32}`: covariance of the position (x,y,z), stored as lower triangle matrix. i.e. cov(x,x) , cov(y,x) , cov(y,y) , cov(z,x) , cov(z,y) , cov(z,z)
- `rawHits::PVector{edm4hep!ObjectID}`: raw data hits. Check getType to get actual data type.
# Methods
- `setRawHits(object::edm4hep!TrackerHitPlane, v::AbstractVector{edm4hep!ObjectID})`: assign a set of values to the `rawHits` vector member
"""
struct edm4hep!TrackerHitPlane <: POD
    index::ObjectID{edm4hep!TrackerHitPlane}  # ObjectID of himself
    #---Data Members
    cellID::UInt64                   # ID of the sensor that created this hit
    type::Int32                      # type of raw data hit, either one of edm4hep::RawTimeSeries, edm4hep::SIMTRACKERHIT - see collection parameters "TrackerHitTypeNames" and "TrackerHitTypeValues".
    quality::Int32                   # quality bit flag of the hit.
    time::Float32                    # time of the hit [ns].
    eDep::Float32                    # energy deposited on the hit [GeV].
    eDepError::Float32               # error measured on EDep [GeV].
    u::edm4hep!Vector2f              # measurement direction vector, u lies in the x-y plane
    v::edm4hep!Vector2f              # measurement direction vector, v is along z
    du::Float32                      # measurement error along the direction
    dv::Float32                      # measurement error along the direction
    position::edm4hep!Vector3d       # hit position in [mm].
    covMatrix::SVector{6,Float32}    # covariance of the position (x,y,z), stored as lower triangle matrix. i.e. cov(x,x) , cov(y,x) , cov(y,y) , cov(z,x) , cov(z,y) , cov(z,z)
    #---VectorMembers
    rawHits::PVector{edm4hep!TrackerHitPlane,edm4hep!ObjectID,1}  # raw data hits. Check getType to get actual data type.
end

function edm4hep!TrackerHitPlane(;cellID=0, type=0, quality=0, time=0, eDep=0, eDepError=0, u=edm4hep!Vector2f(), v=edm4hep!Vector2f(), du=0, dv=0, position=edm4hep!Vector3d(), covMatrix=zero(SVector{6,Float32}), rawHits=PVector{edm4hep!TrackerHitPlane,edm4hep!ObjectID,1}())
    edm4hep!TrackerHitPlane(-1, cellID, type, quality, time, eDep, eDepError, u, v, du, dv, position, covMatrix, rawHits)
end

function setRawHits(o::edm4hep!TrackerHitPlane, v::AbstractVector{edm4hep!ObjectID})
    iszero(o.index) && (o = register(o))
    o = @set o.rawHits = v
    update(o)
end
"""
Simulated tracker hit
- Author: F.Gaede, DESY
# Fields
- `cellID::UInt64`: ID of the sensor that created this hit
- `EDep::Float32`: energy deposited in the hit [GeV].
- `time::Float32`: proper time of the hit in the lab frame in [ns].
- `pathLength::Float32`: path length of the particle in the sensitive material that resulted in this hit.
- `quality::Int32`: quality bit flag.
- `position::edm4hep!Vector3d`: the hit position in [mm].
- `momentum::edm4hep!Vector3f`: the 3-momentum of the particle at the hits position in [GeV]
# Relations
- `mcparticle::edm4hep!MCParticle`: MCParticle that caused the hit.
"""
struct edm4hep!SimTrackerHit <: POD
    index::ObjectID{edm4hep!SimTrackerHit}  # ObjectID of himself
    #---Data Members
    cellID::UInt64                   # ID of the sensor that created this hit
    EDep::Float32                    # energy deposited in the hit [GeV].
    time::Float32                    # proper time of the hit in the lab frame in [ns].
    pathLength::Float32              # path length of the particle in the sensitive material that resulted in this hit.
    quality::Int32                   # quality bit flag.
    position::edm4hep!Vector3d       # the hit position in [mm].
    momentum::edm4hep!Vector3f       # the 3-momentum of the particle at the hits position in [GeV]
    #---OneToOneRelations
    mcparticle_idx::ObjectID{edm4hep!MCParticle}  # MCParticle that caused the hit.
end

function edm4hep!SimTrackerHit(;cellID=0, EDep=0, time=0, pathLength=0, quality=0, position=edm4hep!Vector3d(), momentum=edm4hep!Vector3f(), mcparticle=-1)
    edm4hep!SimTrackerHit(-1, cellID, EDep, time, pathLength, quality, position, momentum, mcparticle)
end

function Base.getproperty(obj::edm4hep!SimTrackerHit, sym::Symbol)
    if sym == :mcparticle
        idx = getfield(obj, :mcparticle_idx)
        return iszero(idx) ? nothing : convert(edm4hep!MCParticle, idx)
    else # fallback to getfield
        return getfield(obj, sym)
    end
end
"""
EIC PMT hit
- Author: S. Joosten, C. Peng
# Fields
- `cellID::UInt64`: The detector specific (geometrical) cell id.
- `npe::Float32`: Estimated number of photo-electrons [#]
- `time::Float32`: Time [ns]
- `timeError::Float32`: Error on the time [ns]
- `position::edm4hep!Vector3f`: PMT hit position [mm]
- `dimension::edm4hep!Vector3f`: The dimension information of the pixel [mm].
- `sector::Int32`: The sector this hit occurred in
- `local_::edm4hep!Vector3f`: The local position of the hit in detector coordinates (relative to the sector) [mm] renamed from local due to reserved word
"""
struct edm4eic!PMTHit <: POD
    index::ObjectID{edm4eic!PMTHit}  # ObjectID of himself
    #---Data Members
    cellID::UInt64                   # The detector specific (geometrical) cell id.
    npe::Float32                     # Estimated number of photo-electrons [#]
    time::Float32                    # Time [ns]
    timeError::Float32               # Error on the time [ns]
    position::edm4hep!Vector3f       # PMT hit position [mm]
    dimension::edm4hep!Vector3f      # The dimension information of the pixel [mm].
    sector::Int32                    # The sector this hit occurred in
    local_::edm4hep!Vector3f         # The local position of the hit in detector coordinates (relative to the sector) [mm] renamed from local due to reserved word
end

function edm4eic!PMTHit(;cellID=0, npe=0, time=0, timeError=0, position=edm4hep!Vector3f(), dimension=edm4hep!Vector3f(), sector=0, local_=edm4hep!Vector3f())
    edm4eic!PMTHit(-1, cellID, npe, time, timeError, position, dimension, sector, local_)
end

"""
ParticleID
- Author: F.Gaede, DESY
# Fields
- `type::Int32`: userdefined type
- `PDG::Int32`: PDG code of this id - ( 999999 ) if unknown.
- `algorithmType::Int32`: type of the algorithm/module that created this hypothesis
- `likelihood::Float32`: likelihood of this hypothesis - in a user defined normalization.
- `parameters::PVector{Float32}`: parameters associated with this hypothesis. Check/set collection parameters ParameterNames_PIDAlgorithmTypeName for decoding the indices.
# Methods
- `setParameters(object::edm4hep!ParticleID, v::AbstractVector{Float32})`: assign a set of values to the `parameters` vector member
"""
struct edm4hep!ParticleID <: POD
    index::ObjectID{edm4hep!ParticleID}  # ObjectID of himself
    #---Data Members
    type::Int32                      # userdefined type
    PDG::Int32                       # PDG code of this id - ( 999999 ) if unknown.
    algorithmType::Int32             # type of the algorithm/module that created this hypothesis
    likelihood::Float32              # likelihood of this hypothesis - in a user defined normalization.
    #---VectorMembers
    parameters::PVector{edm4hep!ParticleID,Float32,1}  # parameters associated with this hypothesis. Check/set collection parameters ParameterNames_PIDAlgorithmTypeName for decoding the indices.
end

function edm4hep!ParticleID(;type=0, PDG=0, algorithmType=0, likelihood=0, parameters=PVector{edm4hep!ParticleID,Float32,1}())
    edm4hep!ParticleID(-1, type, PDG, algorithmType, likelihood, parameters)
end

function setParameters(o::edm4hep!ParticleID, v::AbstractVector{Float32})
    iszero(o.index) && (o = register(o))
    o = @set o.parameters = v
    update(o)
end
"""
Calorimeter Hit Cluster
- Author: F.Gaede, DESY
# Fields
- `type::Int32`: flagword that defines the type of cluster. Bits 16-31 are used internally.
- `energy::Float32`: energy of the cluster [GeV]
- `energyError::Float32`: error on the energy
- `position::edm4hep!Vector3f`: position of the cluster [mm]
- `positionError::SVector{6,Float32}`: covariance matrix of the position (6 Parameters)
- `iTheta::Float32`: intrinsic direction of cluster at position  Theta. Not to be confused with direction cluster is seen from IP.
- `phi::Float32`: intrinsic direction of cluster at position - Phi. Not to be confused with direction cluster is seen from IP.
- `directionError::edm4hep!Vector3f`: covariance matrix of the direction (3 Parameters) [mm^2]
- `shapeParameters::PVector{Float32}`: shape parameters - check/set collection parameter ClusterShapeParameters for size and names of parameters.
- `subdetectorEnergies::PVector{Float32}`: energy observed in a particular subdetector. Check/set collection parameter ClusterSubdetectorNames for decoding the indices of the array.
# Relations
- `clusters::edm4hep!Cluster`: clusters that have been combined to this cluster.
- `hits::edm4hep!CalorimeterHit`: hits that have been combined to this cluster.
- `particleIDs::edm4hep!ParticleID`: particle IDs (sorted by their likelihood)
# Methods
- `setShapeParameters(object::edm4hep!Cluster, v::AbstractVector{Float32})`: assign a set of values to the `shapeParameters` vector member
- `setSubdetectorEnergies(object::edm4hep!Cluster, v::AbstractVector{Float32})`: assign a set of values to the `subdetectorEnergies` vector member
- `pushToClusters(obj::edm4hep!Cluster, robj::edm4hep!Cluster)`: push related object to the `clusters` relation
- `popFromClusters(obj::edm4hep!Cluster)`: pop last related object from `clusters` relation
- `pushToHits(obj::edm4hep!Cluster, robj::edm4hep!CalorimeterHit)`: push related object to the `hits` relation
- `popFromHits(obj::edm4hep!Cluster)`: pop last related object from `hits` relation
- `pushToParticleIDs(obj::edm4hep!Cluster, robj::edm4hep!ParticleID)`: push related object to the `particleIDs` relation
- `popFromParticleIDs(obj::edm4hep!Cluster)`: pop last related object from `particleIDs` relation
"""
struct edm4hep!Cluster <: POD
    index::ObjectID{edm4hep!Cluster} # ObjectID of himself
    #---Data Members
    type::Int32                      # flagword that defines the type of cluster. Bits 16-31 are used internally.
    energy::Float32                  # energy of the cluster [GeV]
    energyError::Float32             # error on the energy
    position::edm4hep!Vector3f       # position of the cluster [mm]
    positionError::SVector{6,Float32}  # covariance matrix of the position (6 Parameters)
    iTheta::Float32                  # intrinsic direction of cluster at position  Theta. Not to be confused with direction cluster is seen from IP.
    phi::Float32                     # intrinsic direction of cluster at position - Phi. Not to be confused with direction cluster is seen from IP.
    directionError::edm4hep!Vector3f # covariance matrix of the direction (3 Parameters) [mm^2]
    #---VectorMembers
    shapeParameters::PVector{edm4hep!Cluster,Float32,1}  # shape parameters - check/set collection parameter ClusterShapeParameters for size and names of parameters.
    subdetectorEnergies::PVector{edm4hep!Cluster,Float32,2}  # energy observed in a particular subdetector. Check/set collection parameter ClusterSubdetectorNames for decoding the indices of the array.
    #---OneToManyRelations
    clusters::Relation{edm4hep!Cluster,edm4hep!Cluster,1}  # clusters that have been combined to this cluster.
    hits::Relation{edm4hep!Cluster,edm4hep!CalorimeterHit,2}  # hits that have been combined to this cluster.
    particleIDs::Relation{edm4hep!Cluster,edm4hep!ParticleID,3}  # particle IDs (sorted by their likelihood)
end

function edm4hep!Cluster(;type=0, energy=0, energyError=0, position=edm4hep!Vector3f(), positionError=zero(SVector{6,Float32}), iTheta=0, phi=0, directionError=edm4hep!Vector3f(), shapeParameters=PVector{edm4hep!Cluster,Float32,1}(), subdetectorEnergies=PVector{edm4hep!Cluster,Float32,2}(), clusters=Relation{edm4hep!Cluster,edm4hep!Cluster,1}(), hits=Relation{edm4hep!Cluster,edm4hep!CalorimeterHit,2}(), particleIDs=Relation{edm4hep!Cluster,edm4hep!ParticleID,3}())
    edm4hep!Cluster(-1, type, energy, energyError, position, positionError, iTheta, phi, directionError, shapeParameters, subdetectorEnergies, clusters, hits, particleIDs)
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
function pushToParticleIDs(c::edm4hep!Cluster, o::edm4hep!ParticleID)
    iszero(c.index) && (c = register(c))
    c = @set c.particleIDs = push(c.particleIDs, o)
    update(c)
end
function popFromParticleIDs(c::edm4hep!Cluster)
    iszero(c.index) && (c = register(c))
    c = @set c.particleIDs = pop(c.particleIDs)
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
Association between a Cluster and a MCParticle
- Author: Placido Fernandez Declara
# Fields
- `weight::Float32`: weight of this association
# Relations
- `rec::edm4hep!Cluster`: reference to the cluster
- `sim::edm4hep!MCParticle`: reference to the Monte-Carlo particle
"""
struct edm4hep!MCRecoClusterParticleAssociation <: POD
    index::ObjectID{edm4hep!MCRecoClusterParticleAssociation}  # ObjectID of himself
    #---Data Members
    weight::Float32                  # weight of this association
    #---OneToOneRelations
    rec_idx::ObjectID{edm4hep!Cluster}  # reference to the cluster
    sim_idx::ObjectID{edm4hep!MCParticle}  # reference to the Monte-Carlo particle
end

function edm4hep!MCRecoClusterParticleAssociation(;weight=0, rec=-1, sim=-1)
    edm4hep!MCRecoClusterParticleAssociation(-1, weight, rec, sim)
end

function Base.getproperty(obj::edm4hep!MCRecoClusterParticleAssociation, sym::Symbol)
    if sym == :rec
        idx = getfield(obj, :rec_idx)
        return iszero(idx) ? nothing : convert(edm4hep!Cluster, idx)
    elseif sym == :sim
        idx = getfield(obj, :sim_idx)
        return iszero(idx) ? nothing : convert(edm4hep!MCParticle, idx)
    else # fallback to getfield
        return getfield(obj, sym)
    end
end
"""
Association between a TrackerHit and the corresponding simulated TrackerHit
- Author: C. Bernet, B. Hegner
# Fields
- `weight::Float32`: weight of this association
# Relations
- `rec::edm4hep!TrackerHit`: reference to the reconstructed hit
- `sim::edm4hep!SimTrackerHit`: reference to the simulated hit
"""
struct edm4hep!MCRecoTrackerAssociation <: POD
    index::ObjectID{edm4hep!MCRecoTrackerAssociation}  # ObjectID of himself
    #---Data Members
    weight::Float32                  # weight of this association
    #---OneToOneRelations
    rec_idx::ObjectID{edm4hep!TrackerHit}  # reference to the reconstructed hit
    sim_idx::ObjectID{edm4hep!SimTrackerHit}  # reference to the simulated hit
end

function edm4hep!MCRecoTrackerAssociation(;weight=0, rec=-1, sim=-1)
    edm4hep!MCRecoTrackerAssociation(-1, weight, rec, sim)
end

function Base.getproperty(obj::edm4hep!MCRecoTrackerAssociation, sym::Symbol)
    if sym == :rec
        idx = getfield(obj, :rec_idx)
        return iszero(idx) ? nothing : convert(edm4hep!TrackerHit, idx)
    elseif sym == :sim
        idx = getfield(obj, :sim_idx)
        return iszero(idx) ? nothing : convert(edm4hep!SimTrackerHit, idx)
    else # fallback to getfield
        return getfield(obj, sym)
    end
end
"""
Calibrated Detector Data
- Author: Wenxing Fang, IHEP
# Fields
- `cellID::UInt64`: cell id.
- `time::Float32`: begin time [ns].
- `interval::Float32`: interval of each sampling [ns].
- `amplitude::PVector{Float32}`: calibrated detector data.
# Methods
- `setAmplitude(object::edm4hep!TimeSeries, v::AbstractVector{Float32})`: assign a set of values to the `amplitude` vector member
"""
struct edm4hep!TimeSeries <: POD
    index::ObjectID{edm4hep!TimeSeries}  # ObjectID of himself
    #---Data Members
    cellID::UInt64                   # cell id.
    time::Float32                    # begin time [ns].
    interval::Float32                # interval of each sampling [ns].
    #---VectorMembers
    amplitude::PVector{edm4hep!TimeSeries,Float32,1}  # calibrated detector data.
end

function edm4hep!TimeSeries(;cellID=0, time=0, interval=0, amplitude=PVector{edm4hep!TimeSeries,Float32,1}())
    edm4hep!TimeSeries(-1, cellID, time, interval, amplitude)
end

function setAmplitude(o::edm4hep!TimeSeries, v::AbstractVector{Float32})
    iszero(o.index) && (o = register(o))
    o = @set o.amplitude = v
    update(o)
end
"""
Reconstructed Tracker Pulse
- Author: Wenxing Fang, IHEP
# Fields
- `cellID::UInt64`: cell id.
- `time::Float32`: time [ns].
- `charge::Float32`: charge [fC].
- `quality::Int16`: quality.
- `covMatrix::SVector{3,Float32}`: lower triangle covariance matrix of the charge(c) and time(t) measurements.
# Relations
- `timeSeries::edm4hep!TimeSeries`: Optionally, the timeSeries that has been used to create the pulse can be stored with the pulse.
"""
struct edm4hep!TrackerPulse <: POD
    index::ObjectID{edm4hep!TrackerPulse}  # ObjectID of himself
    #---Data Members
    cellID::UInt64                   # cell id.
    time::Float32                    # time [ns].
    charge::Float32                  # charge [fC].
    quality::Int16                   # quality.
    covMatrix::SVector{3,Float32}    # lower triangle covariance matrix of the charge(c) and time(t) measurements.
    #---OneToOneRelations
    timeSeries_idx::ObjectID{edm4hep!TimeSeries}  # Optionally, the timeSeries that has been used to create the pulse can be stored with the pulse.
end

function edm4hep!TrackerPulse(;cellID=0, time=0, charge=0, quality=0, covMatrix=zero(SVector{3,Float32}), timeSeries=-1)
    edm4hep!TrackerPulse(-1, cellID, time, charge, quality, covMatrix, timeSeries)
end

function Base.getproperty(obj::edm4hep!TrackerPulse, sym::Symbol)
    if sym == :timeSeries
        idx = getfield(obj, :timeSeries_idx)
        return iszero(idx) ? nothing : convert(edm4hep!TimeSeries, idx)
    else # fallback to getfield
        return getfield(obj, sym)
    end
end
"""
Reconstructed Ionization Cluster
- Author: Wenxing Fang, IHEP
# Fields
- `cellID::UInt64`: cell id.
- `significance::Float32`: significance.
- `type::Int16`: type.
# Relations
- `trackerPulse::edm4hep!TrackerPulse`: the TrackerPulse used to create the ionization cluster.
# Methods
- `pushToTrackerPulse(obj::edm4hep!RecIonizationCluster, robj::edm4hep!TrackerPulse)`: push related object to the `trackerPulse` relation
- `popFromTrackerPulse(obj::edm4hep!RecIonizationCluster)`: pop last related object from `trackerPulse` relation
"""
struct edm4hep!RecIonizationCluster <: POD
    index::ObjectID{edm4hep!RecIonizationCluster}  # ObjectID of himself
    #---Data Members
    cellID::UInt64                   # cell id.
    significance::Float32            # significance.
    type::Int16                      # type.
    #---OneToManyRelations
    trackerPulse::Relation{edm4hep!RecIonizationCluster,edm4hep!TrackerPulse,1}  # the TrackerPulse used to create the ionization cluster.
end

function edm4hep!RecIonizationCluster(;cellID=0, significance=0, type=0, trackerPulse=Relation{edm4hep!RecIonizationCluster,edm4hep!TrackerPulse,1}())
    edm4hep!RecIonizationCluster(-1, cellID, significance, type, trackerPulse)
end

function pushToTrackerPulse(c::edm4hep!RecIonizationCluster, o::edm4hep!TrackerPulse)
    iszero(c.index) && (c = register(c))
    c = @set c.trackerPulse = push(c.trackerPulse, o)
    update(c)
end
function popFromTrackerPulse(c::edm4hep!RecIonizationCluster)
    iszero(c.index) && (c = register(c))
    c = @set c.trackerPulse = pop(c.trackerPulse)
    update(c)
end
"""
Calorimeter hit
- Author: W. Armstrong, S. Joosten
# Fields
- `cellID::UInt64`: The detector specific (geometrical) cell id.
- `energy::Float32`: The energy for this hit in [GeV].
- `energyError::Float32`: Error on energy [GeV].
- `time::Float32`: The time of the hit in [ns].
- `timeError::Float32`: Error on the time
- `position::edm4hep!Vector3f`: The global position of the hit in world coordinates [mm].
- `dimension::edm4hep!Vector3f`: The dimension information of the cell [mm].
- `sector::Int32`: Sector that this hit occurred in
- `layer::Int32`: Layer that the hit occurred in
- `local_::edm4hep!Vector3f`: The local coordinates of the hit in the detector segment [mm]. renamed from local due to reserved word
# Relations
- `rawHit::edm4hep!RawCalorimeterHit`: Related raw calorimeter hit
"""
struct edm4eic!CalorimeterHit <: POD
    index::ObjectID{edm4eic!CalorimeterHit}  # ObjectID of himself
    #---Data Members
    cellID::UInt64                   # The detector specific (geometrical) cell id.
    energy::Float32                  # The energy for this hit in [GeV].
    energyError::Float32             # Error on energy [GeV].
    time::Float32                    # The time of the hit in [ns].
    timeError::Float32               # Error on the time
    position::edm4hep!Vector3f       # The global position of the hit in world coordinates [mm].
    dimension::edm4hep!Vector3f      # The dimension information of the cell [mm].
    sector::Int32                    # Sector that this hit occurred in
    layer::Int32                     # Layer that the hit occurred in
    local_::edm4hep!Vector3f         # The local coordinates of the hit in the detector segment [mm]. renamed from local due to reserved word
    #---OneToOneRelations
    rawHit_idx::ObjectID{edm4hep!RawCalorimeterHit}  # Related raw calorimeter hit
end

function edm4eic!CalorimeterHit(;cellID=0, energy=0, energyError=0, time=0, timeError=0, position=edm4hep!Vector3f(), dimension=edm4hep!Vector3f(), sector=0, layer=0, local_=edm4hep!Vector3f(), rawHit=-1)
    edm4eic!CalorimeterHit(-1, cellID, energy, energyError, time, timeError, position, dimension, sector, layer, local_, rawHit)
end

function Base.getproperty(obj::edm4eic!CalorimeterHit, sym::Symbol)
    if sym == :rawHit
        idx = getfield(obj, :rawHit_idx)
        return iszero(idx) ? nothing : convert(edm4hep!RawCalorimeterHit, idx)
    else # fallback to getfield
        return getfield(obj, sym)
    end
end
"""
EIC hit cluster, reworked to more closely resemble EDM4hep
- Author: W. Armstrong, S. Joosten, C.Peng
# Fields
- `type::Int32`: Flag-word that defines the type of the cluster
- `energy::Float32`: Reconstructed energy of the cluster [GeV].
- `energyError::Float32`: Error on the cluster energy [GeV]
- `time::Float32`: [ns]
- `timeError::Float32`: Error on the cluster time
- `nhits::UInt32`: Number of hits in the cluster.
- `position::edm4hep!Vector3f`: Global position of the cluster [mm].
- `positionError::edm4eic!Cov3f`: Covariance matrix of the position (6 Parameters).
- `intrinsicTheta::Float32`: Intrinsic cluster propagation direction polar angle [rad]
- `intrinsicPhi::Float32`: Intrinsic cluster propagation direction azimuthal angle [rad]
- `intrinsicDirectionError::edm4eic!Cov2f`: Error on the intrinsic cluster propagation direction
- `shapeParameters::PVector{Float32}`: Should be set in metadata, for now it's a list of -- radius [mm], dispersion [mm], 2 entries for theta-phi widths [rad], 3 entries for x-y-z widths [mm].
- `hitContributions::PVector{Float32}`: Energy contributions of the hits. Runs parallel to ::hits()
- `subdetectorEnergies::PVector{Float32}`: Energies observed in each subdetector used for this cluster.
# Relations
- `clusters::edm4eic!Cluster`: Clusters that have been combined to form this cluster
- `hits::edm4eic!CalorimeterHit`: Hits that have been combined to form this cluster
- `particleIDs::edm4hep!ParticleID`: Particle IDs sorted by likelihood
# Methods
- `setShapeParameters(object::edm4eic!Cluster, v::AbstractVector{Float32})`: assign a set of values to the `shapeParameters` vector member
- `setHitContributions(object::edm4eic!Cluster, v::AbstractVector{Float32})`: assign a set of values to the `hitContributions` vector member
- `setSubdetectorEnergies(object::edm4eic!Cluster, v::AbstractVector{Float32})`: assign a set of values to the `subdetectorEnergies` vector member
- `pushToClusters(obj::edm4eic!Cluster, robj::edm4eic!Cluster)`: push related object to the `clusters` relation
- `popFromClusters(obj::edm4eic!Cluster)`: pop last related object from `clusters` relation
- `pushToHits(obj::edm4eic!Cluster, robj::edm4eic!CalorimeterHit)`: push related object to the `hits` relation
- `popFromHits(obj::edm4eic!Cluster)`: pop last related object from `hits` relation
- `pushToParticleIDs(obj::edm4eic!Cluster, robj::edm4hep!ParticleID)`: push related object to the `particleIDs` relation
- `popFromParticleIDs(obj::edm4eic!Cluster)`: pop last related object from `particleIDs` relation
"""
struct edm4eic!Cluster <: POD
    index::ObjectID{edm4eic!Cluster} # ObjectID of himself
    #---Data Members
    type::Int32                      # Flag-word that defines the type of the cluster
    energy::Float32                  # Reconstructed energy of the cluster [GeV].
    energyError::Float32             # Error on the cluster energy [GeV]
    time::Float32                    # [ns]
    timeError::Float32               # Error on the cluster time
    nhits::UInt32                    # Number of hits in the cluster.
    position::edm4hep!Vector3f       # Global position of the cluster [mm].
    positionError::edm4eic!Cov3f     # Covariance matrix of the position (6 Parameters).
    intrinsicTheta::Float32          # Intrinsic cluster propagation direction polar angle [rad]
    intrinsicPhi::Float32            # Intrinsic cluster propagation direction azimuthal angle [rad]
    intrinsicDirectionError::edm4eic!Cov2f  # Error on the intrinsic cluster propagation direction
    #---VectorMembers
    shapeParameters::PVector{edm4eic!Cluster,Float32,1}  # Should be set in metadata, for now it's a list of -- radius [mm], dispersion [mm], 2 entries for theta-phi widths [rad], 3 entries for x-y-z widths [mm].
    hitContributions::PVector{edm4eic!Cluster,Float32,2}  # Energy contributions of the hits. Runs parallel to ::hits()
    subdetectorEnergies::PVector{edm4eic!Cluster,Float32,3}  # Energies observed in each subdetector used for this cluster.
    #---OneToManyRelations
    clusters::Relation{edm4eic!Cluster,edm4eic!Cluster,1}  # Clusters that have been combined to form this cluster
    hits::Relation{edm4eic!Cluster,edm4eic!CalorimeterHit,2}  # Hits that have been combined to form this cluster
    particleIDs::Relation{edm4eic!Cluster,edm4hep!ParticleID,3}  # Particle IDs sorted by likelihood
end

function edm4eic!Cluster(;type=0, energy=0, energyError=0, time=0, timeError=0, nhits=0, position=edm4hep!Vector3f(), positionError=edm4eic!Cov3f(), intrinsicTheta=0, intrinsicPhi=0, intrinsicDirectionError=edm4eic!Cov2f(), shapeParameters=PVector{edm4eic!Cluster,Float32,1}(), hitContributions=PVector{edm4eic!Cluster,Float32,2}(), subdetectorEnergies=PVector{edm4eic!Cluster,Float32,3}(), clusters=Relation{edm4eic!Cluster,edm4eic!Cluster,1}(), hits=Relation{edm4eic!Cluster,edm4eic!CalorimeterHit,2}(), particleIDs=Relation{edm4eic!Cluster,edm4hep!ParticleID,3}())
    edm4eic!Cluster(-1, type, energy, energyError, time, timeError, nhits, position, positionError, intrinsicTheta, intrinsicPhi, intrinsicDirectionError, shapeParameters, hitContributions, subdetectorEnergies, clusters, hits, particleIDs)
end

function pushToClusters(c::edm4eic!Cluster, o::edm4eic!Cluster)
    iszero(c.index) && (c = register(c))
    c = @set c.clusters = push(c.clusters, o)
    update(c)
end
function popFromClusters(c::edm4eic!Cluster)
    iszero(c.index) && (c = register(c))
    c = @set c.clusters = pop(c.clusters)
    update(c)
end
function pushToHits(c::edm4eic!Cluster, o::edm4eic!CalorimeterHit)
    iszero(c.index) && (c = register(c))
    c = @set c.hits = push(c.hits, o)
    update(c)
end
function popFromHits(c::edm4eic!Cluster)
    iszero(c.index) && (c = register(c))
    c = @set c.hits = pop(c.hits)
    update(c)
end
function pushToParticleIDs(c::edm4eic!Cluster, o::edm4hep!ParticleID)
    iszero(c.index) && (c = register(c))
    c = @set c.particleIDs = push(c.particleIDs, o)
    update(c)
end
function popFromParticleIDs(c::edm4eic!Cluster)
    iszero(c.index) && (c = register(c))
    c = @set c.particleIDs = pop(c.particleIDs)
    update(c)
end
function setShapeParameters(o::edm4eic!Cluster, v::AbstractVector{Float32})
    iszero(o.index) && (o = register(o))
    o = @set o.shapeParameters = v
    update(o)
end
function setHitContributions(o::edm4eic!Cluster, v::AbstractVector{Float32})
    iszero(o.index) && (o = register(o))
    o = @set o.hitContributions = v
    update(o)
end
function setSubdetectorEnergies(o::edm4eic!Cluster, v::AbstractVector{Float32})
    iszero(o.index) && (o = register(o))
    o = @set o.subdetectorEnergies = v
    update(o)
end
"""
Association between a Cluster and a MCParticle
- Author: S. Joosten
# Fields
- `simID::UInt32`: Index of corresponding MCParticle (position in MCParticles array)
- `recID::UInt32`: Index of corresponding Cluster (position in Clusters array)
- `weight::Float32`: weight of this association
# Relations
- `rec::edm4eic!Cluster`: reference to the cluster
- `sim::edm4hep!MCParticle`: reference to the Monte-Carlo particle
"""
struct edm4eic!MCRecoClusterParticleAssociation <: POD
    index::ObjectID{edm4eic!MCRecoClusterParticleAssociation}  # ObjectID of himself
    #---Data Members
    simID::UInt32                    # Index of corresponding MCParticle (position in MCParticles array)
    recID::UInt32                    # Index of corresponding Cluster (position in Clusters array)
    weight::Float32                  # weight of this association
    #---OneToOneRelations
    rec_idx::ObjectID{edm4eic!Cluster}  # reference to the cluster
    sim_idx::ObjectID{edm4hep!MCParticle}  # reference to the Monte-Carlo particle
end

function edm4eic!MCRecoClusterParticleAssociation(;simID=0, recID=0, weight=0, rec=-1, sim=-1)
    edm4eic!MCRecoClusterParticleAssociation(-1, simID, recID, weight, rec, sim)
end

function Base.getproperty(obj::edm4eic!MCRecoClusterParticleAssociation, sym::Symbol)
    if sym == :rec
        idx = getfield(obj, :rec_idx)
        return iszero(idx) ? nothing : convert(edm4eic!Cluster, idx)
    elseif sym == :sim
        idx = getfield(obj, :sim_idx)
        return iszero(idx) ? nothing : convert(edm4hep!MCParticle, idx)
    else # fallback to getfield
        return getfield(obj, sym)
    end
end
"""
EIC vertex
- Author: J. Osborn
# Fields
- `type::Int32`: Type flag, to identify what type of vertex it is (e.g. primary, secondary, generated, etc.)
- `chi2::Float32`: Chi-squared of the vertex fit
- `ndf::Int32`: NDF of the vertex fit
- `position::edm4hep!Vector4f`: position [mm] + time t0 [ns] of the vertex. Time is 4th component in vector
- `positionError::edm4eic!Cov4f`: Covariance matrix of the position+time. Time is 4th component, similarly to 4vector
# Relations
- `associatedParticles::edm4hep!POD`: particles associated to this vertex.
# Methods
- `pushToAssociatedParticles(obj::edm4eic!Vertex, robj::edm4hep!POD)`: push related object to the `associatedParticles` relation
- `popFromAssociatedParticles(obj::edm4eic!Vertex)`: pop last related object from `associatedParticles` relation
"""
struct edm4eic!Vertex <: POD
    index::ObjectID{edm4eic!Vertex}  # ObjectID of himself
    #---Data Members
    type::Int32                      # Type flag, to identify what type of vertex it is (e.g. primary, secondary, generated, etc.)
    chi2::Float32                    # Chi-squared of the vertex fit
    ndf::Int32                       # NDF of the vertex fit
    position::edm4hep!Vector4f       # position [mm] + time t0 [ns] of the vertex. Time is 4th component in vector
    positionError::edm4eic!Cov4f     # Covariance matrix of the position+time. Time is 4th component, similarly to 4vector
    #---OneToManyRelations
    associatedParticles::Relation{edm4eic!Vertex,edm4hep!POD,1}  # particles associated to this vertex.
end

function edm4eic!Vertex(;type=0, chi2=0, ndf=0, position=edm4hep!Vector4f(), positionError=edm4eic!Cov4f(), associatedParticles=Relation{edm4eic!Vertex,edm4hep!POD,1}())
    edm4eic!Vertex(-1, type, chi2, ndf, position, positionError, associatedParticles)
end

function pushToAssociatedParticles(c::edm4eic!Vertex, o::edm4hep!POD)
    iszero(c.index) && (c = register(c))
    c = @set c.associatedParticles = push(c.associatedParticles, o)
    update(c)
end
function popFromAssociatedParticles(c::edm4eic!Vertex)
    iszero(c.index) && (c = register(c))
    c = @set c.associatedParticles = pop(c.associatedParticles)
    update(c)
end
"""
Association between a Vertex and a MCParticle
- Author: S. Joosten
# Fields
- `simID::UInt32`: Index of corresponding MCParticle (position in MCParticles array)
- `recID::UInt32`: Index of corresponding Vertex (position in Vertices array)
- `weight::Float32`: weight of this association
# Relations
- `rec::edm4eic!Vertex`: reference to the vertex
- `sim::edm4hep!MCParticle`: reference to the Monte-Carlo particle
"""
struct edm4eic!MCRecoVertexParticleAssociation <: POD
    index::ObjectID{edm4eic!MCRecoVertexParticleAssociation}  # ObjectID of himself
    #---Data Members
    simID::UInt32                    # Index of corresponding MCParticle (position in MCParticles array)
    recID::UInt32                    # Index of corresponding Vertex (position in Vertices array)
    weight::Float32                  # weight of this association
    #---OneToOneRelations
    rec_idx::ObjectID{edm4eic!Vertex}  # reference to the vertex
    sim_idx::ObjectID{edm4hep!MCParticle}  # reference to the Monte-Carlo particle
end

function edm4eic!MCRecoVertexParticleAssociation(;simID=0, recID=0, weight=0, rec=-1, sim=-1)
    edm4eic!MCRecoVertexParticleAssociation(-1, simID, recID, weight, rec, sim)
end

function Base.getproperty(obj::edm4eic!MCRecoVertexParticleAssociation, sym::Symbol)
    if sym == :rec
        idx = getfield(obj, :rec_idx)
        return iszero(idx) ? nothing : convert(edm4eic!Vertex, idx)
    elseif sym == :sim
        idx = getfield(obj, :sim_idx)
        return iszero(idx) ? nothing : convert(edm4hep!MCParticle, idx)
    else # fallback to getfield
        return getfield(obj, sym)
    end
end
"""
Tensor type for use in training in inference of ML models
- Author: D. Kalinkin
# Fields
- `elementType::Int32`: Data type in the same encoding as "ONNXTensorElementDataType", 1 - float, 7 - int64
- `shape::PVector{Int64}`: Vector of tensor lengths along its axes
- `floatData::PVector{Float32}`: Iff elementType==1, values are stored here
- `int64Data::PVector{Int64}`: Iff elementType==7, values are stored here
# Methods
- `setShape(object::edm4eic!Tensor, v::AbstractVector{Int64})`: assign a set of values to the `shape` vector member
- `setFloatData(object::edm4eic!Tensor, v::AbstractVector{Float32})`: assign a set of values to the `floatData` vector member
- `setInt64Data(object::edm4eic!Tensor, v::AbstractVector{Int64})`: assign a set of values to the `int64Data` vector member
"""
struct edm4eic!Tensor <: POD
    index::ObjectID{edm4eic!Tensor}  # ObjectID of himself
    #---Data Members
    elementType::Int32               # Data type in the same encoding as "ONNXTensorElementDataType", 1 - float, 7 - int64
    #---VectorMembers
    shape::PVector{edm4eic!Tensor,Int64,1}  # Vector of tensor lengths along its axes
    floatData::PVector{edm4eic!Tensor,Float32,2}  # Iff elementType==1, values are stored here
    int64Data::PVector{edm4eic!Tensor,Int64,3}  # Iff elementType==7, values are stored here
end

function edm4eic!Tensor(;elementType=0, shape=PVector{edm4eic!Tensor,Int64,1}(), floatData=PVector{edm4eic!Tensor,Float32,2}(), int64Data=PVector{edm4eic!Tensor,Int64,3}())
    edm4eic!Tensor(-1, elementType, shape, floatData, int64Data)
end

function setShape(o::edm4eic!Tensor, v::AbstractVector{Int64})
    iszero(o.index) && (o = register(o))
    o = @set o.shape = v
    update(o)
end
function setFloatData(o::edm4eic!Tensor, v::AbstractVector{Float32})
    iszero(o.index) && (o = register(o))
    o = @set o.floatData = v
    update(o)
end
function setInt64Data(o::edm4eic!Tensor, v::AbstractVector{Int64})
    iszero(o.index) && (o = register(o))
    o = @set o.int64Data = v
    update(o)
end
"""
Raw hit from an HGCROC chip
- Author: D. Anderson, S. Joosten, T. Protzman, N. Novitzky, D. Kalinkin
# Fields
- `cellID::UInt64`: Detector specific (geometrical) cell id
- `samplePhase::Int32`: Phase of samples in [# samples], for synchronizing across chips
- `timeStamp::Int32`: [TDC counts]
- `samples::PVector{edm4eic!HGCROCSample}`: ADC, Time of Arrival (TOA), and Time over Threshold (TOT) values for each sample read out
# Methods
- `setSamples(object::edm4eic!RawHGCROCHit, v::AbstractVector{edm4eic!HGCROCSample})`: assign a set of values to the `samples` vector member
"""
struct edm4eic!RawHGCROCHit <: POD
    index::ObjectID{edm4eic!RawHGCROCHit}  # ObjectID of himself
    #---Data Members
    cellID::UInt64                   # Detector specific (geometrical) cell id
    samplePhase::Int32               # Phase of samples in [# samples], for synchronizing across chips
    timeStamp::Int32                 # [TDC counts]
    #---VectorMembers
    samples::PVector{edm4eic!RawHGCROCHit,edm4eic!HGCROCSample,1}  # ADC, Time of Arrival (TOA), and Time over Threshold (TOT) values for each sample read out
end

function edm4eic!RawHGCROCHit(;cellID=0, samplePhase=0, timeStamp=0, samples=PVector{edm4eic!RawHGCROCHit,edm4eic!HGCROCSample,1}())
    edm4eic!RawHGCROCHit(-1, cellID, samplePhase, timeStamp, samples)
end

function setSamples(o::edm4eic!RawHGCROCHit, v::AbstractVector{edm4eic!HGCROCSample})
    iszero(o.index) && (o = register(o))
    o = @set o.samples = v
    update(o)
end
"""
Association between a CalorimeterHit and a MCParticle
- Author: Placido Fernandez Declara
# Fields
- `weight::Float32`: weight of this association
# Relations
- `rec::edm4hep!CalorimeterHit`: reference to the reconstructed hit
- `sim::edm4hep!MCParticle`: reference to the Monte-Carlo particle
"""
struct edm4hep!MCRecoCaloParticleAssociation <: POD
    index::ObjectID{edm4hep!MCRecoCaloParticleAssociation}  # ObjectID of himself
    #---Data Members
    weight::Float32                  # weight of this association
    #---OneToOneRelations
    rec_idx::ObjectID{edm4hep!CalorimeterHit}  # reference to the reconstructed hit
    sim_idx::ObjectID{edm4hep!MCParticle}  # reference to the Monte-Carlo particle
end

function edm4hep!MCRecoCaloParticleAssociation(;weight=0, rec=-1, sim=-1)
    edm4hep!MCRecoCaloParticleAssociation(-1, weight, rec, sim)
end

function Base.getproperty(obj::edm4hep!MCRecoCaloParticleAssociation, sym::Symbol)
    if sym == :rec
        idx = getfield(obj, :rec_idx)
        return iszero(idx) ? nothing : convert(edm4hep!CalorimeterHit, idx)
    elseif sym == :sim
        idx = getfield(obj, :sim_idx)
        return iszero(idx) ? nothing : convert(edm4hep!MCParticle, idx)
    else # fallback to getfield
        return getfield(obj, sym)
    end
end
"""
Raw (digitized) tracker hit
- Author: W. Armstrong, S. Joosten
# Fields
- `cellID::UInt64`: The detector specific (geometrical) cell id.
- `charge::Int32`: ADC value
- `timeStamp::Int32`: TDC value.
"""
struct edm4eic!RawTrackerHit <: POD
    index::ObjectID{edm4eic!RawTrackerHit}  # ObjectID of himself
    #---Data Members
    cellID::UInt64                   # The detector specific (geometrical) cell id.
    charge::Int32                    # ADC value
    timeStamp::Int32                 # TDC value.
end

function edm4eic!RawTrackerHit(;cellID=0, charge=0, timeStamp=0)
    edm4eic!RawTrackerHit(-1, cellID, charge, timeStamp)
end

"""
Association between a RawTrackerHit and a SimTrackerHit
- Author: C. Dilks, W. Deconinck
# Fields
- `weight::Float32`: weight of this association
# Relations
- `rawHit::edm4eic!RawTrackerHit`: reference to the digitized hit
- `simHit::edm4hep!SimTrackerHit`: reference to the simulated hit
"""
struct edm4eic!MCRecoTrackerHitAssociation <: POD
    index::ObjectID{edm4eic!MCRecoTrackerHitAssociation}  # ObjectID of himself
    #---Data Members
    weight::Float32                  # weight of this association
    #---OneToOneRelations
    rawHit_idx::ObjectID{edm4eic!RawTrackerHit}  # reference to the digitized hit
    simHit_idx::ObjectID{edm4hep!SimTrackerHit}  # reference to the simulated hit
end

function edm4eic!MCRecoTrackerHitAssociation(;weight=0, rawHit=-1, simHit=-1)
    edm4eic!MCRecoTrackerHitAssociation(-1, weight, rawHit, simHit)
end

function Base.getproperty(obj::edm4eic!MCRecoTrackerHitAssociation, sym::Symbol)
    if sym == :rawHit
        idx = getfield(obj, :rawHit_idx)
        return iszero(idx) ? nothing : convert(edm4eic!RawTrackerHit, idx)
    elseif sym == :simHit
        idx = getfield(obj, :simHit_idx)
        return iszero(idx) ? nothing : convert(edm4hep!SimTrackerHit, idx)
    else # fallback to getfield
        return getfield(obj, sym)
    end
end
"""
EIC Ring Image Cluster
- Author: S. Joosten, C. Peng
# Fields
- `npe::Float32`: Number of photo-electrons [#]
- `position::edm4hep!Vector3f`: Global position of the cluster [mm]
- `positionError::edm4hep!Vector3f`: Error on the position
- `theta::Float32`: Opening angle of the ring [rad, 0->pi]
- `thetaError::Float32`: Error on the opening angle
- `radius::Float32`: Radius of the best fit ring [mm]
- `radiusError::Float32`: Estimated error from the fit [mm]
"""
struct edm4eic!RingImage <: POD
    index::ObjectID{edm4eic!RingImage}  # ObjectID of himself
    #---Data Members
    npe::Float32                     # Number of photo-electrons [#]
    position::edm4hep!Vector3f       # Global position of the cluster [mm]
    positionError::edm4hep!Vector3f  # Error on the position
    theta::Float32                   # Opening angle of the ring [rad, 0->pi]
    thetaError::Float32              # Error on the opening angle
    radius::Float32                  # Radius of the best fit ring [mm]
    radiusError::Float32             # Estimated error from the fit [mm]
end

function edm4eic!RingImage(;npe=0, position=edm4hep!Vector3f(), positionError=edm4hep!Vector3f(), theta=0, thetaError=0, radius=0, radiusError=0)
    edm4eic!RingImage(-1, npe, position, positionError, theta, thetaError, radius, radiusError)
end

"""
Association between a TrackerHitPlane and the corresponding simulated TrackerHit
- Author: Placido Fernandez Declara
# Fields
- `weight::Float32`: weight of this association
# Relations
- `rec::edm4hep!TrackerHitPlane`: reference to the reconstructed hit
- `sim::edm4hep!SimTrackerHit`: reference to the simulated hit
"""
struct edm4hep!MCRecoTrackerHitPlaneAssociation <: POD
    index::ObjectID{edm4hep!MCRecoTrackerHitPlaneAssociation}  # ObjectID of himself
    #---Data Members
    weight::Float32                  # weight of this association
    #---OneToOneRelations
    rec_idx::ObjectID{edm4hep!TrackerHitPlane}  # reference to the reconstructed hit
    sim_idx::ObjectID{edm4hep!SimTrackerHit}  # reference to the simulated hit
end

function edm4hep!MCRecoTrackerHitPlaneAssociation(;weight=0, rec=-1, sim=-1)
    edm4hep!MCRecoTrackerHitPlaneAssociation(-1, weight, rec, sim)
end

function Base.getproperty(obj::edm4hep!MCRecoTrackerHitPlaneAssociation, sym::Symbol)
    if sym == :rec
        idx = getfield(obj, :rec_idx)
        return iszero(idx) ? nothing : convert(edm4hep!TrackerHitPlane, idx)
    elseif sym == :sim
        idx = getfield(obj, :sim_idx)
        return iszero(idx) ? nothing : convert(edm4hep!SimTrackerHit, idx)
    else # fallback to getfield
        return getfield(obj, sym)
    end
end
"""
Monte Carlo contribution to SimCalorimeterHit
- Author: F.Gaede, DESY
# Fields
- `PDG::Int32`: PDG code of the shower particle that caused this contribution.
- `energy::Float32`: energy in [GeV] of the this contribution
- `time::Float32`: time in [ns] of this contribution
- `stepPosition::edm4hep!Vector3f`: position of this energy deposition (step) [mm]
# Relations
- `particle::edm4hep!MCParticle`: primary MCParticle that caused the shower responsible for this contribution to the hit.
"""
struct edm4hep!CaloHitContribution <: POD
    index::ObjectID{edm4hep!CaloHitContribution}  # ObjectID of himself
    #---Data Members
    PDG::Int32                       # PDG code of the shower particle that caused this contribution.
    energy::Float32                  # energy in [GeV] of the this contribution
    time::Float32                    # time in [ns] of this contribution
    stepPosition::edm4hep!Vector3f   # position of this energy deposition (step) [mm]
    #---OneToOneRelations
    particle_idx::ObjectID{edm4hep!MCParticle}  # primary MCParticle that caused the shower responsible for this contribution to the hit.
end

function edm4hep!CaloHitContribution(;PDG=0, energy=0, time=0, stepPosition=edm4hep!Vector3f(), particle=-1)
    edm4hep!CaloHitContribution(-1, PDG, energy, time, stepPosition, particle)
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
- Author: F.Gaede, DESY
# Fields
- `cellID::UInt64`: ID of the sensor that created this hit
- `energy::Float32`: energy of the hit in [GeV].
- `position::edm4hep!Vector3f`: position of the hit in world coordinates in [mm].
# Relations
- `contributions::edm4hep!CaloHitContribution`: Monte Carlo step contribution - parallel to particle
# Methods
- `pushToContributions(obj::edm4hep!SimCalorimeterHit, robj::edm4hep!CaloHitContribution)`: push related object to the `contributions` relation
- `popFromContributions(obj::edm4hep!SimCalorimeterHit)`: pop last related object from `contributions` relation
"""
struct edm4hep!SimCalorimeterHit <: POD
    index::ObjectID{edm4hep!SimCalorimeterHit}  # ObjectID of himself
    #---Data Members
    cellID::UInt64                   # ID of the sensor that created this hit
    energy::Float32                  # energy of the hit in [GeV].
    position::edm4hep!Vector3f       # position of the hit in world coordinates in [mm].
    #---OneToManyRelations
    contributions::Relation{edm4hep!SimCalorimeterHit,edm4hep!CaloHitContribution,1}  # Monte Carlo step contribution - parallel to particle
end

function edm4hep!SimCalorimeterHit(;cellID=0, energy=0, position=edm4hep!Vector3f(), contributions=Relation{edm4hep!SimCalorimeterHit,edm4hep!CaloHitContribution,1}())
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
Association between a RawCalorimeterHit and a SimCalorimeterHit
- Author: S. Rahman
# Fields
- `weight::Float32`: weight of this association
# Relations
- `rawHit::edm4hep!RawCalorimeterHit`: reference to the digitized calorimeter hit
- `simHit::edm4hep!SimCalorimeterHit`: reference to the simulated calorimeter hit
"""
struct edm4eic!MCRecoCalorimeterHitAssociation <: POD
    index::ObjectID{edm4eic!MCRecoCalorimeterHitAssociation}  # ObjectID of himself
    #---Data Members
    weight::Float32                  # weight of this association
    #---OneToOneRelations
    rawHit_idx::ObjectID{edm4hep!RawCalorimeterHit}  # reference to the digitized calorimeter hit
    simHit_idx::ObjectID{edm4hep!SimCalorimeterHit}  # reference to the simulated calorimeter hit
end

function edm4eic!MCRecoCalorimeterHitAssociation(;weight=0, rawHit=-1, simHit=-1)
    edm4eic!MCRecoCalorimeterHitAssociation(-1, weight, rawHit, simHit)
end

function Base.getproperty(obj::edm4eic!MCRecoCalorimeterHitAssociation, sym::Symbol)
    if sym == :rawHit
        idx = getfield(obj, :rawHit_idx)
        return iszero(idx) ? nothing : convert(edm4hep!RawCalorimeterHit, idx)
    elseif sym == :simHit
        idx = getfield(obj, :simHit_idx)
        return iszero(idx) ? nothing : convert(edm4hep!SimCalorimeterHit, idx)
    else # fallback to getfield
        return getfield(obj, sym)
    end
end
"""
Simulated pulse prior to digitization.
- Author: D. Anderson, S. Gardner, S. Joosten., D. Kalinkin
# Fields
- `cellID::UInt64`: ID of the readout cell for this pulse.
- `integral::Float32`: Total pulse integral in relevant units.
- `position::edm4hep!Vector3f`: Position the pulse is evaluated in world coordinates [mm].
- `time::Float32`: Start time for the pulse in [ns].
- `interval::Float32`: Time interval between amplitude values [ns].
- `amplitude::PVector{Float32}`: Pulse amplitude in relevant units, sum of amplitude values equals integral
# Relations
- `calorimeterHits::edm4hep!SimCalorimeterHit`: SimCalorimeterHits used to create this pulse
- `trackerHits::edm4hep!SimTrackerHit`: SimTrackerHits used to create this pulse
- `pulses::edm4eic!SimPulse`: SimPulses used to create this pulse
- `particles::edm4hep!MCParticle`: MCParticle that caused the pulse
# Methods
- `setAmplitude(object::edm4eic!SimPulse, v::AbstractVector{Float32})`: assign a set of values to the `amplitude` vector member
- `pushToCalorimeterHits(obj::edm4eic!SimPulse, robj::edm4hep!SimCalorimeterHit)`: push related object to the `calorimeterHits` relation
- `popFromCalorimeterHits(obj::edm4eic!SimPulse)`: pop last related object from `calorimeterHits` relation
- `pushToTrackerHits(obj::edm4eic!SimPulse, robj::edm4hep!SimTrackerHit)`: push related object to the `trackerHits` relation
- `popFromTrackerHits(obj::edm4eic!SimPulse)`: pop last related object from `trackerHits` relation
- `pushToPulses(obj::edm4eic!SimPulse, robj::edm4eic!SimPulse)`: push related object to the `pulses` relation
- `popFromPulses(obj::edm4eic!SimPulse)`: pop last related object from `pulses` relation
- `pushToParticles(obj::edm4eic!SimPulse, robj::edm4hep!MCParticle)`: push related object to the `particles` relation
- `popFromParticles(obj::edm4eic!SimPulse)`: pop last related object from `particles` relation
"""
struct edm4eic!SimPulse <: POD
    index::ObjectID{edm4eic!SimPulse}  # ObjectID of himself
    #---Data Members
    cellID::UInt64                   # ID of the readout cell for this pulse.
    integral::Float32                # Total pulse integral in relevant units.
    position::edm4hep!Vector3f       # Position the pulse is evaluated in world coordinates [mm].
    time::Float32                    # Start time for the pulse in [ns].
    interval::Float32                # Time interval between amplitude values [ns].
    #---VectorMembers
    amplitude::PVector{edm4eic!SimPulse,Float32,1}  # Pulse amplitude in relevant units, sum of amplitude values equals integral
    #---OneToManyRelations
    calorimeterHits::Relation{edm4eic!SimPulse,edm4hep!SimCalorimeterHit,1}  # SimCalorimeterHits used to create this pulse
    trackerHits::Relation{edm4eic!SimPulse,edm4hep!SimTrackerHit,2}  # SimTrackerHits used to create this pulse
    pulses::Relation{edm4eic!SimPulse,edm4eic!SimPulse,3}  # SimPulses used to create this pulse
    particles::Relation{edm4eic!SimPulse,edm4hep!MCParticle,4}  # MCParticle that caused the pulse
end

function edm4eic!SimPulse(;cellID=0, integral=0, position=edm4hep!Vector3f(), time=0, interval=0, amplitude=PVector{edm4eic!SimPulse,Float32,1}(), calorimeterHits=Relation{edm4eic!SimPulse,edm4hep!SimCalorimeterHit,1}(), trackerHits=Relation{edm4eic!SimPulse,edm4hep!SimTrackerHit,2}(), pulses=Relation{edm4eic!SimPulse,edm4eic!SimPulse,3}(), particles=Relation{edm4eic!SimPulse,edm4hep!MCParticle,4}())
    edm4eic!SimPulse(-1, cellID, integral, position, time, interval, amplitude, calorimeterHits, trackerHits, pulses, particles)
end

function pushToCalorimeterHits(c::edm4eic!SimPulse, o::edm4hep!SimCalorimeterHit)
    iszero(c.index) && (c = register(c))
    c = @set c.calorimeterHits = push(c.calorimeterHits, o)
    update(c)
end
function popFromCalorimeterHits(c::edm4eic!SimPulse)
    iszero(c.index) && (c = register(c))
    c = @set c.calorimeterHits = pop(c.calorimeterHits)
    update(c)
end
function pushToTrackerHits(c::edm4eic!SimPulse, o::edm4hep!SimTrackerHit)
    iszero(c.index) && (c = register(c))
    c = @set c.trackerHits = push(c.trackerHits, o)
    update(c)
end
function popFromTrackerHits(c::edm4eic!SimPulse)
    iszero(c.index) && (c = register(c))
    c = @set c.trackerHits = pop(c.trackerHits)
    update(c)
end
function pushToPulses(c::edm4eic!SimPulse, o::edm4eic!SimPulse)
    iszero(c.index) && (c = register(c))
    c = @set c.pulses = push(c.pulses, o)
    update(c)
end
function popFromPulses(c::edm4eic!SimPulse)
    iszero(c.index) && (c = register(c))
    c = @set c.pulses = pop(c.pulses)
    update(c)
end
function pushToParticles(c::edm4eic!SimPulse, o::edm4hep!MCParticle)
    iszero(c.index) && (c = register(c))
    c = @set c.particles = push(c.particles, o)
    update(c)
end
function popFromParticles(c::edm4eic!SimPulse)
    iszero(c.index) && (c = register(c))
    c = @set c.particles = pop(c.particles)
    update(c)
end
function setAmplitude(o::edm4eic!SimPulse, v::AbstractVector{Float32})
    iszero(o.index) && (o = register(o))
    o = @set o.amplitude = v
    update(o)
end
"""
Tracker hit (reconstructed from Raw)
- Author: W. Armstrong, S. Joosten
# Fields
- `cellID::UInt64`: The detector specific (geometrical) cell id.
- `position::edm4hep!Vector3f`: Hit (cell) position [mm]
- `positionError::edm4eic!CovDiag3f`: Covariance Matrix
- `time::Float32`: Hit time [ns]
- `timeError::Float32`: Error on the time
- `edep::Float32`: Energy deposit in this hit [GeV]
- `edepError::Float32`: Error on the energy deposit [GeV]
# Relations
- `rawHit::edm4eic!RawTrackerHit`: Related raw tracker hit
"""
struct edm4eic!TrackerHit <: POD
    index::ObjectID{edm4eic!TrackerHit}  # ObjectID of himself
    #---Data Members
    cellID::UInt64                   # The detector specific (geometrical) cell id.
    position::edm4hep!Vector3f       # Hit (cell) position [mm]
    positionError::edm4eic!CovDiag3f # Covariance Matrix
    time::Float32                    # Hit time [ns]
    timeError::Float32               # Error on the time
    edep::Float32                    # Energy deposit in this hit [GeV]
    edepError::Float32               # Error on the energy deposit [GeV]
    #---OneToOneRelations
    rawHit_idx::ObjectID{edm4eic!RawTrackerHit}  # Related raw tracker hit
end

function edm4eic!TrackerHit(;cellID=0, position=edm4hep!Vector3f(), positionError=edm4eic!CovDiag3f(), time=0, timeError=0, edep=0, edepError=0, rawHit=-1)
    edm4eic!TrackerHit(-1, cellID, position, positionError, time, timeError, edep, edepError, rawHit)
end

function Base.getproperty(obj::edm4eic!TrackerHit, sym::Symbol)
    if sym == :rawHit
        idx = getfield(obj, :rawHit_idx)
        return iszero(idx) ? nothing : convert(edm4eic!RawTrackerHit, idx)
    else # fallback to getfield
        return getfield(obj, sym)
    end
end
"""
2D measurement (on an arbitrary surface)
- Author: W. Deconinck
# Fields
- `surface::UInt64`: Surface for bound coordinates (geometryID)
- `loc::edm4hep!Vector2f`: 2D location on surface
- `time::Float32`: Measurement time
- `covariance::edm4eic!Cov3f`: Covariance on location and time
- `weights::PVector{Float32}`: Weight for each of the hits, mirrors hits array
# Relations
- `hits::edm4eic!TrackerHit`: Hits in this measurement (single or clustered)
# Methods
- `setWeights(object::edm4eic!Measurement2D, v::AbstractVector{Float32})`: assign a set of values to the `weights` vector member
- `pushToHits(obj::edm4eic!Measurement2D, robj::edm4eic!TrackerHit)`: push related object to the `hits` relation
- `popFromHits(obj::edm4eic!Measurement2D)`: pop last related object from `hits` relation
"""
struct edm4eic!Measurement2D <: POD
    index::ObjectID{edm4eic!Measurement2D}  # ObjectID of himself
    #---Data Members
    surface::UInt64                  # Surface for bound coordinates (geometryID)
    loc::edm4hep!Vector2f            # 2D location on surface
    time::Float32                    # Measurement time
    covariance::edm4eic!Cov3f        # Covariance on location and time
    #---VectorMembers
    weights::PVector{edm4eic!Measurement2D,Float32,1}  # Weight for each of the hits, mirrors hits array
    #---OneToManyRelations
    hits::Relation{edm4eic!Measurement2D,edm4eic!TrackerHit,1}  # Hits in this measurement (single or clustered)
end

function edm4eic!Measurement2D(;surface=0, loc=edm4hep!Vector2f(), time=0, covariance=edm4eic!Cov3f(), weights=PVector{edm4eic!Measurement2D,Float32,1}(), hits=Relation{edm4eic!Measurement2D,edm4eic!TrackerHit,1}())
    edm4eic!Measurement2D(-1, surface, loc, time, covariance, weights, hits)
end

function pushToHits(c::edm4eic!Measurement2D, o::edm4eic!TrackerHit)
    iszero(c.index) && (c = register(c))
    c = @set c.hits = push(c.hits, o)
    update(c)
end
function popFromHits(c::edm4eic!Measurement2D)
    iszero(c.index) && (c = register(c))
    c = @set c.hits = pop(c.hits)
    update(c)
end
function setWeights(o::edm4eic!Measurement2D, v::AbstractVector{Float32})
    iszero(o.index) && (o = register(o))
    o = @set o.weights = v
    update(o)
end
"""
Seed info from the realistic seed finder
- Author: S. Li, B. Schmookler, J. Osborn
# Fields
- `perigee::edm4hep!Vector3f`: Vector for the perigee (line surface)
# Relations
- `params::edm4eic!TrackParameters`: Initial track parameters
- `hits::edm4eic!TrackerHit`: Tracker hits triplet for seeding
# Methods
- `pushToHits(obj::edm4eic!TrackSeed, robj::edm4eic!TrackerHit)`: push related object to the `hits` relation
- `popFromHits(obj::edm4eic!TrackSeed)`: pop last related object from `hits` relation
"""
struct edm4eic!TrackSeed <: POD
    index::ObjectID{edm4eic!TrackSeed}  # ObjectID of himself
    #---Data Members
    perigee::edm4hep!Vector3f        # Vector for the perigee (line surface)
    #---OneToManyRelations
    hits::Relation{edm4eic!TrackSeed,edm4eic!TrackerHit,1}  # Tracker hits triplet for seeding
    #---OneToOneRelations
    params_idx::ObjectID{edm4eic!TrackParameters}  # Initial track parameters
end

function edm4eic!TrackSeed(;perigee=edm4hep!Vector3f(), hits=Relation{edm4eic!TrackSeed,edm4eic!TrackerHit,1}(), params=-1)
    edm4eic!TrackSeed(-1, perigee, hits, params)
end

function Base.getproperty(obj::edm4eic!TrackSeed, sym::Symbol)
    if sym == :params
        idx = getfield(obj, :params_idx)
        return iszero(idx) ? nothing : convert(edm4eic!TrackParameters, idx)
    else # fallback to getfield
        return getfield(obj, sym)
    end
end
function pushToHits(c::edm4eic!TrackSeed, o::edm4eic!TrackerHit)
    iszero(c.index) && (c = register(c))
    c = @set c.hits = push(c.hits, o)
    update(c)
end
function popFromHits(c::edm4eic!TrackSeed)
    iszero(c.index) && (c = register(c))
    c = @set c.hits = pop(c.hits)
    update(c)
end
"""
Raw trajectory from the tracking algorithm. What is called hit here is 2d measurement indeed.
- Author: S. Joosten, S. Li
# Fields
- `type::UInt32`: 0 (does not have good track fit), 1 (has good track fit)
- `nStates::UInt32`: Number of tracking steps
- `nMeasurements::UInt32`: Number of hits used
- `nOutliers::UInt32`: Number of hits not considered
- `nHoles::UInt32`: Number of missing hits
- `nSharedHits::UInt32`: Number of shared hits with other trajectories
- `measurementChi2::PVector{Float32}`: Chi2 for each of the measurements
- `outlierChi2::PVector{Float32}`: Chi2 for each of the outliers
# Relations
- `seed::edm4eic!TrackSeed`: Corresponding track seed
- `trackParameters::edm4eic!TrackParameters`: Associated track parameters, if any
- `measurements_deprecated::edm4eic!Measurement2D`: Measurements that were used for this track. Will move this to the edm4eic::Track
- `outliers_deprecated::edm4eic!Measurement2D`: Measurements that were not used for this track. Will move this to the edm4eic::Track
# Methods
- `setMeasurementChi2(object::edm4eic!Trajectory, v::AbstractVector{Float32})`: assign a set of values to the `measurementChi2` vector member
- `setOutlierChi2(object::edm4eic!Trajectory, v::AbstractVector{Float32})`: assign a set of values to the `outlierChi2` vector member
- `pushToTrackParameters(obj::edm4eic!Trajectory, robj::edm4eic!TrackParameters)`: push related object to the `trackParameters` relation
- `popFromTrackParameters(obj::edm4eic!Trajectory)`: pop last related object from `trackParameters` relation
- `pushToMeasurements_deprecated(obj::edm4eic!Trajectory, robj::edm4eic!Measurement2D)`: push related object to the `measurements_deprecated` relation
- `popFromMeasurements_deprecated(obj::edm4eic!Trajectory)`: pop last related object from `measurements_deprecated` relation
- `pushToOutliers_deprecated(obj::edm4eic!Trajectory, robj::edm4eic!Measurement2D)`: push related object to the `outliers_deprecated` relation
- `popFromOutliers_deprecated(obj::edm4eic!Trajectory)`: pop last related object from `outliers_deprecated` relation
"""
struct edm4eic!Trajectory <: POD
    index::ObjectID{edm4eic!Trajectory}  # ObjectID of himself
    #---Data Members
    type::UInt32                     # 0 (does not have good track fit), 1 (has good track fit)
    nStates::UInt32                  # Number of tracking steps
    nMeasurements::UInt32            # Number of hits used
    nOutliers::UInt32                # Number of hits not considered
    nHoles::UInt32                   # Number of missing hits
    nSharedHits::UInt32              # Number of shared hits with other trajectories
    #---VectorMembers
    measurementChi2::PVector{edm4eic!Trajectory,Float32,1}  # Chi2 for each of the measurements
    outlierChi2::PVector{edm4eic!Trajectory,Float32,2}  # Chi2 for each of the outliers
    #---OneToManyRelations
    trackParameters::Relation{edm4eic!Trajectory,edm4eic!TrackParameters,1}  # Associated track parameters, if any
    measurements_deprecated::Relation{edm4eic!Trajectory,edm4eic!Measurement2D,2}  # Measurements that were used for this track. Will move this to the edm4eic::Track
    outliers_deprecated::Relation{edm4eic!Trajectory,edm4eic!Measurement2D,3}  # Measurements that were not used for this track. Will move this to the edm4eic::Track
    #---OneToOneRelations
    seed_idx::ObjectID{edm4eic!TrackSeed}  # Corresponding track seed
end

function edm4eic!Trajectory(;type=0, nStates=0, nMeasurements=0, nOutliers=0, nHoles=0, nSharedHits=0, measurementChi2=PVector{edm4eic!Trajectory,Float32,1}(), outlierChi2=PVector{edm4eic!Trajectory,Float32,2}(), trackParameters=Relation{edm4eic!Trajectory,edm4eic!TrackParameters,1}(), measurements_deprecated=Relation{edm4eic!Trajectory,edm4eic!Measurement2D,2}(), outliers_deprecated=Relation{edm4eic!Trajectory,edm4eic!Measurement2D,3}(), seed=-1)
    edm4eic!Trajectory(-1, type, nStates, nMeasurements, nOutliers, nHoles, nSharedHits, measurementChi2, outlierChi2, trackParameters, measurements_deprecated, outliers_deprecated, seed)
end

function Base.getproperty(obj::edm4eic!Trajectory, sym::Symbol)
    if sym == :seed
        idx = getfield(obj, :seed_idx)
        return iszero(idx) ? nothing : convert(edm4eic!TrackSeed, idx)
    else # fallback to getfield
        return getfield(obj, sym)
    end
end
function pushToTrackParameters(c::edm4eic!Trajectory, o::edm4eic!TrackParameters)
    iszero(c.index) && (c = register(c))
    c = @set c.trackParameters = push(c.trackParameters, o)
    update(c)
end
function popFromTrackParameters(c::edm4eic!Trajectory)
    iszero(c.index) && (c = register(c))
    c = @set c.trackParameters = pop(c.trackParameters)
    update(c)
end
function pushToMeasurements_deprecated(c::edm4eic!Trajectory, o::edm4eic!Measurement2D)
    iszero(c.index) && (c = register(c))
    c = @set c.measurements_deprecated = push(c.measurements_deprecated, o)
    update(c)
end
function popFromMeasurements_deprecated(c::edm4eic!Trajectory)
    iszero(c.index) && (c = register(c))
    c = @set c.measurements_deprecated = pop(c.measurements_deprecated)
    update(c)
end
function pushToOutliers_deprecated(c::edm4eic!Trajectory, o::edm4eic!Measurement2D)
    iszero(c.index) && (c = register(c))
    c = @set c.outliers_deprecated = push(c.outliers_deprecated, o)
    update(c)
end
function popFromOutliers_deprecated(c::edm4eic!Trajectory)
    iszero(c.index) && (c = register(c))
    c = @set c.outliers_deprecated = pop(c.outliers_deprecated)
    update(c)
end
function setMeasurementChi2(o::edm4eic!Trajectory, v::AbstractVector{Float32})
    iszero(o.index) && (o = register(o))
    o = @set o.measurementChi2 = v
    update(o)
end
function setOutlierChi2(o::edm4eic!Trajectory, v::AbstractVector{Float32})
    iszero(o.index) && (o = register(o))
    o = @set o.outlierChi2 = v
    update(o)
end
"""
Track information at the vertex
- Author: S. Joosten, J. Osborn
# Fields
- `type::Int32`: Flag that defines the type of track
- `position::edm4hep!Vector3f`: Track 3-position at the vertex
- `momentum::edm4hep!Vector3f`: Track 3-momentum at the vertex [GeV]
- `positionMomentumCovariance::edm4eic!Cov6f`: Covariance matrix in basis [x,y,z,px,py,pz]
- `time::Float32`: Track time at the vertex [ns]
- `timeError::Float32`: Error on the track vertex time
- `charge::Float32`: Particle charge
- `chi2::Float32`: Total chi2
- `ndf::UInt32`: Number of degrees of freedom
- `pdg::Int32`: PDG particle ID hypothesis
# Relations
- `trajectory::edm4eic!Trajectory`: Trajectory of this track
- `measurements::edm4eic!Measurement2D`: Measurements that were used for this track
- `tracks::edm4eic!Track`: Tracks (segments) that have been combined to create this track
# Methods
- `pushToMeasurements(obj::edm4eic!Track, robj::edm4eic!Measurement2D)`: push related object to the `measurements` relation
- `popFromMeasurements(obj::edm4eic!Track)`: pop last related object from `measurements` relation
- `pushToTracks(obj::edm4eic!Track, robj::edm4eic!Track)`: push related object to the `tracks` relation
- `popFromTracks(obj::edm4eic!Track)`: pop last related object from `tracks` relation
"""
struct edm4eic!Track <: POD
    index::ObjectID{edm4eic!Track}   # ObjectID of himself
    #---Data Members
    type::Int32                      # Flag that defines the type of track
    position::edm4hep!Vector3f       # Track 3-position at the vertex
    momentum::edm4hep!Vector3f       # Track 3-momentum at the vertex [GeV]
    positionMomentumCovariance::edm4eic!Cov6f  # Covariance matrix in basis [x,y,z,px,py,pz]
    time::Float32                    # Track time at the vertex [ns]
    timeError::Float32               # Error on the track vertex time
    charge::Float32                  # Particle charge
    chi2::Float32                    # Total chi2
    ndf::UInt32                      # Number of degrees of freedom
    pdg::Int32                       # PDG particle ID hypothesis
    #---OneToManyRelations
    measurements::Relation{edm4eic!Track,edm4eic!Measurement2D,1}  # Measurements that were used for this track
    tracks::Relation{edm4eic!Track,edm4eic!Track,2}  # Tracks (segments) that have been combined to create this track
    #---OneToOneRelations
    trajectory_idx::ObjectID{edm4eic!Trajectory}  # Trajectory of this track
end

function edm4eic!Track(;type=0, position=edm4hep!Vector3f(), momentum=edm4hep!Vector3f(), positionMomentumCovariance=edm4eic!Cov6f(), time=0, timeError=0, charge=0, chi2=0, ndf=0, pdg=0, measurements=Relation{edm4eic!Track,edm4eic!Measurement2D,1}(), tracks=Relation{edm4eic!Track,edm4eic!Track,2}(), trajectory=-1)
    edm4eic!Track(-1, type, position, momentum, positionMomentumCovariance, time, timeError, charge, chi2, ndf, pdg, measurements, tracks, trajectory)
end

function Base.getproperty(obj::edm4eic!Track, sym::Symbol)
    if sym == :trajectory
        idx = getfield(obj, :trajectory_idx)
        return iszero(idx) ? nothing : convert(edm4eic!Trajectory, idx)
    else # fallback to getfield
        return getfield(obj, sym)
    end
end
function pushToMeasurements(c::edm4eic!Track, o::edm4eic!Measurement2D)
    iszero(c.index) && (c = register(c))
    c = @set c.measurements = push(c.measurements, o)
    update(c)
end
function popFromMeasurements(c::edm4eic!Track)
    iszero(c.index) && (c = register(c))
    c = @set c.measurements = pop(c.measurements)
    update(c)
end
function pushToTracks(c::edm4eic!Track, o::edm4eic!Track)
    iszero(c.index) && (c = register(c))
    c = @set c.tracks = push(c.tracks, o)
    update(c)
end
function popFromTracks(c::edm4eic!Track)
    iszero(c.index) && (c = register(c))
    c = @set c.tracks = pop(c.tracks)
    update(c)
end
"""
EIC Reconstructed Particle
- Author: W. Armstrong, S. Joosten, F. Gaede
# Fields
- `type::Int32`: type of reconstructed particle. Check/set collection parameters ReconstructedParticleTypeNames and ReconstructedParticleTypeValues.
- `energy::Float32`: [GeV] energy of the reconstructed particle. Four momentum state is not kept consistent internally.
- `momentum::edm4hep!Vector3f`: [GeV] particle momentum. Four momentum state is not kept consistent internally.
- `referencePoint::edm4hep!Vector3f`: [mm] reference, i.e. where the particle has been measured
- `charge::Float32`: charge of the reconstructed particle.
- `mass::Float32`: [GeV] mass of the reconstructed particle, set independently from four vector. Four momentum state is not kept consistent internally.
- `goodnessOfPID::Float32`: overall goodness of the PID on a scale of [0;1]
- `covMatrix::edm4eic!Cov4f`: covariance matrix of the reconstructed particle 4vector (10 parameters).
- `PDG::Int32`: PDG code for this particle
# Relations
- `startVertex::edm4eic!Vertex`: Start vertex associated to this particle
- `particleIDUsed::edm4hep!ParticleID`: particle ID used for the kinematics of this particle
- `clusters::edm4eic!Cluster`: Clusters used for this particle
- `tracks::edm4eic!Track`: Tracks used for this particle
- `particles::edm4eic!ReconstructedParticle`: Reconstructed particles that have been combined to this particle
- `particleIDs::edm4hep!ParticleID`: All associated particle IDs for this particle (not sorted by likelihood)
# Methods
- `pushToClusters(obj::edm4eic!ReconstructedParticle, robj::edm4eic!Cluster)`: push related object to the `clusters` relation
- `popFromClusters(obj::edm4eic!ReconstructedParticle)`: pop last related object from `clusters` relation
- `pushToTracks(obj::edm4eic!ReconstructedParticle, robj::edm4eic!Track)`: push related object to the `tracks` relation
- `popFromTracks(obj::edm4eic!ReconstructedParticle)`: pop last related object from `tracks` relation
- `pushToParticles(obj::edm4eic!ReconstructedParticle, robj::edm4eic!ReconstructedParticle)`: push related object to the `particles` relation
- `popFromParticles(obj::edm4eic!ReconstructedParticle)`: pop last related object from `particles` relation
- `pushToParticleIDs(obj::edm4eic!ReconstructedParticle, robj::edm4hep!ParticleID)`: push related object to the `particleIDs` relation
- `popFromParticleIDs(obj::edm4eic!ReconstructedParticle)`: pop last related object from `particleIDs` relation
"""
struct edm4eic!ReconstructedParticle <: POD
    index::ObjectID{edm4eic!ReconstructedParticle}  # ObjectID of himself
    #---Data Members
    type::Int32                      # type of reconstructed particle. Check/set collection parameters ReconstructedParticleTypeNames and ReconstructedParticleTypeValues.
    energy::Float32                  # [GeV] energy of the reconstructed particle. Four momentum state is not kept consistent internally.
    momentum::edm4hep!Vector3f       # [GeV] particle momentum. Four momentum state is not kept consistent internally.
    referencePoint::edm4hep!Vector3f # [mm] reference, i.e. where the particle has been measured
    charge::Float32                  # charge of the reconstructed particle.
    mass::Float32                    # [GeV] mass of the reconstructed particle, set independently from four vector. Four momentum state is not kept consistent internally.
    goodnessOfPID::Float32           # overall goodness of the PID on a scale of [0;1]
    covMatrix::edm4eic!Cov4f         # covariance matrix of the reconstructed particle 4vector (10 parameters).
    PDG::Int32                       # PDG code for this particle
    #---OneToManyRelations
    clusters::Relation{edm4eic!ReconstructedParticle,edm4eic!Cluster,1}  # Clusters used for this particle
    tracks::Relation{edm4eic!ReconstructedParticle,edm4eic!Track,2}  # Tracks used for this particle
    particles::Relation{edm4eic!ReconstructedParticle,edm4eic!ReconstructedParticle,3}  # Reconstructed particles that have been combined to this particle
    particleIDs::Relation{edm4eic!ReconstructedParticle,edm4hep!ParticleID,4}  # All associated particle IDs for this particle (not sorted by likelihood)
    #---OneToOneRelations
    startVertex_idx::ObjectID{edm4eic!Vertex}  # Start vertex associated to this particle
    particleIDUsed_idx::ObjectID{edm4hep!ParticleID}  # particle ID used for the kinematics of this particle
end

function edm4eic!ReconstructedParticle(;type=0, energy=0, momentum=edm4hep!Vector3f(), referencePoint=edm4hep!Vector3f(), charge=0, mass=0, goodnessOfPID=0, covMatrix=edm4eic!Cov4f(), PDG=0, clusters=Relation{edm4eic!ReconstructedParticle,edm4eic!Cluster,1}(), tracks=Relation{edm4eic!ReconstructedParticle,edm4eic!Track,2}(), particles=Relation{edm4eic!ReconstructedParticle,edm4eic!ReconstructedParticle,3}(), particleIDs=Relation{edm4eic!ReconstructedParticle,edm4hep!ParticleID,4}(), startVertex=-1, particleIDUsed=-1)
    edm4eic!ReconstructedParticle(-1, type, energy, momentum, referencePoint, charge, mass, goodnessOfPID, covMatrix, PDG, clusters, tracks, particles, particleIDs, startVertex, particleIDUsed)
end

function Base.getproperty(obj::edm4eic!ReconstructedParticle, sym::Symbol)
    if sym == :startVertex
        idx = getfield(obj, :startVertex_idx)
        return iszero(idx) ? nothing : convert(edm4eic!Vertex, idx)
    elseif sym == :particleIDUsed
        idx = getfield(obj, :particleIDUsed_idx)
        return iszero(idx) ? nothing : convert(edm4hep!ParticleID, idx)
    else # fallback to getfield
        return getfield(obj, sym)
    end
end
function pushToClusters(c::edm4eic!ReconstructedParticle, o::edm4eic!Cluster)
    iszero(c.index) && (c = register(c))
    c = @set c.clusters = push(c.clusters, o)
    update(c)
end
function popFromClusters(c::edm4eic!ReconstructedParticle)
    iszero(c.index) && (c = register(c))
    c = @set c.clusters = pop(c.clusters)
    update(c)
end
function pushToTracks(c::edm4eic!ReconstructedParticle, o::edm4eic!Track)
    iszero(c.index) && (c = register(c))
    c = @set c.tracks = push(c.tracks, o)
    update(c)
end
function popFromTracks(c::edm4eic!ReconstructedParticle)
    iszero(c.index) && (c = register(c))
    c = @set c.tracks = pop(c.tracks)
    update(c)
end
function pushToParticles(c::edm4eic!ReconstructedParticle, o::edm4eic!ReconstructedParticle)
    iszero(c.index) && (c = register(c))
    c = @set c.particles = push(c.particles, o)
    update(c)
end
function popFromParticles(c::edm4eic!ReconstructedParticle)
    iszero(c.index) && (c = register(c))
    c = @set c.particles = pop(c.particles)
    update(c)
end
function pushToParticleIDs(c::edm4eic!ReconstructedParticle, o::edm4hep!ParticleID)
    iszero(c.index) && (c = register(c))
    c = @set c.particleIDs = push(c.particleIDs, o)
    update(c)
end
function popFromParticleIDs(c::edm4eic!ReconstructedParticle)
    iszero(c.index) && (c = register(c))
    c = @set c.particleIDs = pop(c.particleIDs)
    update(c)
end
"""
Kinematic variables for DIS events
- Author: S. Joosten, W. Deconinck
# Fields
- `x::Float32`: Bjorken x (Q2/2P.q)
- `Q2::Float32`: Four-momentum transfer squared [GeV^2]
- `W::Float32`: Invariant mass of final state [GeV]
- `y::Float32`: Inelasticity (P.q/P.k)
- `nu::Float32`: Energy transfer P.q/M [GeV]
# Relations
- `scat::edm4eic!ReconstructedParticle`: Associated scattered electron (if identified)
"""
struct edm4eic!InclusiveKinematics <: POD
    index::ObjectID{edm4eic!InclusiveKinematics}  # ObjectID of himself
    #---Data Members
    x::Float32                       # Bjorken x (Q2/2P.q)
    Q2::Float32                      # Four-momentum transfer squared [GeV^2]
    W::Float32                       # Invariant mass of final state [GeV]
    y::Float32                       # Inelasticity (P.q/P.k)
    nu::Float32                      # Energy transfer P.q/M [GeV]
    #---OneToOneRelations
    scat_idx::ObjectID{edm4eic!ReconstructedParticle}  # Associated scattered electron (if identified)
end

function edm4eic!InclusiveKinematics(;x=0, Q2=0, W=0, y=0, nu=0, scat=-1)
    edm4eic!InclusiveKinematics(-1, x, Q2, W, y, nu, scat)
end

function Base.getproperty(obj::edm4eic!InclusiveKinematics, sym::Symbol)
    if sym == :scat
        idx = getfield(obj, :scat_idx)
        return iszero(idx) ? nothing : convert(edm4eic!ReconstructedParticle, idx)
    else # fallback to getfield
        return getfield(obj, sym)
    end
end
"""
Match between a Cluster and a Track
- Author: D. Anderson, D. Brandenburg, D. Kalinkin, S. Joosten
# Fields
- `weight::Float32`: weight of this association
# Relations
- `cluster::edm4eic!Cluster`: reference to the cluster
- `track::edm4eic!Track`: reference to the track
"""
struct edm4eic!TrackClusterMatch <: POD
    index::ObjectID{edm4eic!TrackClusterMatch}  # ObjectID of himself
    #---Data Members
    weight::Float32                  # weight of this association
    #---OneToOneRelations
    cluster_idx::ObjectID{edm4eic!Cluster}  # reference to the cluster
    track_idx::ObjectID{edm4eic!Track}  # reference to the track
end

function edm4eic!TrackClusterMatch(;weight=0, cluster=-1, track=-1)
    edm4eic!TrackClusterMatch(-1, weight, cluster, track)
end

function Base.getproperty(obj::edm4eic!TrackClusterMatch, sym::Symbol)
    if sym == :cluster
        idx = getfield(obj, :cluster_idx)
        return iszero(idx) ? nothing : convert(edm4eic!Cluster, idx)
    elseif sym == :track
        idx = getfield(obj, :track_idx)
        return iszero(idx) ? nothing : convert(edm4eic!Track, idx)
    else # fallback to getfield
        return getfield(obj, sym)
    end
end
"""
Association between a Track and a MCParticle
- Author: S. Joosten
# Fields
- `simID::UInt32`: Index of corresponding MCParticle (position in MCParticles array)
- `recID::UInt32`: Index of corresponding Track (position in Tracks array)
- `weight::Float32`: weight of this association
# Relations
- `rec::edm4eic!Track`: reference to the track
- `sim::edm4hep!MCParticle`: reference to the Monte-Carlo particle
"""
struct edm4eic!MCRecoTrackParticleAssociation <: POD
    index::ObjectID{edm4eic!MCRecoTrackParticleAssociation}  # ObjectID of himself
    #---Data Members
    simID::UInt32                    # Index of corresponding MCParticle (position in MCParticles array)
    recID::UInt32                    # Index of corresponding Track (position in Tracks array)
    weight::Float32                  # weight of this association
    #---OneToOneRelations
    rec_idx::ObjectID{edm4eic!Track} # reference to the track
    sim_idx::ObjectID{edm4hep!MCParticle}  # reference to the Monte-Carlo particle
end

function edm4eic!MCRecoTrackParticleAssociation(;simID=0, recID=0, weight=0, rec=-1, sim=-1)
    edm4eic!MCRecoTrackParticleAssociation(-1, simID, recID, weight, rec, sim)
end

function Base.getproperty(obj::edm4eic!MCRecoTrackParticleAssociation, sym::Symbol)
    if sym == :rec
        idx = getfield(obj, :rec_idx)
        return iszero(idx) ? nothing : convert(edm4eic!Track, idx)
    elseif sym == :sim
        idx = getfield(obj, :sim_idx)
        return iszero(idx) ? nothing : convert(edm4hep!MCParticle, idx)
    else # fallback to getfield
        return getfield(obj, sym)
    end
end
"""
A track segment defined by one or more points along a track.
- Author: S. Joosten
# Fields
- `length::Float32`: Pathlength from the first to the last point
- `lengthError::Float32`: Error on the segment length
- `points::PVector{edm4eic!TrackPoint}`: Points where the track parameters were evaluated
# Relations
- `track::edm4eic!Track`: Track used for this projection
# Methods
- `setPoints(object::edm4eic!TrackSegment, v::AbstractVector{edm4eic!TrackPoint})`: assign a set of values to the `points` vector member
"""
struct edm4eic!TrackSegment <: POD
    index::ObjectID{edm4eic!TrackSegment}  # ObjectID of himself
    #---Data Members
    length::Float32                  # Pathlength from the first to the last point
    lengthError::Float32             # Error on the segment length
    #---VectorMembers
    points::PVector{edm4eic!TrackSegment,edm4eic!TrackPoint,1}  # Points where the track parameters were evaluated
    #---OneToOneRelations
    track_idx::ObjectID{edm4eic!Track}  # Track used for this projection
end

function edm4eic!TrackSegment(;length=0, lengthError=0, points=PVector{edm4eic!TrackSegment,edm4eic!TrackPoint,1}(), track=-1)
    edm4eic!TrackSegment(-1, length, lengthError, points, track)
end

function Base.getproperty(obj::edm4eic!TrackSegment, sym::Symbol)
    if sym == :track
        idx = getfield(obj, :track_idx)
        return iszero(idx) ? nothing : convert(edm4eic!Track, idx)
    else # fallback to getfield
        return getfield(obj, sym)
    end
end
function setPoints(o::edm4eic!TrackSegment, v::AbstractVector{edm4eic!TrackPoint})
    iszero(o.index) && (o = register(o))
    o = @set o.points = v
    update(o)
end
"""
Cherenkov detector PID
- Author: A. Kiselev, C. Chatterjee, C. Dilks
# Fields
- `npe::Float32`: Overall photoelectron count
- `refractiveIndex::Float32`: Average refractive index at the Cherenkov photons' vertices
- `photonEnergy::Float32`: Average energy for these Cherenkov photons [GeV]
- `hypotheses::PVector{edm4eic!CherenkovParticleIDHypothesis}`: Evaluated PDG hypotheses
- `thetaPhiPhotons::PVector{edm4hep!Vector2f}`: estimated (theta,phi) for each Cherenkov photon
# Relations
- `chargedParticle::edm4eic!TrackSegment`: reconstructed charged particle
- `rawHitAssociations::edm4eic!MCRecoTrackerHitAssociation`: raw sensor hits, associated with MC hits
# Methods
- `setHypotheses(object::edm4eic!CherenkovParticleID, v::AbstractVector{edm4eic!CherenkovParticleIDHypothesis})`: assign a set of values to the `hypotheses` vector member
- `setThetaPhiPhotons(object::edm4eic!CherenkovParticleID, v::AbstractVector{edm4hep!Vector2f})`: assign a set of values to the `thetaPhiPhotons` vector member
- `pushToRawHitAssociations(obj::edm4eic!CherenkovParticleID, robj::edm4eic!MCRecoTrackerHitAssociation)`: push related object to the `rawHitAssociations` relation
- `popFromRawHitAssociations(obj::edm4eic!CherenkovParticleID)`: pop last related object from `rawHitAssociations` relation
"""
struct edm4eic!CherenkovParticleID <: POD
    index::ObjectID{edm4eic!CherenkovParticleID}  # ObjectID of himself
    #---Data Members
    npe::Float32                     # Overall photoelectron count
    refractiveIndex::Float32         # Average refractive index at the Cherenkov photons' vertices
    photonEnergy::Float32            # Average energy for these Cherenkov photons [GeV]
    #---VectorMembers
    hypotheses::PVector{edm4eic!CherenkovParticleID,edm4eic!CherenkovParticleIDHypothesis,1}  # Evaluated PDG hypotheses
    thetaPhiPhotons::PVector{edm4eic!CherenkovParticleID,edm4hep!Vector2f,2}  # estimated (theta,phi) for each Cherenkov photon
    #---OneToManyRelations
    rawHitAssociations::Relation{edm4eic!CherenkovParticleID,edm4eic!MCRecoTrackerHitAssociation,1}  # raw sensor hits, associated with MC hits
    #---OneToOneRelations
    chargedParticle_idx::ObjectID{edm4eic!TrackSegment}  # reconstructed charged particle
end

function edm4eic!CherenkovParticleID(;npe=0, refractiveIndex=0, photonEnergy=0, hypotheses=PVector{edm4eic!CherenkovParticleID,edm4eic!CherenkovParticleIDHypothesis,1}(), thetaPhiPhotons=PVector{edm4eic!CherenkovParticleID,edm4hep!Vector2f,2}(), rawHitAssociations=Relation{edm4eic!CherenkovParticleID,edm4eic!MCRecoTrackerHitAssociation,1}(), chargedParticle=-1)
    edm4eic!CherenkovParticleID(-1, npe, refractiveIndex, photonEnergy, hypotheses, thetaPhiPhotons, rawHitAssociations, chargedParticle)
end

function Base.getproperty(obj::edm4eic!CherenkovParticleID, sym::Symbol)
    if sym == :chargedParticle
        idx = getfield(obj, :chargedParticle_idx)
        return iszero(idx) ? nothing : convert(edm4eic!TrackSegment, idx)
    else # fallback to getfield
        return getfield(obj, sym)
    end
end
function pushToRawHitAssociations(c::edm4eic!CherenkovParticleID, o::edm4eic!MCRecoTrackerHitAssociation)
    iszero(c.index) && (c = register(c))
    c = @set c.rawHitAssociations = push(c.rawHitAssociations, o)
    update(c)
end
function popFromRawHitAssociations(c::edm4eic!CherenkovParticleID)
    iszero(c.index) && (c = register(c))
    c = @set c.rawHitAssociations = pop(c.rawHitAssociations)
    update(c)
end
function setHypotheses(o::edm4eic!CherenkovParticleID, v::AbstractVector{edm4eic!CherenkovParticleIDHypothesis})
    iszero(o.index) && (o = register(o))
    o = @set o.hypotheses = v
    update(o)
end
function setThetaPhiPhotons(o::edm4eic!CherenkovParticleID, v::AbstractVector{edm4hep!Vector2f})
    iszero(o.index) && (o = register(o))
    o = @set o.thetaPhiPhotons = v
    update(o)
end
"""
Association between a CaloHit and the corresponding simulated CaloHit
- Author: C. Bernet, B. Hegner
# Fields
- `weight::Float32`: weight of this association
# Relations
- `rec::edm4hep!CalorimeterHit`: reference to the reconstructed hit
- `sim::edm4hep!SimCalorimeterHit`: reference to the simulated hit
"""
struct edm4hep!MCRecoCaloAssociation <: POD
    index::ObjectID{edm4hep!MCRecoCaloAssociation}  # ObjectID of himself
    #---Data Members
    weight::Float32                  # weight of this association
    #---OneToOneRelations
    rec_idx::ObjectID{edm4hep!CalorimeterHit}  # reference to the reconstructed hit
    sim_idx::ObjectID{edm4hep!SimCalorimeterHit}  # reference to the simulated hit
end

function edm4hep!MCRecoCaloAssociation(;weight=0, rec=-1, sim=-1)
    edm4hep!MCRecoCaloAssociation(-1, weight, rec, sim)
end

function Base.getproperty(obj::edm4hep!MCRecoCaloAssociation, sym::Symbol)
    if sym == :rec
        idx = getfield(obj, :rec_idx)
        return iszero(idx) ? nothing : convert(edm4hep!CalorimeterHit, idx)
    elseif sym == :sim
        idx = getfield(obj, :sim_idx)
        return iszero(idx) ? nothing : convert(edm4hep!SimCalorimeterHit, idx)
    else # fallback to getfield
        return getfield(obj, sym)
    end
end
"""
Summed quantities of the hadronic final state
- Author: T. Kutz
# Fields
- `sigma::Float32`: Longitudinal energy-momentum balance (aka E - pz)
- `pT::Float32`: Transverse momentum
- `gamma::Float32`: Hadronic angle
# Relations
- `hadrons::edm4eic!ReconstructedParticle`: Reconstructed hadrons used in calculation
# Methods
- `pushToHadrons(obj::edm4eic!HadronicFinalState, robj::edm4eic!ReconstructedParticle)`: push related object to the `hadrons` relation
- `popFromHadrons(obj::edm4eic!HadronicFinalState)`: pop last related object from `hadrons` relation
"""
struct edm4eic!HadronicFinalState <: POD
    index::ObjectID{edm4eic!HadronicFinalState}  # ObjectID of himself
    #---Data Members
    sigma::Float32                   # Longitudinal energy-momentum balance (aka E - pz)
    pT::Float32                      # Transverse momentum
    gamma::Float32                   # Hadronic angle
    #---OneToManyRelations
    hadrons::Relation{edm4eic!HadronicFinalState,edm4eic!ReconstructedParticle,1}  # Reconstructed hadrons used in calculation
end

function edm4eic!HadronicFinalState(;sigma=0, pT=0, gamma=0, hadrons=Relation{edm4eic!HadronicFinalState,edm4eic!ReconstructedParticle,1}())
    edm4eic!HadronicFinalState(-1, sigma, pT, gamma, hadrons)
end

function pushToHadrons(c::edm4eic!HadronicFinalState, o::edm4eic!ReconstructedParticle)
    iszero(c.index) && (c = register(c))
    c = @set c.hadrons = push(c.hadrons, o)
    update(c)
end
function popFromHadrons(c::edm4eic!HadronicFinalState)
    iszero(c.index) && (c = register(c))
    c = @set c.hadrons = pop(c.hadrons)
    update(c)
end
"""
Event Header. Additional parameters are assumed to go into the metadata tree.
- Author: F.Gaede
# Fields
- `eventNumber::Int32`: event number
- `runNumber::Int32`: run number
- `timeStamp::UInt64`: time stamp
- `weight::Float32`: event weight
"""
struct edm4hep!EventHeader <: POD
    index::ObjectID{edm4hep!EventHeader}  # ObjectID of himself
    #---Data Members
    eventNumber::Int32               # event number
    runNumber::Int32                 # run number
    timeStamp::UInt64                # time stamp
    weight::Float32                  # event weight
end

function edm4hep!EventHeader(;eventNumber=0, runNumber=0, timeStamp=0, weight=0)
    edm4hep!EventHeader(-1, eventNumber, runNumber, timeStamp, weight)
end

"""
Vertex
- Author: F.Gaede, DESY
# Fields
- `primary::Int32`: boolean flag, if vertex is the primary vertex of the event
- `chi2::Float32`: chi-squared of the vertex fit
- `probability::Float32`: probability of the vertex fit
- `position::edm4hep!Vector3f`: [mm] position of the vertex.
- `covMatrix::SVector{6,Float32}`: covariance matrix of the position (stored as lower triangle matrix, i.e. cov(xx),cov(y,x),cov(z,x),cov(y,y),... )
- `algorithmType::Int32`: type code for the algorithm that has been used to create the vertex - check/set the collection parameters AlgorithmName and AlgorithmType.
- `parameters::PVector{Float32}`: additional parameters related to this vertex - check/set the collection parameter "VertexParameterNames" for the parameters meaning.
# Relations
- `associatedParticle::edm4hep!POD`: reconstructed particle associated to this vertex.
# Methods
- `setParameters(object::edm4hep!Vertex, v::AbstractVector{Float32})`: assign a set of values to the `parameters` vector member
"""
struct edm4hep!Vertex <: POD
    index::ObjectID{edm4hep!Vertex}  # ObjectID of himself
    #---Data Members
    primary::Int32                   # boolean flag, if vertex is the primary vertex of the event
    chi2::Float32                    # chi-squared of the vertex fit
    probability::Float32             # probability of the vertex fit
    position::edm4hep!Vector3f       # [mm] position of the vertex.
    covMatrix::SVector{6,Float32}    # covariance matrix of the position (stored as lower triangle matrix, i.e. cov(xx),cov(y,x),cov(z,x),cov(y,y),... )
    algorithmType::Int32             # type code for the algorithm that has been used to create the vertex - check/set the collection parameters AlgorithmName and AlgorithmType.
    #---VectorMembers
    parameters::PVector{edm4hep!Vertex,Float32,1}  # additional parameters related to this vertex - check/set the collection parameter "VertexParameterNames" for the parameters meaning.
    #---OneToOneRelations
    associatedParticle_idx::ObjectID{edm4hep!POD}  # reconstructed particle associated to this vertex.
end

function edm4hep!Vertex(;primary=0, chi2=0, probability=0, position=edm4hep!Vector3f(), covMatrix=zero(SVector{6,Float32}), algorithmType=0, parameters=PVector{edm4hep!Vertex,Float32,1}(), associatedParticle=-1)
    edm4hep!Vertex(-1, primary, chi2, probability, position, covMatrix, algorithmType, parameters, associatedParticle)
end

function Base.getproperty(obj::edm4hep!Vertex, sym::Symbol)
    if sym == :associatedParticle
        idx = getfield(obj, :associatedParticle_idx)
        return iszero(idx) ? nothing : convert(edm4hep!POD, idx)
    else # fallback to getfield
        return getfield(obj, sym)
    end
end
function setParameters(o::edm4hep!Vertex, v::AbstractVector{Float32})
    iszero(o.index) && (o = register(o))
    o = @set o.parameters = v
    update(o)
end
"""
Reconstructed track
- Author: F.Gaede, DESY
# Fields
- `type::Int32`: flagword that defines the type of track.Bits 16-31 are used internally
- `chi2::Float32`: Chi^2 of the track fit
- `ndf::Int32`: number of degrees of freedom of the track fit
- `dEdx::Float32`: dEdx of the track.
- `dEdxError::Float32`: error of dEdx.
- `radiusOfInnermostHit::Float32`: radius of the innermost hit that has been used in the track fit
- `subdetectorHitNumbers::PVector{Int32}`: number of hits in particular subdetectors.Check/set collection variable TrackSubdetectorNames for decoding the indices
- `trackStates::PVector{edm4hep!TrackState}`: track states
- `dxQuantities::PVector{edm4hep!Quantity}`: different measurements of dx quantities
# Relations
- `trackerHits::edm4hep!TrackerHit`: hits that have been used to create this track
- `tracks::edm4hep!Track`: tracks (segments) that have been combined to create this track
# Methods
- `setSubdetectorHitNumbers(object::edm4hep!Track, v::AbstractVector{Int32})`: assign a set of values to the `subdetectorHitNumbers` vector member
- `setTrackStates(object::edm4hep!Track, v::AbstractVector{edm4hep!TrackState})`: assign a set of values to the `trackStates` vector member
- `setDxQuantities(object::edm4hep!Track, v::AbstractVector{edm4hep!Quantity})`: assign a set of values to the `dxQuantities` vector member
- `pushToTrackerHits(obj::edm4hep!Track, robj::edm4hep!TrackerHit)`: push related object to the `trackerHits` relation
- `popFromTrackerHits(obj::edm4hep!Track)`: pop last related object from `trackerHits` relation
- `pushToTracks(obj::edm4hep!Track, robj::edm4hep!Track)`: push related object to the `tracks` relation
- `popFromTracks(obj::edm4hep!Track)`: pop last related object from `tracks` relation
"""
struct edm4hep!Track <: POD
    index::ObjectID{edm4hep!Track}   # ObjectID of himself
    #---Data Members
    type::Int32                      # flagword that defines the type of track.Bits 16-31 are used internally
    chi2::Float32                    # Chi^2 of the track fit
    ndf::Int32                       # number of degrees of freedom of the track fit
    dEdx::Float32                    # dEdx of the track.
    dEdxError::Float32               # error of dEdx.
    radiusOfInnermostHit::Float32    # radius of the innermost hit that has been used in the track fit
    #---VectorMembers
    subdetectorHitNumbers::PVector{edm4hep!Track,Int32,1}  # number of hits in particular subdetectors.Check/set collection variable TrackSubdetectorNames for decoding the indices
    trackStates::PVector{edm4hep!Track,edm4hep!TrackState,2}  # track states
    dxQuantities::PVector{edm4hep!Track,edm4hep!Quantity,3}  # different measurements of dx quantities
    #---OneToManyRelations
    trackerHits::Relation{edm4hep!Track,edm4hep!TrackerHit,1}  # hits that have been used to create this track
    tracks::Relation{edm4hep!Track,edm4hep!Track,2}  # tracks (segments) that have been combined to create this track
end

function edm4hep!Track(;type=0, chi2=0, ndf=0, dEdx=0, dEdxError=0, radiusOfInnermostHit=0, subdetectorHitNumbers=PVector{edm4hep!Track,Int32,1}(), trackStates=PVector{edm4hep!Track,edm4hep!TrackState,2}(), dxQuantities=PVector{edm4hep!Track,edm4hep!Quantity,3}(), trackerHits=Relation{edm4hep!Track,edm4hep!TrackerHit,1}(), tracks=Relation{edm4hep!Track,edm4hep!Track,2}())
    edm4hep!Track(-1, type, chi2, ndf, dEdx, dEdxError, radiusOfInnermostHit, subdetectorHitNumbers, trackStates, dxQuantities, trackerHits, tracks)
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
function setTrackStates(o::edm4hep!Track, v::AbstractVector{edm4hep!TrackState})
    iszero(o.index) && (o = register(o))
    o = @set o.trackStates = v
    update(o)
end
function setDxQuantities(o::edm4hep!Track, v::AbstractVector{edm4hep!Quantity})
    iszero(o.index) && (o = register(o))
    o = @set o.dxQuantities = v
    update(o)
end
"""
Association between a Track and a MCParticle
- Author: Placido Fernandez Declara
# Fields
- `weight::Float32`: weight of this association
# Relations
- `rec::edm4hep!Track`: reference to the track
- `sim::edm4hep!MCParticle`: reference to the Monte-Carlo particle
"""
struct edm4hep!MCRecoTrackParticleAssociation <: POD
    index::ObjectID{edm4hep!MCRecoTrackParticleAssociation}  # ObjectID of himself
    #---Data Members
    weight::Float32                  # weight of this association
    #---OneToOneRelations
    rec_idx::ObjectID{edm4hep!Track} # reference to the track
    sim_idx::ObjectID{edm4hep!MCParticle}  # reference to the Monte-Carlo particle
end

function edm4hep!MCRecoTrackParticleAssociation(;weight=0, rec=-1, sim=-1)
    edm4hep!MCRecoTrackParticleAssociation(-1, weight, rec, sim)
end

function Base.getproperty(obj::edm4hep!MCRecoTrackParticleAssociation, sym::Symbol)
    if sym == :rec
        idx = getfield(obj, :rec_idx)
        return iszero(idx) ? nothing : convert(edm4hep!Track, idx)
    elseif sym == :sim
        idx = getfield(obj, :sim_idx)
        return iszero(idx) ? nothing : convert(edm4hep!MCParticle, idx)
    else # fallback to getfield
        return getfield(obj, sym)
    end
end
"""
Reconstructed Particle
- Author: F.Gaede, DESY
# Fields
- `type::Int32`: type of reconstructed particle. Check/set collection parameters ReconstructedParticleTypeNames and ReconstructedParticleTypeValues.
- `energy::Float32`: [GeV] energy of the reconstructed particle. Four momentum state is not kept consistent internally.
- `momentum::edm4hep!Vector3f`: [GeV] particle momentum. Four momentum state is not kept consistent internally.
- `referencePoint::edm4hep!Vector3f`: [mm] reference, i.e. where the particle has been measured
- `charge::Float32`: charge of the reconstructed particle.
- `mass::Float32`: [GeV] mass of the reconstructed particle, set independently from four vector. Four momentum state is not kept consistent internally.
- `goodnessOfPID::Float32`: overall goodness of the PID on a scale of [0;1]
- `covMatrix::SVector{10,Float32}`: cvariance matrix of the reconstructed particle 4vector (10 parameters). Stored as lower triangle matrix of the four momentum (px,py,pz,E), i.e. cov(px,px), cov(py,##
# Relations
- `startVertex::edm4hep!Vertex`: start vertex associated to this particle
- `particleIDUsed::edm4hep!ParticleID`: particle Id used for the kinematics of this particle
- `clusters::edm4hep!Cluster`: clusters that have been used for this particle.
- `tracks::edm4hep!Track`: tracks that have been used for this particle.
- `particles::edm4hep!ReconstructedParticle`: reconstructed particles that have been combined to this particle.
- `particleIDs::edm4hep!ParticleID`: particle Ids (not sorted by their likelihood)
# Methods
- `pushToClusters(obj::edm4hep!ReconstructedParticle, robj::edm4hep!Cluster)`: push related object to the `clusters` relation
- `popFromClusters(obj::edm4hep!ReconstructedParticle)`: pop last related object from `clusters` relation
- `pushToTracks(obj::edm4hep!ReconstructedParticle, robj::edm4hep!Track)`: push related object to the `tracks` relation
- `popFromTracks(obj::edm4hep!ReconstructedParticle)`: pop last related object from `tracks` relation
- `pushToParticles(obj::edm4hep!ReconstructedParticle, robj::edm4hep!ReconstructedParticle)`: push related object to the `particles` relation
- `popFromParticles(obj::edm4hep!ReconstructedParticle)`: pop last related object from `particles` relation
- `pushToParticleIDs(obj::edm4hep!ReconstructedParticle, robj::edm4hep!ParticleID)`: push related object to the `particleIDs` relation
- `popFromParticleIDs(obj::edm4hep!ReconstructedParticle)`: pop last related object from `particleIDs` relation
"""
struct edm4hep!ReconstructedParticle <: POD
    index::ObjectID{edm4hep!ReconstructedParticle}  # ObjectID of himself
    #---Data Members
    type::Int32                      # type of reconstructed particle. Check/set collection parameters ReconstructedParticleTypeNames and ReconstructedParticleTypeValues.
    energy::Float32                  # [GeV] energy of the reconstructed particle. Four momentum state is not kept consistent internally.
    momentum::edm4hep!Vector3f       # [GeV] particle momentum. Four momentum state is not kept consistent internally.
    referencePoint::edm4hep!Vector3f # [mm] reference, i.e. where the particle has been measured
    charge::Float32                  # charge of the reconstructed particle.
    mass::Float32                    # [GeV] mass of the reconstructed particle, set independently from four vector. Four momentum state is not kept consistent internally.
    goodnessOfPID::Float32           # overall goodness of the PID on a scale of [0;1]
    covMatrix::SVector{10,Float32}   # cvariance matrix of the reconstructed particle 4vector (10 parameters). Stored as lower triangle matrix of the four momentum (px,py,pz,E), i.e. cov(px,px), cov(py,##
    #---OneToManyRelations
    clusters::Relation{edm4hep!ReconstructedParticle,edm4hep!Cluster,1}  # clusters that have been used for this particle.
    tracks::Relation{edm4hep!ReconstructedParticle,edm4hep!Track,2}  # tracks that have been used for this particle.
    particles::Relation{edm4hep!ReconstructedParticle,edm4hep!ReconstructedParticle,3}  # reconstructed particles that have been combined to this particle.
    particleIDs::Relation{edm4hep!ReconstructedParticle,edm4hep!ParticleID,4}  # particle Ids (not sorted by their likelihood)
    #---OneToOneRelations
    startVertex_idx::ObjectID{edm4hep!Vertex}  # start vertex associated to this particle
    particleIDUsed_idx::ObjectID{edm4hep!ParticleID}  # particle Id used for the kinematics of this particle
end

function edm4hep!ReconstructedParticle(;type=0, energy=0, momentum=edm4hep!Vector3f(), referencePoint=edm4hep!Vector3f(), charge=0, mass=0, goodnessOfPID=0, covMatrix=zero(SVector{10,Float32}), clusters=Relation{edm4hep!ReconstructedParticle,edm4hep!Cluster,1}(), tracks=Relation{edm4hep!ReconstructedParticle,edm4hep!Track,2}(), particles=Relation{edm4hep!ReconstructedParticle,edm4hep!ReconstructedParticle,3}(), particleIDs=Relation{edm4hep!ReconstructedParticle,edm4hep!ParticleID,4}(), startVertex=-1, particleIDUsed=-1)
    edm4hep!ReconstructedParticle(-1, type, energy, momentum, referencePoint, charge, mass, goodnessOfPID, covMatrix, clusters, tracks, particles, particleIDs, startVertex, particleIDUsed)
end

function Base.getproperty(obj::edm4hep!ReconstructedParticle, sym::Symbol)
    if sym == :startVertex
        idx = getfield(obj, :startVertex_idx)
        return iszero(idx) ? nothing : convert(edm4hep!Vertex, idx)
    elseif sym == :particleIDUsed
        idx = getfield(obj, :particleIDUsed_idx)
        return iszero(idx) ? nothing : convert(edm4hep!ParticleID, idx)
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
function pushToParticleIDs(c::edm4hep!ReconstructedParticle, o::edm4hep!ParticleID)
    iszero(c.index) && (c = register(c))
    c = @set c.particleIDs = push(c.particleIDs, o)
    update(c)
end
function popFromParticleIDs(c::edm4hep!ReconstructedParticle)
    iszero(c.index) && (c = register(c))
    c = @set c.particleIDs = pop(c.particleIDs)
    update(c)
end
"""
Used to keep track of the correspondence between MC and reconstructed particles
- Author: C. Bernet, B. Hegner
# Fields
- `weight::Float32`: weight of this association
# Relations
- `rec::edm4hep!ReconstructedParticle`: reference to the reconstructed particle
- `sim::edm4hep!MCParticle`: reference to the Monte-Carlo particle
"""
struct edm4hep!MCRecoParticleAssociation <: POD
    index::ObjectID{edm4hep!MCRecoParticleAssociation}  # ObjectID of himself
    #---Data Members
    weight::Float32                  # weight of this association
    #---OneToOneRelations
    rec_idx::ObjectID{edm4hep!ReconstructedParticle}  # reference to the reconstructed particle
    sim_idx::ObjectID{edm4hep!MCParticle}  # reference to the Monte-Carlo particle
end

function edm4hep!MCRecoParticleAssociation(;weight=0, rec=-1, sim=-1)
    edm4hep!MCRecoParticleAssociation(-1, weight, rec, sim)
end

function Base.getproperty(obj::edm4hep!MCRecoParticleAssociation, sym::Symbol)
    if sym == :rec
        idx = getfield(obj, :rec_idx)
        return iszero(idx) ? nothing : convert(edm4hep!ReconstructedParticle, idx)
    elseif sym == :sim
        idx = getfield(obj, :sim_idx)
        return iszero(idx) ? nothing : convert(edm4hep!MCParticle, idx)
    else # fallback to getfield
        return getfield(obj, sym)
    end
end
"""
Association between a Reconstructed Particle and a Vertex
- Author: Placido Fernandez Declara
# Fields
- `weight::Float32`: weight of this association
# Relations
- `rec::edm4hep!ReconstructedParticle`: reference to the reconstructed particle
- `vertex::edm4hep!Vertex`: reference to the vertex
"""
struct edm4hep!RecoParticleVertexAssociation <: POD
    index::ObjectID{edm4hep!RecoParticleVertexAssociation}  # ObjectID of himself
    #---Data Members
    weight::Float32                  # weight of this association
    #---OneToOneRelations
    rec_idx::ObjectID{edm4hep!ReconstructedParticle}  # reference to the reconstructed particle
    vertex_idx::ObjectID{edm4hep!Vertex}  # reference to the vertex
end

function edm4hep!RecoParticleVertexAssociation(;weight=0, rec=-1, vertex=-1)
    edm4hep!RecoParticleVertexAssociation(-1, weight, rec, vertex)
end

function Base.getproperty(obj::edm4hep!RecoParticleVertexAssociation, sym::Symbol)
    if sym == :rec
        idx = getfield(obj, :rec_idx)
        return iszero(idx) ? nothing : convert(edm4hep!ReconstructedParticle, idx)
    elseif sym == :vertex
        idx = getfield(obj, :vertex_idx)
        return iszero(idx) ? nothing : convert(edm4hep!Vertex, idx)
    else # fallback to getfield
        return getfield(obj, sym)
    end
end
"""
dN/dx or dE/dx info of Track.
- Author: Wenxing Fang, IHEP
# Fields
- `dQdx::edm4hep!Quantity`: the reconstructed dEdx or dNdx and its error
- `particleType::Int16`: particle type, e(0),mu(1),pi(2),K(3),p(4).
- `type::Int16`: type.
- `hypotheses::SVector{5,edm4hep!Hypothesis}`: 5 particle hypothesis
- `hitData::PVector{edm4hep!HitLevelData}`: hit level data
# Relations
- `track::edm4hep!Track`: the corresponding track.
# Methods
- `setHitData(object::edm4hep!RecDqdx, v::AbstractVector{edm4hep!HitLevelData})`: assign a set of values to the `hitData` vector member
"""
struct edm4hep!RecDqdx <: POD
    index::ObjectID{edm4hep!RecDqdx} # ObjectID of himself
    #---Data Members
    dQdx::edm4hep!Quantity           # the reconstructed dEdx or dNdx and its error
    particleType::Int16              # particle type, e(0),mu(1),pi(2),K(3),p(4).
    type::Int16                      # type.
    hypotheses::SVector{5,edm4hep!Hypothesis}  # 5 particle hypothesis
    #---VectorMembers
    hitData::PVector{edm4hep!RecDqdx,edm4hep!HitLevelData,1}  # hit level data
    #---OneToOneRelations
    track_idx::ObjectID{edm4hep!Track}  # the corresponding track.
end

function edm4hep!RecDqdx(;dQdx=edm4hep!Quantity(), particleType=0, type=0, hypotheses=zero(SVector{5,edm4hep!Hypothesis}), hitData=PVector{edm4hep!RecDqdx,edm4hep!HitLevelData,1}(), track=-1)
    edm4hep!RecDqdx(-1, dQdx, particleType, type, hypotheses, hitData, track)
end

function Base.getproperty(obj::edm4hep!RecDqdx, sym::Symbol)
    if sym == :track
        idx = getfield(obj, :track_idx)
        return iszero(idx) ? nothing : convert(edm4hep!Track, idx)
    else # fallback to getfield
        return getfield(obj, sym)
    end
end
function setHitData(o::edm4hep!RecDqdx, v::AbstractVector{edm4hep!HitLevelData})
    iszero(o.index) && (o = register(o))
    o = @set o.hitData = v
    update(o)
end
"""
Used to keep track of the correspondence between MC and reconstructed particles
- Author: S. Joosten
# Fields
- `simID::UInt32`: Index of corresponding MCParticle (position in MCParticles array)
- `recID::UInt32`: Index of corresponding ReconstructedParticle (position in ReconstructedParticles array)
- `weight::Float32`: weight of this association
# Relations
- `rec::edm4eic!ReconstructedParticle`: reference to the reconstructed particle
- `sim::edm4hep!MCParticle`: reference to the Monte-Carlo particle
"""
struct edm4eic!MCRecoParticleAssociation <: POD
    index::ObjectID{edm4eic!MCRecoParticleAssociation}  # ObjectID of himself
    #---Data Members
    simID::UInt32                    # Index of corresponding MCParticle (position in MCParticles array)
    recID::UInt32                    # Index of corresponding ReconstructedParticle (position in ReconstructedParticles array)
    weight::Float32                  # weight of this association
    #---OneToOneRelations
    rec_idx::ObjectID{edm4eic!ReconstructedParticle}  # reference to the reconstructed particle
    sim_idx::ObjectID{edm4hep!MCParticle}  # reference to the Monte-Carlo particle
end

function edm4eic!MCRecoParticleAssociation(;simID=0, recID=0, weight=0, rec=-1, sim=-1)
    edm4eic!MCRecoParticleAssociation(-1, simID, recID, weight, rec, sim)
end

function Base.getproperty(obj::edm4eic!MCRecoParticleAssociation, sym::Symbol)
    if sym == :rec
        idx = getfield(obj, :rec_idx)
        return iszero(idx) ? nothing : convert(edm4eic!ReconstructedParticle, idx)
    elseif sym == :sim
        idx = getfield(obj, :sim_idx)
        return iszero(idx) ? nothing : convert(edm4hep!MCParticle, idx)
    else # fallback to getfield
        return getfield(obj, sym)
    end
end
"""
Collection of hits identified by the clustering algorithm to belong together
- Author: S. Joosten
# Fields
- `weights::PVector{Float32}`: Weight for each of the hits, mirrors hits array
# Relations
- `hits::edm4eic!CalorimeterHit`: Hits associated with this cluster
# Methods
- `setWeights(object::edm4eic!ProtoCluster, v::AbstractVector{Float32})`: assign a set of values to the `weights` vector member
- `pushToHits(obj::edm4eic!ProtoCluster, robj::edm4eic!CalorimeterHit)`: push related object to the `hits` relation
- `popFromHits(obj::edm4eic!ProtoCluster)`: pop last related object from `hits` relation
"""
struct edm4eic!ProtoCluster <: POD
    index::ObjectID{edm4eic!ProtoCluster}  # ObjectID of himself
    #---Data Members
    #---VectorMembers
    weights::PVector{edm4eic!ProtoCluster,Float32,1}  # Weight for each of the hits, mirrors hits array
    #---OneToManyRelations
    hits::Relation{edm4eic!ProtoCluster,edm4eic!CalorimeterHit,1}  # Hits associated with this cluster
end

function edm4eic!ProtoCluster(;weights=PVector{edm4eic!ProtoCluster,Float32,1}(), hits=Relation{edm4eic!ProtoCluster,edm4eic!CalorimeterHit,1}())
    edm4eic!ProtoCluster(-1, weights, hits)
end

function pushToHits(c::edm4eic!ProtoCluster, o::edm4eic!CalorimeterHit)
    iszero(c.index) && (c = register(c))
    c = @set c.hits = push(c.hits, o)
    update(c)
end
function popFromHits(c::edm4eic!ProtoCluster)
    iszero(c.index) && (c = register(c))
    c = @set c.hits = pop(c.hits)
    update(c)
end
function setWeights(o::edm4eic!ProtoCluster, v::AbstractVector{Float32})
    iszero(o.index) && (o = register(o))
    o = @set o.weights = v
    update(o)
end
const ProtoCluster = edm4eic!ProtoCluster
const MCRecoParticleAssociation = edm4eic!MCRecoParticleAssociation
const Track = edm4eic!Track
const Vertex = edm4eic!Vertex
const EventHeader = edm4hep!EventHeader
const MCRecoTrackParticleAssociation = edm4eic!MCRecoTrackParticleAssociation
const HadronicFinalState = edm4eic!HadronicFinalState
const MCRecoCaloAssociation = edm4hep!MCRecoCaloAssociation
const TrackerHit = edm4eic!TrackerHit
const ReconstructedParticle = edm4eic!ReconstructedParticle
const CaloHitContribution = edm4hep!CaloHitContribution
const MCRecoTrackerHitPlaneAssociation = edm4hep!MCRecoTrackerHitPlaneAssociation
const RingImage = edm4eic!RingImage
const RawTrackerHit = edm4eic!RawTrackerHit
const MCRecoCaloParticleAssociation = edm4hep!MCRecoCaloParticleAssociation
const Measurement2D = edm4eic!Measurement2D
const SimCalorimeterHit = edm4hep!SimCalorimeterHit
const RawHGCROCHit = edm4eic!RawHGCROCHit
const InclusiveKinematics = edm4eic!InclusiveKinematics
const Tensor = edm4eic!Tensor
const TrackClusterMatch = edm4eic!TrackClusterMatch
const RecoParticleVertexAssociation = edm4hep!RecoParticleVertexAssociation
const MCRecoClusterParticleAssociation = edm4eic!MCRecoClusterParticleAssociation
const CalorimeterHit = edm4eic!CalorimeterHit
const TimeSeries = edm4hep!TimeSeries
const MCRecoTrackerAssociation = edm4hep!MCRecoTrackerAssociation
const ParticleID = edm4hep!ParticleID
const PMTHit = edm4eic!PMTHit
const SimTrackerHit = edm4hep!SimTrackerHit
const TrackerHitPlane = edm4hep!TrackerHitPlane
const RecIonizationCluster = edm4hep!RecIonizationCluster
const RawCalorimeterHit = edm4hep!RawCalorimeterHit
const CherenkovParticleID = edm4eic!CherenkovParticleID
const TrackSeed = edm4eic!TrackSeed
const TrackParameters = edm4eic!TrackParameters
const TrackerPulse = edm4hep!TrackerPulse
const MCRecoCalorimeterHitAssociation = edm4eic!MCRecoCalorimeterHitAssociation
const RawTimeSeries = edm4hep!RawTimeSeries
const Cluster = edm4eic!Cluster
const MCRecoVertexParticleAssociation = edm4eic!MCRecoVertexParticleAssociation
const Trajectory = edm4eic!Trajectory
const MCParticle = edm4hep!MCParticle
const SimPrimaryIonizationCluster = edm4hep!SimPrimaryIonizationCluster
const SimPulse = edm4eic!SimPulse
const TrackSegment = edm4eic!TrackSegment
const RecDqdx = edm4hep!RecDqdx
const MCRecoTrackerHitAssociation = edm4eic!MCRecoTrackerHitAssociation
export pushToParents, popFromParents, pushToDaughters, popFromDaughters, setElectronCellID, setElectronTime, setElectronPosition, setPulseTime, setPulseAmplitude, setAdcCounts, setRawHits, setParameters, pushToClusters, popFromClusters, pushToHits, popFromHits, pushToParticleIDs, popFromParticleIDs, setShapeParameters, setSubdetectorEnergies, setAmplitude, pushToTrackerPulse, popFromTrackerPulse, setHitContributions, pushToAssociatedParticles, popFromAssociatedParticles, setShape, setFloatData, setInt64Data, setSamples, pushToContributions, popFromContributions, pushToCalorimeterHits, popFromCalorimeterHits, pushToTrackerHits, popFromTrackerHits, pushToPulses, popFromPulses, pushToParticles, popFromParticles, setWeights, pushToTrackParameters, popFromTrackParameters, pushToMeasurements_deprecated, popFromMeasurements_deprecated, pushToOutliers_deprecated, popFromOutliers_deprecated, setMeasurementChi2, setOutlierChi2, pushToMeasurements, popFromMeasurements, pushToTracks, popFromTracks, setPoints, pushToRawHitAssociations, popFromRawHitAssociations, setHypotheses, setThetaPhiPhotons, pushToHadrons, popFromHadrons, setSubdetectorHitNumbers, setTrackStates, setDxQuantities, setHitData, ProtoCluster, MCRecoParticleAssociation, Track, Vertex, EventHeader, MCRecoTrackParticleAssociation, HadronicFinalState, MCRecoCaloAssociation, TrackerHit, ReconstructedParticle, CaloHitContribution, MCRecoTrackerHitPlaneAssociation, RingImage, RawTrackerHit, MCRecoCaloParticleAssociation, Measurement2D, SimCalorimeterHit, RawHGCROCHit, InclusiveKinematics, Tensor, TrackClusterMatch, RecoParticleVertexAssociation, MCRecoClusterParticleAssociation, CalorimeterHit, TimeSeries, MCRecoTrackerAssociation, ParticleID, PMTHit, SimTrackerHit, TrackerHitPlane, RecIonizationCluster, RawCalorimeterHit, CherenkovParticleID, TrackSeed, TrackParameters, TrackerPulse, MCRecoCalorimeterHitAssociation, RawTimeSeries, Cluster, MCRecoVertexParticleAssociation, Trajectory, MCParticle, SimPrimaryIonizationCluster, SimPulse, TrackSegment, RecDqdx, MCRecoTrackerHitAssociation
