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
const SimTrackerHit = edm4hep!SimTrackerHit
const TrackerHitPlane = edm4hep!TrackerHitPlane
const Track = edm4hep!Track
const Vertex = edm4hep!Vertex
const RecIonizationCluster = edm4hep!RecIonizationCluster
const RawCalorimeterHit = edm4hep!RawCalorimeterHit
const TrackerHit = edm4hep!TrackerHit
const EventHeader = edm4hep!EventHeader
const MCRecoTrackParticleAssociation = edm4hep!MCRecoTrackParticleAssociation
const TrackerPulse = edm4hep!TrackerPulse
const MCRecoParticleAssociation = edm4hep!MCRecoParticleAssociation
const MCRecoCaloAssociation = edm4hep!MCRecoCaloAssociation
const RawTimeSeries = edm4hep!RawTimeSeries
const CaloHitContribution = edm4hep!CaloHitContribution
const MCRecoTrackerHitPlaneAssociation = edm4hep!MCRecoTrackerHitPlaneAssociation
const MCRecoCaloParticleAssociation = edm4hep!MCRecoCaloParticleAssociation
const MCParticle = edm4hep!MCParticle
const ReconstructedParticle = edm4hep!ReconstructedParticle
const SimPrimaryIonizationCluster = edm4hep!SimPrimaryIonizationCluster
const SimCalorimeterHit = edm4hep!SimCalorimeterHit
const Cluster = edm4hep!Cluster
const RecoParticleVertexAssociation = edm4hep!RecoParticleVertexAssociation
const RecDqdx = edm4hep!RecDqdx
const CalorimeterHit = edm4hep!CalorimeterHit
const TimeSeries = edm4hep!TimeSeries
const MCRecoTrackerAssociation = edm4hep!MCRecoTrackerAssociation
const ParticleID = edm4hep!ParticleID
const MCRecoClusterParticleAssociation = edm4hep!MCRecoClusterParticleAssociation
export setParameters, setAmplitude, pushToClusters, popFromClusters, pushToHits, popFromHits, pushToParticleIDs, popFromParticleIDs, setShapeParameters, setSubdetectorEnergies, pushToParents, popFromParents, pushToDaughters, popFromDaughters, setElectronCellID, setElectronTime, setElectronPosition, setPulseTime, setPulseAmplitude, pushToContributions, popFromContributions, setAdcCounts, setRawHits, pushToTrackerPulse, popFromTrackerPulse, pushToTrackerHits, popFromTrackerHits, pushToTracks, popFromTracks, setSubdetectorHitNumbers, setTrackStates, setDxQuantities, pushToParticles, popFromParticles, setHitData, SimTrackerHit, TrackerHitPlane, Track, Vertex, RecIonizationCluster, RawCalorimeterHit, TrackerHit, EventHeader, MCRecoTrackParticleAssociation, TrackerPulse, MCRecoParticleAssociation, MCRecoCaloAssociation, RawTimeSeries, CaloHitContribution, MCRecoTrackerHitPlaneAssociation, MCRecoCaloParticleAssociation, MCParticle, ReconstructedParticle, SimPrimaryIonizationCluster, SimCalorimeterHit, Cluster, RecoParticleVertexAssociation, RecDqdx, CalorimeterHit, TimeSeries, MCRecoTrackerAssociation, ParticleID, MCRecoClusterParticleAssociation
