var documenterSearchIndex = {"docs":
[{"location":"api/#Public-Documentation","page":"Public APIs","title":"Public Documentation","text":"","category":"section"},{"location":"api/","page":"Public APIs","title":"Public APIs","text":"Documentation for EDM4hep.jl public interface.","category":"page"},{"location":"api/#Modules","page":"Public APIs","title":"Modules","text":"","category":"section"},{"location":"api/","page":"Public APIs","title":"Public APIs","text":"Modules = [EDM4hep, EDM4hep.RootIO]\nOrder = [:module]","category":"page"},{"location":"api/#EDM4hep.EDM4hep","page":"Public APIs","title":"EDM4hep.EDM4hep","text":"Main module for EDM4hep.jl – Key4hep Event Data Model for Julia.\n\nAll data model types are exported from this module for public use\n\nExports\n\n\n\n\n\n","category":"module"},{"location":"api/#EDM4hep.RootIO","page":"Public APIs","title":"EDM4hep.RootIO","text":"ROOT I/O module for EDM4hep.jl\n\nIt supports both formats: TTree and RNTuple\n\n\n\n\n\n","category":"module"},{"location":"api/#Types","page":"Public APIs","title":"Types","text":"","category":"section"},{"location":"api/","page":"Public APIs","title":"Public APIs","text":"This is the list of all types defined for EDM4hep using the PODIO yaml file.","category":"page"},{"location":"api/","page":"Public APIs","title":"Public APIs","text":"Modules = [EDM4hep, EDM4hep.RootIO]\nOrder = [:type]","category":"page"},{"location":"api/#EDM4hep.CaloHitContribution","page":"Public APIs","title":"EDM4hep.CaloHitContribution","text":"Monte Carlo contribution to SimCalorimeterHit\n\nAuthor: F.Gaede, DESY\n\nFields\n\nPDG::Int32: PDG code of the shower particle that caused this contribution.\nenergy::Float32: energy in [GeV] of the this contribution\ntime::Float32: time in [ns] of this contribution\nstepPosition::Vector3f: position of this energy deposition (step) [mm]\n\nRelations\n\nparticle::MCParticle: primary MCParticle that caused the shower responsible for this contribution to the hit.\n\n\n\n\n\n","category":"type"},{"location":"api/#EDM4hep.CalorimeterHit","page":"Public APIs","title":"EDM4hep.CalorimeterHit","text":"Calorimeter hit\n\nAuthor: F.Gaede, DESY\n\nFields\n\ncellID::UInt64: detector specific (geometrical) cell id.\nenergy::Float32: energy of the hit in [GeV].\nenergyError::Float32: error of the hit energy in [GeV].\ntime::Float32: time of the hit in [ns].\nposition::Vector3f: position of the hit in world coordinates in [mm].\ntype::Int32: type of hit. Mapping of integer types to names via collection parameters \"CalorimeterHitTypeNames\" and \"CalorimeterHitTypeValues\".\n\n\n\n\n\n","category":"type"},{"location":"api/#EDM4hep.Cluster","page":"Public APIs","title":"EDM4hep.Cluster","text":"Calorimeter Hit Cluster\n\nAuthor: F.Gaede, DESY\n\nFields\n\ntype::Int32: flagword that defines the type of cluster. Bits 16-31 are used internally.\nenergy::Float32: energy of the cluster [GeV]\nenergyError::Float32: error on the energy\nposition::Vector3f: position of the cluster [mm]\npositionError::SVector{6,Float32}: covariance matrix of the position (6 Parameters)\niTheta::Float32: intrinsic direction of cluster at position  Theta. Not to be confused with direction cluster is seen from IP.\nphi::Float32: intrinsic direction of cluster at position - Phi. Not to be confused with direction cluster is seen from IP.\ndirectionError::Vector3f: covariance matrix of the direction (3 Parameters) [mm^2]\nshapeParameters::Float32: shape parameters - check/set collection parameter ClusterShapeParameters for size and names of parameters.\nsubdetectorEnergies::Float32: energy observed in a particular subdetector. Check/set collection parameter ClusterSubdetectorNames for decoding the indices of the array.\n\nRelations\n\nclusters::Cluster: clusters that have been combined to this cluster.\nhits::CalorimeterHit: hits that have been combined to this cluster.\nparticleIDs::ParticleID: particle IDs (sorted by their likelihood)\n\n\n\n\n\n","category":"type"},{"location":"api/#EDM4hep.EventHeader","page":"Public APIs","title":"EDM4hep.EventHeader","text":"Event Header. Additional parameters are assumed to go into the metadata tree.\n\nAuthor: F.Gaede\n\nFields\n\neventNumber::Int32: event number\nrunNumber::Int32: run number\ntimeStamp::UInt64: time stamp\nweight::Float32: event weight\n\n\n\n\n\n","category":"type"},{"location":"api/#EDM4hep.HitLevelData","page":"Public APIs","title":"EDM4hep.HitLevelData","text":"Fields\n\ncellID::UInt64: cell id\nN::UInt32: number of reconstructed ionization cluster.\neDep::Float32: reconstructed energy deposit [GeV].\npathLength::Float32: track path length [mm].\n\n\n\n\n\n","category":"type"},{"location":"api/#EDM4hep.Hypothesis","page":"Public APIs","title":"EDM4hep.Hypothesis","text":"Fields\n\nchi2::Float32: chi2\nexpected::Float32: expected value\nsigma::Float32: sigma value\n\n\n\n\n\n","category":"type"},{"location":"api/#EDM4hep.MCParticle","page":"Public APIs","title":"EDM4hep.MCParticle","text":"The Monte Carlo particle - based on the lcio::MCParticle.\n\nAuthor: F.Gaede, DESY\n\nFields\n\nPDG::Int32: PDG code of the particle\ngeneratorStatus::Int32: status of the particle as defined by the generator\nsimulatorStatus::Int32: status of the particle from the simulation program - use BIT constants below\ncharge::Float32: particle charge\ntime::Float32: creation time of the particle in [ns] wrt. the event, e.g. for preassigned decays or decays in flight from the simulator.\nmass::Float64: mass of the particle in [GeV]\nvertex::Vector3d: production vertex of the particle in [mm].\nendpoint::Vector3d: endpoint of the particle in [mm]\nmomentum::Vector3d: particle 3-momentum at the production vertex in [GeV]\nmomentumAtEndpoint::Vector3d: particle 3-momentum at the endpoint in [GeV]\nspin::Vector3f: spin (helicity) vector of the particle.\ncolorFlow::Vector2i: color flow as defined by the generator\n\nRelations\n\nparents::MCParticle: The parents of this particle.\ndaughters::MCParticle: The daughters this particle.\n\n\n\n\n\n","category":"type"},{"location":"api/#EDM4hep.MCRecoCaloAssociation","page":"Public APIs","title":"EDM4hep.MCRecoCaloAssociation","text":"Association between a CaloHit and the corresponding simulated CaloHit\n\nAuthor: C. Bernet, B. Hegner\n\nFields\n\nweight::Float32: weight of this association\n\nRelations\n\nrec::CalorimeterHit: reference to the reconstructed hit\nsim::SimCalorimeterHit: reference to the simulated hit\n\n\n\n\n\n","category":"type"},{"location":"api/#EDM4hep.MCRecoCaloParticleAssociation","page":"Public APIs","title":"EDM4hep.MCRecoCaloParticleAssociation","text":"Association between a CalorimeterHit and a MCParticle\n\nAuthor: Placido Fernandez Declara\n\nFields\n\nweight::Float32: weight of this association\n\nRelations\n\nrec::CalorimeterHit: reference to the reconstructed hit\nsim::MCParticle: reference to the Monte-Carlo particle\n\n\n\n\n\n","category":"type"},{"location":"api/#EDM4hep.MCRecoClusterParticleAssociation","page":"Public APIs","title":"EDM4hep.MCRecoClusterParticleAssociation","text":"Association between a Cluster and a MCParticle\n\nAuthor: Placido Fernandez Declara\n\nFields\n\nweight::Float32: weight of this association\n\nRelations\n\nrec::Cluster: reference to the cluster\nsim::MCParticle: reference to the Monte-Carlo particle\n\n\n\n\n\n","category":"type"},{"location":"api/#EDM4hep.MCRecoParticleAssociation","page":"Public APIs","title":"EDM4hep.MCRecoParticleAssociation","text":"Used to keep track of the correspondence between MC and reconstructed particles\n\nAuthor: C. Bernet, B. Hegner\n\nFields\n\nweight::Float32: weight of this association\n\nRelations\n\nrec::ReconstructedParticle: reference to the reconstructed particle\nsim::MCParticle: reference to the Monte-Carlo particle\n\n\n\n\n\n","category":"type"},{"location":"api/#EDM4hep.MCRecoTrackParticleAssociation","page":"Public APIs","title":"EDM4hep.MCRecoTrackParticleAssociation","text":"Association between a Track and a MCParticle\n\nAuthor: Placido Fernandez Declara\n\nFields\n\nweight::Float32: weight of this association\n\nRelations\n\nrec::Track: reference to the track\nsim::MCParticle: reference to the Monte-Carlo particle\n\n\n\n\n\n","category":"type"},{"location":"api/#EDM4hep.MCRecoTrackerAssociation","page":"Public APIs","title":"EDM4hep.MCRecoTrackerAssociation","text":"Association between a TrackerHit and the corresponding simulated TrackerHit\n\nAuthor: C. Bernet, B. Hegner\n\nFields\n\nweight::Float32: weight of this association\n\nRelations\n\nrec::TrackerHit: reference to the reconstructed hit\nsim::SimTrackerHit: reference to the simulated hit\n\n\n\n\n\n","category":"type"},{"location":"api/#EDM4hep.MCRecoTrackerHitPlaneAssociation","page":"Public APIs","title":"EDM4hep.MCRecoTrackerHitPlaneAssociation","text":"Association between a TrackerHitPlane and the corresponding simulated TrackerHit\n\nAuthor: Placido Fernandez Declara\n\nFields\n\nweight::Float32: weight of this association\n\nRelations\n\nrec::TrackerHitPlane: reference to the reconstructed hit\nsim::SimTrackerHit: reference to the simulated hit\n\n\n\n\n\n","category":"type"},{"location":"api/#EDM4hep.ParticleID","page":"Public APIs","title":"EDM4hep.ParticleID","text":"ParticleID\n\nAuthor: F.Gaede, DESY\n\nFields\n\ntype::Int32: userdefined type\nPDG::Int32: PDG code of this id - ( 999999 ) if unknown.\nalgorithmType::Int32: type of the algorithm/module that created this hypothesis\nlikelihood::Float32: likelihood of this hypothesis - in a user defined normalization.\nparameters::Float32: parameters associated with this hypothesis. Check/set collection parameters ParameterNames_PIDAlgorithmTypeName for decoding the indices.\n\n\n\n\n\n","category":"type"},{"location":"api/#EDM4hep.Quantity","page":"Public APIs","title":"EDM4hep.Quantity","text":"Fields\n\ntype::Int16: flag identifying how to interpret the quantity\nvalue::Float32: value of the quantity\nerror::Float32: error on the value of the quantity\n\n\n\n\n\n","category":"type"},{"location":"api/#EDM4hep.RawCalorimeterHit","page":"Public APIs","title":"EDM4hep.RawCalorimeterHit","text":"Raw calorimeter hit\n\nAuthor: F.Gaede, DESY\n\nFields\n\ncellID::UInt64: detector specific (geometrical) cell id.\namplitude::Int32: amplitude of the hit in ADC counts.\ntimeStamp::Int32: time stamp for the hit.\n\n\n\n\n\n","category":"type"},{"location":"api/#EDM4hep.RawTimeSeries","page":"Public APIs","title":"EDM4hep.RawTimeSeries","text":"Raw data of a detector readout\n\nAuthor: F.Gaede, DESY\n\nFields\n\ncellID::UInt64: detector specific cell id.\nquality::Int32: quality flag for the hit.\ntime::Float32: time of the hit [ns].\ncharge::Float32: integrated charge of the hit [fC].\ninterval::Float32: interval of each sampling [ns].\nadcCounts::Int32: raw data (32-bit) word at i.\n\n\n\n\n\n","category":"type"},{"location":"api/#EDM4hep.RecDqdx","page":"Public APIs","title":"EDM4hep.RecDqdx","text":"dN/dx or dE/dx info of Track.\n\nAuthor: Wenxing Fang, IHEP\n\nFields\n\ndQdx::Quantity: the reconstructed dEdx or dNdx and its error\nparticleType::Int16: particle type, e(0),mu(1),pi(2),K(3),p(4).\ntype::Int16: type.\nhypotheses::SVector{5,Hypothesis}: 5 particle hypothesis\nhitData::HitLevelData: hit level data\n\nRelations\n\ntrack::Track: the corresponding track.\n\n\n\n\n\n","category":"type"},{"location":"api/#EDM4hep.RecIonizationCluster","page":"Public APIs","title":"EDM4hep.RecIonizationCluster","text":"Reconstructed Ionization Cluster\n\nAuthor: Wenxing Fang, IHEP\n\nFields\n\ncellID::UInt64: cell id.\nsignificance::Float32: significance.\ntype::Int16: type.\n\nRelations\n\ntrackerPulse::TrackerPulse: the TrackerPulse used to create the ionization cluster.\n\n\n\n\n\n","category":"type"},{"location":"api/#EDM4hep.RecoParticleVertexAssociation","page":"Public APIs","title":"EDM4hep.RecoParticleVertexAssociation","text":"Association between a Reconstructed Particle and a Vertex\n\nAuthor: Placido Fernandez Declara\n\nFields\n\nweight::Float32: weight of this association\n\nRelations\n\nrec::ReconstructedParticle: reference to the reconstructed particle\nvertex::Vertex: reference to the vertex\n\n\n\n\n\n","category":"type"},{"location":"api/#EDM4hep.ReconstructedParticle","page":"Public APIs","title":"EDM4hep.ReconstructedParticle","text":"Reconstructed Particle\n\nAuthor: F.Gaede, DESY\n\nFields\n\ntype::Int32: type of reconstructed particle. Check/set collection parameters ReconstructedParticleTypeNames and ReconstructedParticleTypeValues.\nenergy::Float32: [GeV] energy of the reconstructed particle. Four momentum state is not kept consistent internally.\nmomentum::Vector3f: [GeV] particle momentum. Four momentum state is not kept consistent internally.\nreferencePoint::Vector3f: [mm] reference, i.e. where the particle has been measured\ncharge::Float32: charge of the reconstructed particle.\nmass::Float32: [GeV] mass of the reconstructed particle, set independently from four vector. Four momentum state is not kept consistent internally.\ngoodnessOfPID::Float32: overall goodness of the PID on a scale of [0;1]\ncovMatrix::SVector{10,Float32}: cvariance matrix of the reconstructed particle 4vector (10 parameters). Stored as lower triangle matrix of the four momentum (px,py,pz,E), i.e. cov(px,px), cov(py,##\n\nRelations\n\nstartVertex::Vertex: start vertex associated to this particle\nparticleIDUsed::ParticleID: particle Id used for the kinematics of this particle\nclusters::Cluster: clusters that have been used for this particle.\ntracks::Track: tracks that have been used for this particle.\nparticles::ReconstructedParticle: reconstructed particles that have been combined to this particle.\nparticleIDs::ParticleID: particle Ids (not sorted by their likelihood)\n\n\n\n\n\n","category":"type"},{"location":"api/#EDM4hep.SimCalorimeterHit","page":"Public APIs","title":"EDM4hep.SimCalorimeterHit","text":"Simulated calorimeter hit\n\nAuthor: F.Gaede, DESY\n\nFields\n\ncellID::UInt64: ID of the sensor that created this hit\nenergy::Float32: energy of the hit in [GeV].\nposition::Vector3f: position of the hit in world coordinates in [mm].\n\nRelations\n\ncontributions::CaloHitContribution: Monte Carlo step contribution - parallel to particle\n\n\n\n\n\n","category":"type"},{"location":"api/#EDM4hep.SimPrimaryIonizationCluster","page":"Public APIs","title":"EDM4hep.SimPrimaryIonizationCluster","text":"Simulated Primary Ionization\n\nAuthor: Wenxing Fang, IHEP\n\nFields\n\ncellID::UInt64: cell id.\ntime::Float32: the primary ionization's time in the lab frame [ns].\nposition::Vector3d: the primary ionization's position [mm].\ntype::Int16: type.\nelectronCellID::UInt64: cell id.\nelectronTime::Float32: the time in the lab frame [ns].\nelectronPosition::Vector3d: the position in the lab frame [mm].\npulseTime::Float32: the pulse's time in the lab frame [ns].\npulseAmplitude::Float32: the pulse's amplitude [fC].\n\nRelations\n\nmcparticle::MCParticle: the particle that caused the ionizing collisions.\n\n\n\n\n\n","category":"type"},{"location":"api/#EDM4hep.SimTrackerHit","page":"Public APIs","title":"EDM4hep.SimTrackerHit","text":"Simulated tracker hit\n\nAuthor: F.Gaede, DESY\n\nFields\n\ncellID::UInt64: ID of the sensor that created this hit\nEDep::Float32: energy deposited in the hit [GeV].\ntime::Float32: proper time of the hit in the lab frame in [ns].\npathLength::Float32: path length of the particle in the sensitive material that resulted in this hit.\nquality::Int32: quality bit flag.\nposition::Vector3d: the hit position in [mm].\nmomentum::Vector3f: the 3-momentum of the particle at the hits position in [GeV]\n\nRelations\n\nmcparticle::MCParticle: MCParticle that caused the hit.\n\n\n\n\n\n","category":"type"},{"location":"api/#EDM4hep.TimeSeries","page":"Public APIs","title":"EDM4hep.TimeSeries","text":"Calibrated Detector Data\n\nAuthor: Wenxing Fang, IHEP\n\nFields\n\ncellID::UInt64: cell id.\ntime::Float32: begin time [ns].\ninterval::Float32: interval of each sampling [ns].\namplitude::Float32: calibrated detector data.\n\n\n\n\n\n","category":"type"},{"location":"api/#EDM4hep.Track","page":"Public APIs","title":"EDM4hep.Track","text":"Reconstructed track\n\nAuthor: F.Gaede, DESY\n\nFields\n\ntype::Int32: flagword that defines the type of track.Bits 16-31 are used internally\nchi2::Float32: Chi^2 of the track fit\nndf::Int32: number of degrees of freedom of the track fit\ndEdx::Float32: dEdx of the track.\ndEdxError::Float32: error of dEdx.\nradiusOfInnermostHit::Float32: radius of the innermost hit that has been used in the track fit\nsubdetectorHitNumbers::Int32: number of hits in particular subdetectors.Check/set collection variable TrackSubdetectorNames for decoding the indices\ntrackStates::TrackState: track states\ndxQuantities::Quantity: different measurements of dx quantities\n\nRelations\n\ntrackerHits::TrackerHit: hits that have been used to create this track\ntracks::Track: tracks (segments) that have been combined to create this track\n\n\n\n\n\n","category":"type"},{"location":"api/#EDM4hep.TrackState","page":"Public APIs","title":"EDM4hep.TrackState","text":"Fields\n\nlocation::Int32: for use with At{Other|IP|FirstHit|LastHit|Calorimeter|Vertex}|LastLocation\nD0::Float32: transverse impact parameter\nphi::Float32: azimuthal angle\nomega::Float32: is the signed curvature of the track in [1/mm].\nZ0::Float32: longitudinal impact parameter\ntanLambda::Float32: lambda is the dip angle of the track in r-z\ntime::Float32: time of the track at this trackstate\nreferencePoint::Vector3f: Reference point of the track parameters, e.g. the origin at the IP, or the position  of the first/last hits or the entry point into the calorimeter. [mm]\ncovMatrix::SVector{21,Float32}: lower triangular covariance matrix of the track parameters.  the order of parameters is  d0, phi, omega, z0, tan(lambda), time. the array is a row-major flattening of the matrix.\n\n\n\n\n\n","category":"type"},{"location":"api/#EDM4hep.TrackerHit","page":"Public APIs","title":"EDM4hep.TrackerHit","text":"Tracker hit\n\nAuthor: F.Gaede, DESY\n\nFields\n\ncellID::UInt64: ID of the sensor that created this hit\ntype::Int32: type of raw data hit, either one of edm4hep::RawTimeSeries, edm4hep::SIMTRACKERHIT - see collection parameters \"TrackerHitTypeNames\" and \"TrackerHitTypeValues\".\nquality::Int32: quality bit flag of the hit.\ntime::Float32: time of the hit [ns].\neDep::Float32: energy deposited on the hit [GeV].\neDepError::Float32: error measured on EDep [GeV].\nposition::Vector3d: hit position in [mm].\ncovMatrix::SVector{6,Float32}: covariance of the position (x,y,z), stored as lower triangle matrix. i.e. cov(x,x) , cov(y,x) , cov(y,y) , cov(z,x) , cov(z,y) , cov(z,z)\nrawHits::ObjectID: raw data hits. Check getType to get actual data type.\n\n\n\n\n\n","category":"type"},{"location":"api/#EDM4hep.TrackerHitPlane","page":"Public APIs","title":"EDM4hep.TrackerHitPlane","text":"Tracker hit plane\n\nAuthor: Placido Fernandez Declara, CERN\n\nFields\n\ncellID::UInt64: ID of the sensor that created this hit\ntype::Int32: type of raw data hit, either one of edm4hep::RawTimeSeries, edm4hep::SIMTRACKERHIT - see collection parameters \"TrackerHitTypeNames\" and \"TrackerHitTypeValues\".\nquality::Int32: quality bit flag of the hit.\ntime::Float32: time of the hit [ns].\neDep::Float32: energy deposited on the hit [GeV].\neDepError::Float32: error measured on EDep [GeV].\nu::Vector2f: measurement direction vector, u lies in the x-y plane\nv::Vector2f: measurement direction vector, v is along z\ndu::Float32: measurement error along the direction\ndv::Float32: measurement error along the direction\nposition::Vector3d: hit position in [mm].\ncovMatrix::SVector{6,Float32}: covariance of the position (x,y,z), stored as lower triangle matrix. i.e. cov(x,x) , cov(y,x) , cov(y,y) , cov(z,x) , cov(z,y) , cov(z,z)\nrawHits::ObjectID: raw data hits. Check getType to get actual data type.\n\n\n\n\n\n","category":"type"},{"location":"api/#EDM4hep.TrackerPulse","page":"Public APIs","title":"EDM4hep.TrackerPulse","text":"Reconstructed Tracker Pulse\n\nAuthor: Wenxing Fang, IHEP\n\nFields\n\ncellID::UInt64: cell id.\ntime::Float32: time [ns].\ncharge::Float32: charge [fC].\nquality::Int16: quality.\ncovMatrix::SVector{3,Float32}: lower triangle covariance matrix of the charge(c) and time(t) measurements.\n\nRelations\n\ntimeSeries::TimeSeries: Optionally, the timeSeries that has been used to create the pulse can be stored with the pulse.\n\n\n\n\n\n","category":"type"},{"location":"api/#EDM4hep.Vector2f","page":"Public APIs","title":"EDM4hep.Vector2f","text":"Fields\n\na::Float32: \nb::Float32: \n\n\n\n\n\n","category":"type"},{"location":"api/#EDM4hep.Vector2i","page":"Public APIs","title":"EDM4hep.Vector2i","text":"Fields\n\na::Int32: \nb::Int32: \n\n\n\n\n\n","category":"type"},{"location":"api/#EDM4hep.Vector3d","page":"Public APIs","title":"EDM4hep.Vector3d","text":"Fields\n\nx::Float64: \ny::Float64: \nz::Float64: \n\n\n\n\n\n","category":"type"},{"location":"api/#EDM4hep.Vector3f","page":"Public APIs","title":"EDM4hep.Vector3f","text":"Fields\n\nx::Float32: \ny::Float32: \nz::Float32: \n\n\n\n\n\n","category":"type"},{"location":"api/#EDM4hep.Vector4f","page":"Public APIs","title":"EDM4hep.Vector4f","text":"Generic vector for storing classical 4D coordinates in memory. Four momentum helper functions are in edm4hep::utils\n\nFields\n\nx::Float32: \ny::Float32: \nz::Float32: \nt::Float32: \n\n\n\n\n\n","category":"type"},{"location":"api/#EDM4hep.Vertex","page":"Public APIs","title":"EDM4hep.Vertex","text":"Vertex\n\nAuthor: F.Gaede, DESY\n\nFields\n\nprimary::Int32: boolean flag, if vertex is the primary vertex of the event\nchi2::Float32: chi-squared of the vertex fit\nprobability::Float32: probability of the vertex fit\nposition::Vector3f: [mm] position of the vertex.\ncovMatrix::SVector{6,Float32}: covariance matrix of the position (stored as lower triangle matrix, i.e. cov(xx),cov(y,x),cov(z,x),cov(y,y),... )\nalgorithmType::Int32: type code for the algorithm that has been used to create the vertex - check/set the collection parameters AlgorithmName and AlgorithmType.\nparameters::Float32: additional parameters related to this vertex - check/set the collection parameter \"VertexParameterNames\" for the parameters meaning.\n\nRelations\n\nassociatedParticle::POD: reconstructed particle associated to this vertex.\n\n\n\n\n\n","category":"type"},{"location":"api/#EDM4hep.RootIO.Reader","page":"Public APIs","title":"EDM4hep.RootIO.Reader","text":"The Reader struture keeps a reference to the UnROOT LazyTree and caches already built 'layouts' of the EDM4hep types. The layouts maps a set of columns in the LazyTree into an object.\n\n\n\n\n\n","category":"type"},{"location":"api/#Functions","page":"Public APIs","title":"Functions","text":"","category":"section"},{"location":"api/","page":"Public APIs","title":"Public APIs","text":"Modules = [EDM4hep, EDM4hep.RootIO]\nOrder = [:function]","category":"page"},{"location":"api/#EDM4hep.RootIO.get-Tuple{EDM4hep.RootIO.Reader, String}","page":"Public APIs","title":"EDM4hep.RootIO.get","text":"get(reader::Reader, treename::String)\n\nOpens a 'TTree' in the ROOT file (typically the events tree).  It returns a 'LazyTree' that allows the user to iterate over events. \n\n\n\n\n\n","category":"method"},{"location":"api/#EDM4hep.RootIO.get-Tuple{EDM4hep.RootIO.Reader, UnROOT.LazyEvent, String}","page":"Public APIs","title":"EDM4hep.RootIO.get","text":"get(reader::Reader, evt::UnROOT.LazyEvent, bname::String; btype::Type=Any, register=true)\n\nGets an object collection by its name, with the possibility to overwrite the mapping Julia type or use the  type known in the ROOT file (C++ class name). The optonal key parameter register indicates is the collection needs to be registered to the EDStore.\n\n\n\n\n\n","category":"method"},{"location":"#EDM4hep-in-Julia","page":"Introduction","title":"EDM4hep in Julia","text":"","category":"section"},{"location":"","page":"Introduction","title":"Introduction","text":"Prototype of the EDM4hep (generic Event Data Model for HEP experiments part of Key4hep) for Julia with the goal to have very simple structures (isbits) with the purpose to evaluate its ergonomic design and implementation performance.","category":"page"},{"location":"#PODIO-generation","page":"Introduction","title":"PODIO generation","text":"","category":"section"},{"location":"","page":"Introduction","title":"Introduction","text":"The Julia POD structs should be generated from the the edm4hep.yaml yaml file using PODIO scripts.","category":"page"},{"location":"","page":"Introduction","title":"Introduction","text":"For the time being, for this evaluation, the Julia structs are generated by a local Julia script in ./podio/generate.jl. The files genComponents/jl and genDatatypes.jl are fully generated. ","category":"page"},{"location":"#Main-design-features","page":"Introduction","title":"Main design features","text":"","category":"section"},{"location":"","page":"Introduction","title":"Introduction","text":"All EDM4hep entities are immutable structs and containing basic types and structs. Including the relationships (one-to-one and one-to-many) and vector members. Objects attributes cannot be changed. This makes the objects with the isbits(obj) == true.\nObjects are created by default not being registered, they are free floating. The user can register them with register(obj).\nNote that operations like register, setting  relationships (add_daughter,...), etc. will automatically create a new instance. The typical pattern is to overwrite the user variable with the new instance, e.g.:\np1 = MCParticle(...)\np1 = register(p1)\np1, d1 = add_daughter(p1, MCParticle(...))\nThe main goal for reading EDM4hep containers from a ROOT file is to obtain as result a StructArray(see StructArrays.jl documentation). This provides a very efficient access by column and the same time provide a convenient views as object instances. For example if you want to sum the momentum of a range of MCParticles the user should be able to write:\njulia> mcparticles = get(reader,...)\njulia> mcparticles[5:8].momentum\n4-element StructArray(::Vector{Float32}, ::Vector{Float32}, ::Vector{Float32}) with eltype Vector3f:\n  (0.8740664,-0.002116337,124.84335)\n  (0.8602309,-0.056633994,-124.632545)\n  (-0.00012483617,0.0021162117,0.0026654897)\n  (0.014539731,0.05663412,-0.32755017)\n\njulia> sum(mcparticles[5:8].momentum)\n  (1.7487122,0.0,-0.11407688)","category":"page"},{"location":"#Roadmap","page":"Introduction","title":"Roadmap","text":"","category":"section"},{"location":"","page":"Introduction","title":"Introduction","text":"There are a number of issues and problems still to be revolved. We keep track of them in this list:","category":"page"},{"location":"","page":"Introduction","title":"Introduction","text":"~~Need to generate and support VectorMembers. For the time being are ignored.~~\n~~Handle cyclic datatype dependencies. In EDM4hep there is one case that is not yet resolved. Vertex depends on ReconstructedParticle in a one-to-one relation and ReconstructedParticle relates to Vertex. Using the abstract class POD in this case works well to break the cycle.~~ \n~~Better handle collectionID in one-to-many relations~~\n~~Be able to read RNTuple files in addition to TTree files~~\nSupport latest version (RC2) of RNTuple format (waiting for a file being generated)\nGenerate accessors for one-to-many relations, vector members\nGenerate doc string with member information","category":"page"},{"location":"#Tests","page":"Introduction","title":"Tests","text":"","category":"section"},{"location":"","page":"Introduction","title":"Introduction","text":"Unit tests can be run with julia --project=. test/runtests.jl","category":"page"},{"location":"#Examples","page":"Introduction","title":"Examples","text":"","category":"section"},{"location":"#examples/mcparticle_tree.jl","page":"Introduction","title":"examples/mcparticle_tree.jl","text":"","category":"section"},{"location":"","page":"Introduction","title":"Introduction","text":"Creates a collection of MCParticles and a collection of SimTrackerHits in memory, constructing the relations between particle parents and daughters, as well as, the the one-to-one relation between the simulation hit to the originator MCParticle.","category":"page"},{"location":"#examples/read_example.jl","page":"Introduction","title":"examples/read_example.jl","text":"","category":"section"},{"location":"","page":"Introduction","title":"Introduction","text":"An example of reading from a ROOT file created by the C++ implementation of EDM4hep. This is the full code:","category":"page"},{"location":"","page":"Introduction","title":"Introduction","text":"using EDM4hep\nusing EDM4hep.RootIO\n\ncd(@__DIR__)\nf = \"ttbar_edm4hep_digi.root\"\n\nreader = RootIO.Reader(f)\nevents = RootIO.get(reader, \"events\")\n\nevt = events[1];\n\nhits = RootIO.get(reader, evt, \"InnerTrackerBarrelCollection\")\nmcps = RootIO.get(reader, evt, \"MCParticle\")\n\nfor hit in hits\n    println(\"Hit $(hit.index) is related to MCParticle $(hit.mcparticle.index) with name $(hit.mcparticle.name)\")\nend\n\nfor p in mcps\n    println(\"MCParticle $(p.index) $(p.name) with momentum $(p.momentum) and energy $(p.energy) has $(length(p.daughters)) daughters\")\n    for d in p.daughters\n        println(\"   ---> $(d.index) $(d.name) and momentum $(d.momentum) has $(length(d.parents)) parents\")\n        for m in d.parents\n            println(\"      ---> $(m.index) $(m.name)\")\n        end \n    end\nend","category":"page"},{"location":"#EDM4hep-Data-Model","page":"Introduction","title":"EDM4hep Data Model","text":"","category":"section"},{"location":"","page":"Introduction","title":"Introduction","text":"This is the diagram for the EDM4hep datamodel including relationships. ","category":"page"},{"location":"","page":"Introduction","title":"Introduction","text":"(Image: Figure)","category":"page"},{"location":"#Release-Notes","page":"Introduction","title":"Release Notes","text":"","category":"section"},{"location":"#0.1.1-In-preparation","page":"Introduction","title":"0.1.1 - In preparation","text":"","category":"section"}]
}