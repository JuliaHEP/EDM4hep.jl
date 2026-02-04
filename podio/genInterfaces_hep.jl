"""
Tracker hit interface class
- Author: Thomas Madlener, DESY
# Fields
- `cellID::UInt64`:  ID of the sensor that created this hit
- `type::Int32`:  type of the raw data hit
- `quality::Int32`:  quality bit flag of the hit
- `time::Float32`:  time of the hit
- `eDep::Float32`:  energy deposited on the hit
- `eDepError::Float32`:  error measured on eDep
- `position::edm4hep!Vector3d`:  hit position as recorded by the sensor. The exact interpretation will depend on the currently held type of the interface
"""
abstract type edm4hep!TrackerHit <: POD
end


#---Exports
export edm4hep!TrackerHit
