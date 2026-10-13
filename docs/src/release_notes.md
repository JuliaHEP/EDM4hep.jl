
# Release Notes

## Unreleased
- Added support for the EDM4eic data model (schema version 830), ported from the `v0-patches` branch
  - Select it with `EDM4hep.set_edmodel("eic")` (or the `EDModel=eic` environment variable) before loading the package
  - Generated `podio/gen*_eic.jl` files from `podio/edm4eic.yaml` merged with the EDM4hep 1.0 schema
  - The ROOT reader recognises `edm4eic::` classes and skips branches whose type is not in the selected model
  - Added an ePIC example and the `EICModelReader` test; CI runs the tests for both models
- Documentation restructured as in `v0-patches`: tutorials under `docs/src/tutorials`, examples (FCC and ePIC) under `docs/src/examples`, generated with Literate by `docs/literate.jl`. The FCC example pages are the ones rendered with version 0.4, since their input sample uses an EDM4hep schema older than 1.0
- Added the `EDM4hepPlotsExt` package extension with Plots recipes for `H1D` and `H2D` histograms
- `RootIO.create_getter` registers the collection in the `EDStore` (keyword `register=true`) so that relations can be navigated, and handles one-to-one relations in files written with newer podio versions
- Fixed `H1D`, `H2D` and `H3D` constructors for the current FHist version

## 1.0.0
- New EDM4hep schema version 1 not backward compatible
  - Relation between `ReconstructedParticle` to ``ParticleID` inverted
  - Names and comments of attributes have been revised
  - Link entities are now implemented using a parametric type `Link{FROM,TO}`
  - Introduced `GeneratorEventParameters` for MC generators
    covMatrix::CovMatrix2f           #  covariance matrix of the charge and time measurements 
  - New types for covariant matrices (e.g. `CovMatrix2f`)
  - Introduced type interfaces (e.g. `TrackerHit`) that types can inherit from (e.g. `TrackerHit3D`)
- Introduced `EDCollection{ED}` to represent a collection of data of a given EDM4hep type. It is always implemented as a `StructArray` type to provide an effcient "SoA" access to the data. Added, as well, a few functions to provide an iterable and indexing behavior to the collections, and to control their lifetime (e.g. `getEDCollection`, `hasEDCollection`, `initEDCollection`, `emptyEDStore`). It replaces `EDStore` type and their functions.
- Added more tests, in particular reading files exercising all defined types and their relations.
- Added two new tutorials as part of the generated documentation. One to exercise the EDM4hep model (EDM), and one for reading ROOT data file (I/O).
- Adapted to use new `events___CollectionTypeInfo` in `podio_metadata`
- Use the version 1.0 of EDM4hep yaml file
- Prepared for supporting other edm models (e.g. eic). Using Preferences to select the edmodel
- Update latest version of UnROOT

## 0.4.2
- Moved to JuliaHEP organization
- Update to FHist 0.11 series (#9)

## 0.4.1
- Use StructArrays also for EDStore when building the model directly in memory.

## 0.4.0
### New Functionality
- Added function `RootIO.create_getter(reader::Reader, bname::String; selection=nothing)` to create a getter function for a specific branch.
  The optional argument allows to select leaves to be read.
- The overall performance is highly improved (factor 3 with respect previous version)

## 0.3.1
### Bug Fixes
- Legacy podio test fixed

## 0.3.0
### New Functionality
- Optimisations by explicitly generation of StructArrays
- Support for multi-files
- Support for multi-threading
- Support for RNTuple RC2
- Added analysis module ROOT.Analysis (mini-framework for MT analysis)

## 0.2.0
### New Functionality
- Added support for older versions of PODIO file formats
- Added analysis example for FCCAnalysis

## 0.1.0
- First release. It provides the basic functionality to create in memory an event model and be able to read from a ROOT file (supporting both TTree and RNTuple formats).
