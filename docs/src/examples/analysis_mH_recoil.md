```@meta
EditURL = "analysis_mH_recoil.lit"
```

# FCCee Analysis Examples

Using the example `higgs/mH-recoil/mumu` from [FCCAnalyses](https://github.com/HEP-FCC/FCCAnalyses)

!!! note "Note that"
    You can also download this example as a
    [Jupyter notebook](analysis_mH_recoil.ipynb) and a plain
    [Julia source file](analysis_mH_recoil.jl).

#### Table of contents
```@contents
Pages = ["analysis_mH_recoil.md"]
Depth = 2:3
```

## Load the necessary modules

````julia
using EDM4hep
using EDM4hep.RootIO
using EDM4hep.SystemOfUnits
using EDM4hep.Histograms
using Plots; gr()
theme(:boxed)
````

````
Precompiling packages...
   5267.1 ms  ✓ EDM4hep → EDM4hepPlotsExt
  1 dependency successfully precompiled in 6 seconds. 268 already precompiled.

````

## Definition of some analysis functions
These are couple of examples of high-level functions that makes use of `ReconstructedParticle`
objects to build resonances and recoils.
They make use of standard Julia functions to generate combinations, to sort a vector,
and to work with LorentzVectors.

re-using convenient existing packages

````julia
using LorentzVectorHEP
using Combinatorics

"""
    resonanceBuilder(rmass::AbstractFloat, legs::AbstractVector{ReconstructedParticle})

Returns a container with the best resonance of 2 by 2 combinatorics of the `legs` container
sorted by closest to the input `rmass` in absolute value.
"""
function resonanceBuilder(rmass::AbstractFloat, legs::AbstractVector{ReconstructedParticle})
    result = ReconstructedParticle[]
    length(legs) < 2 && return result
    for (a,b) in combinations(legs, 2)
        lv = LorentzVector(a.energy, a.momentum...) + LorentzVector(b.energy, b.momentum...)
        rcharge = a.charge + b.charge
        push!(result, ReconstructedParticle(mass=mass(lv), momentum=(lv.x, lv.y, lv.z), charge=rcharge))
    end
    sort!(result, lt =  (a,b) -> abs(rmass-a.mass) < abs(rmass-b.mass))
    return result[1:1]  # take the best one
end;

"""
    recoilBuilder(comenergy::AbstractFloat, legs::AbstractVector{ReconstructedParticle})

    build the recoil from an arbitrary list of input `ReconstructedParticle`s and the center of mass energy.
"""
function recoilBuilder(comenergy::AbstractFloat, in::AbstractVector{ReconstructedParticle})
    result = ReconstructedParticle[]
    isempty(in) && return result
    recoil_lv = LorentzVector(comenergy, 0, 0, 0)
    for p in in
        recoil_lv -= LorentzVector(p.mass, p.momentum...)
    end
    push!(result, ReconstructedParticle(mass=mass(recoil_lv), momentum=(recoil_lv.x, recoil_lv.y, recoil_lv.z)))
    return result
end;
````

## Defining the histograms
We create a custom structure with all the histograms intialized with their binning,
units and titles. We use and the way of plotting them.
We use the module `Parameters` that allows to create user structures with defaults.

````julia
myhists = (
    mz          = H1D("m_{Z}",125,0,250, unit=:GeV),
    mz_zoom     = H1D("m_{Z} (zoom)",40,80,100, unit=:GeV),
    lr_m        = H1D("Z leptonic recoil", 100, 0, 200, unit=:GeV),
    lr_m_zoom   = H1D("Z leptonic recoil (zoom)", 200, 80, 160, unit=:GeV),
    lr_m_zoom1  = H1D("Z leptonic recoil (zoom)", 100, 120, 140, unit=:GeV),
    lr_m_zoom2  = H1D("Z leptonic recoil (zoom)", 200, 120, 140, unit=:GeV),
    lr_m_zoom3  = H1D("Z leptonic recoil (zoom)", 400, 120, 140, unit=:GeV),
    lr_m_zoom4  = H1D("Z leptonic recoil (zoom)", 800, 120, 140, unit=:GeV),
    lr_m_zoom5  = H1D("Z leptonic recoil (zoom)", 100, 130.3, 132.5, unit=:GeV),
    mz_lr_m     = H2D("m_{Z} vs Z leptonic recoil", 40, 80, 100, 100, 120, 140, units=(:GeV, :GeV)),
);
````

## Open the data file to get the events
- It is using a file in EOS with the `root:` protocol
- The obtained `events` is a `LazyTree` created by the [UnROOT.jl](https://github.com/JuliaHEP/UnROOT.jl) package.
  As the name indicates, the event is actually yet read.

````julia
# f = "root://eospublic.cern.ch//eos/experiment/fcc/ee/generation/DelphesEvents/winter2023/IDEA/p8_ee_ZZ_ecm240/events_000189367.root"
f = joinpath(@__DIR__, "../../../examples/FCC/events_000189367.root")

reader = RootIO.Reader(f);
events = RootIO.get(reader, "events");
````

## Loop over events and fill the histograms

````julia
@time for evt in events
    muids = RootIO.get(reader, evt, "Muon#0")     # get the ObjectIDs of Muons
    length(muids) < 2 && continue                 # need at least 2 muons to build a Z
    recps = RootIO.get(reader, evt, "ReconstructedParticles")
    muons = recps[muids]                          # use the objectIDs to collect the referenced ReconstructedParticles
    sel_muons = filter(x -> pₜ(x) > 10GeV, muons)  # select muons with pT > 10 GeV
    zed_leptonic = resonanceBuilder(91GeV, sel_muons)
    zed_leptonic_recoil = recoilBuilder(240GeV, zed_leptonic)
    if length(zed_leptonic) == 1                  #  Filter to have exactly one Z candidate
        Zcand_m        = zed_leptonic[1].mass
        Zcand_recoil_m = zed_leptonic_recoil[1].mass
        Zcand_q        = zed_leptonic[1].charge
        if 80GeV <= Zcand_m <= 100GeV
            ##---Fill histograms now--------------------------------------
            push!(myhists.mz, Zcand_m)
            push!(myhists.mz_zoom, Zcand_m)
            push!(myhists.lr_m, Zcand_recoil_m)
            push!(myhists.lr_m_zoom, Zcand_recoil_m)
            push!(myhists.lr_m_zoom1, Zcand_recoil_m)
            push!(myhists.lr_m_zoom2, Zcand_recoil_m)
            push!(myhists.lr_m_zoom3, Zcand_recoil_m)
            push!(myhists.lr_m_zoom4, Zcand_recoil_m)
            push!(myhists.lr_m_zoom5, Zcand_recoil_m)
            push!(myhists.mz_lr_m, Zcand_m, Zcand_recoil_m)
        end
    end
end
````

````
 11.443717 seconds (24.92 M allocations: 4.544 GiB, 8.94% gc time, 63.08% compilation time)

````

## Plot the results

````julia
plot((plot(h) for h in myhists)..., layout=(5,2), size=(1200,1500))
````

````


ERROR: FieldError: type EDM4hep.Histograms.H2D has no field `usym`, available fields: `title`, `hist`, `unit`, `uval`
````

---

*This page was generated using [Literate.jl](https://github.com/fredrikekre/Literate.jl).*

