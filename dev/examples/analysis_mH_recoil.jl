using EDM4hep
using EDM4hep.RootIO
using EDM4hep.SystemOfUnits
using EDM4hep.Histograms
using Plots; gr()
theme(:boxed)

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

# f = "root://eospublic.cern.ch//eos/experiment/fcc/ee/generation/DelphesEvents/winter2023/IDEA/p8_ee_ZZ_ecm240/events_000189367.root"
f = joinpath(@__DIR__, "../../../examples/FCC/events_000189367.root")

reader = RootIO.Reader(f);
events = RootIO.get(reader, "events");

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

plot((plot(h) for h in myhists)..., layout=(5,2), size=(1200,1500))
