using EDM4hep
using EDM4hep.RootIO
using EDM4hep.Analysis
using EDM4hep.SystemOfUnits
using Plots
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

using DataFrames

mutable struct AnalysisData <: AbstractAnalysisData
    df::DataFrame
    pevts::Int64
    sevts::Int64
    AnalysisData() = new(DataFrame(Zcand_m = Float32[], Zcand_recoil_m = Float32[], Zcand_q = Int32[], Zcand_recoil_θ = Float32[]), 0, 0)
end
# Need to tell how to merge two DataFramnes
Base.merge!(df1::DataFrame, df2::DataFrame) = append!(df1, df2)

# f = "root://eospublic.cern.ch//eos/experiment/fcc/ee/generation/DelphesEvents/winter2023/IDEA/p8_ee_ZZ_ecm240/events_000189367.root"
f = joinpath(@__DIR__, "../../../examples/FCC/events_000189367.root")

reader = RootIO.Reader(f);
events = RootIO.get(reader, "events");

reader

function myanalysis!(data::AnalysisData, reader, events)
    for evt in events
        data.pevts += 1
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
            Zcand_recoil_θ = zed_leptonic_recoil[1].momentum |> EDM4hep.θ
            Zcand_q        = zed_leptonic[1].charge
            if 80GeV <= Zcand_m <= 100GeV
                push!(data.df, (Zcand_m, Zcand_recoil_m, Zcand_q, Zcand_recoil_θ))
                data.sevts += 1
            end
        end
    end
    return data
end

N = Threads.nthreads()
data = AnalysisData();

elapsed1 = @elapsed do_analysis!(data, myanalysis!, reader, events; mt=false)
println("Serial: total time: $elapsed1, $(data.pevts/elapsed1) events/s. Selected events: $(data.sevts)")

elapsed2 = @elapsed do_analysis!(data, myanalysis!, reader, events; mt=true)
println("MT[$N]: total time: $elapsed2, $(data.pevts/elapsed2) events/s. Selected events: $(data.sevts)")
println("Speeedup: $(elapsed1/elapsed2)")

data = AnalysisData()
elapsed1 = @elapsed do_analysis!(data, myanalysis!, reader, events; mt=false)
println("Serial: total time: $elapsed1, $(data.pevts/elapsed1) events/s. Selected events: $(data.sevts)")

elapsed2 = @elapsed do_analysis!(data, myanalysis!, reader, events; mt=true)
println("MT[$N]: total time: $elapsed2, $(data.pevts/elapsed2) events/s. Selected events: $(data.sevts)")
println("Speeedup: $(elapsed1/elapsed2)")

histogram(data.df.Zcand_m, title="Resonance mass plot",xlabel="GeV")

histogram(data.df.Zcand_recoil_m, title="Z leptonic recoil",xlabel="GeV")
