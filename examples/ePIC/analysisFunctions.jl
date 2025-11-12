using LorentzVectorHEP
using EDM4hep

function invariantmass(p1::ReconstructedParticle, p2::ReconstructedParticle)
    lv1 = LorentzVector(p1.energy, p1.momentum...)
    lv2 = LorentzVector(p2.energy, p2.momentum...)
    return mass(lv1 + lv2)
end
LorentzVectorHEP.eta(p::MCParticle) = 0.5 * log((norm(p.momentum) + p.momentum.z) / (norm(p.momentum) - p.momentum.z))
LorentzVectorHEP.phi(p::MCParticle) = atan(p.momentum.y, p.momentum.x)

function Base.getproperty(s::StructArray{MCParticle}, key::Symbol)
    if key === :eta; return eta.(s)
    elseif key === :energy; return getproperty.(s, :energy)
    else
        return getfield(getfield(s, :component), key)
    end
end