"""
    BurstyGeneExpression(; k_on, k_off, k_tx_burst, k_tl, gamma_mRNA, gamma_protein, promoter0, mRNA0, protein0)

Telegraph-promoter stochastic gene expression: a two-state promoter (`:promoter
∈ {0, 1}`) gates mRNA production at rate `k_tx_burst`, producing visible
mRNA-count bursting (Fano factor > 1). States `[:promoter, :mRNA, :protein]`,
six jump reactions. Used as the v0.0.1 SSA block; complements
[`StochasticGeneExpression`](@ref) (which is constitutive).
"""
struct BurstyGeneExpression <: AbstractSubModel
    params::Vector{InferParameter}
end

function BurstyGeneExpression(;
        k_on=0.2, k_off=0.5, k_tx_burst=10.0, k_tl=2.0,
        gamma_mRNA=0.5, gamma_protein=0.1,
        promoter0=0, mRNA0=0, protein0=0)
    params = [
        InferParameter(k_on,          LogNormal(0, 1), false, :k_on,          :bge, :rate),
        InferParameter(k_off,         LogNormal(0, 1), false, :k_off,         :bge, :rate),
        InferParameter(k_tx_burst,    LogNormal(0, 1), false, :k_tx_burst,    :bge, :rate),
        InferParameter(k_tl,          LogNormal(0, 1), false, :k_tl,          :bge, :rate),
        InferParameter(gamma_mRNA,    LogNormal(0, 1), false, :gamma_mRNA,    :bge, :rate),
        InferParameter(gamma_protein, LogNormal(0, 1), false, :gamma_protein, :bge, :rate),
        InferParameter(Float64(promoter0), Normal(0, 1), true, :promoter0, :bge, :initial_condition),
        InferParameter(Float64(mRNA0),     Normal(0, 1), true, :mRNA0,     :bge, :initial_condition),
        InferParameter(Float64(protein0),  Normal(0, 1), true, :protein0,  :bge, :initial_condition),
    ]
    return BurstyGeneExpression(params)
end

states(::BurstyGeneExpression) = [:promoter, :mRNA, :protein]
parameters(m::BurstyGeneExpression) = m.params
formalism(::BurstyGeneExpression) = :jump
inference_mode(::BurstyGeneExpression) = :simulation

function reactions(::BurstyGeneExpression)
    # Parameter order (free, deduplicated): k_on, k_off, k_tx_burst, k_tl, gamma_mRNA, gamma_protein
    # State order: promoter (u[1]), mRNA (u[2]), protein (u[3])

    # Promoter OFF -> ON: fires only when promoter == 0
    on_rate(u, p, t) = p[1] * (1 - u[1])
    on_affect!(integrator) = (integrator.u[1] = 1)
    rxn_on = ConstantRateJump(on_rate, on_affect!)

    # Promoter ON -> OFF: fires only when promoter == 1
    off_rate(u, p, t) = p[2] * u[1]
    off_affect!(integrator) = (integrator.u[1] = 0)
    rxn_off = ConstantRateJump(off_rate, off_affect!)

    # mRNA production: gated by promoter state
    tx_rate(u, p, t) = p[3] * u[1]
    tx_affect!(integrator) = (integrator.u[2] += 1)
    rxn_tx = ConstantRateJump(tx_rate, tx_affect!)

    # mRNA degradation
    mrna_deg_rate(u, p, t) = p[5] * u[2]
    mrna_deg_affect!(integrator) = (integrator.u[2] -= 1)
    rxn_mrna_deg = ConstantRateJump(mrna_deg_rate, mrna_deg_affect!)

    # Translation
    tl_rate(u, p, t) = p[4] * u[2]
    tl_affect!(integrator) = (integrator.u[3] += 1)
    rxn_tl = ConstantRateJump(tl_rate, tl_affect!)

    # Protein degradation
    prot_deg_rate(u, p, t) = p[6] * u[3]
    prot_deg_affect!(integrator) = (integrator.u[3] -= 1)
    rxn_prot_deg = ConstantRateJump(prot_deg_rate, prot_deg_affect!)

    return [rxn_on, rxn_off, rxn_tx, rxn_mrna_deg, rxn_tl, rxn_prot_deg]
end
