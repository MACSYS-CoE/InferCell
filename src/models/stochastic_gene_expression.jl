"""
    StochasticGeneExpression(; k_tx, k_tl, gamma_mRNA, gamma_protein, mRNA0, protein0)

Constitutive stochastic transcription + translation as a jump process: four
reactions (mRNA production at constant rate `k_tx`, mRNA degradation,
translation, protein degradation). States `[:mRNA, :protein]`. Inference uses
ABC-SMC; this block has no ODE equivalent below the law-of-large-numbers limit.
"""
struct StochasticGeneExpression <: AbstractSubModel
    params::Vector{InferParameter}
end

function StochasticGeneExpression(;
        k_tx=1.0, k_tl=2.0, gamma_mRNA=0.5, gamma_protein=0.1,
        mRNA0=0, protein0=0)
    params = [
        InferParameter(k_tx,          LogNormal(0, 1), false, :k_tx,          :sge, :rate),
        InferParameter(k_tl,          LogNormal(0, 1), false, :k_tl,          :sge, :rate),
        InferParameter(gamma_mRNA,    LogNormal(0, 1), false, :gamma_mRNA,    :sge, :rate),
        InferParameter(gamma_protein, LogNormal(0, 1), false, :gamma_protein, :sge, :rate),
        InferParameter(Float64(mRNA0),    Normal(0, 1), true, :mRNA0,         :sge, :initial_condition),
        InferParameter(Float64(protein0), Normal(0, 1), true, :protein0,      :sge, :initial_condition),
    ]
    return StochasticGeneExpression(params)
end

states(::StochasticGeneExpression) = [:mRNA, :protein]
parameters(m::StochasticGeneExpression) = m.params
formalism(::StochasticGeneExpression) = :jump
inference_mode(::StochasticGeneExpression) = :simulation

function reactions(::StochasticGeneExpression)
    # Local coordinates: u = [mRNA, protein]; p = [k_tx, k_tl, gamma_mRNA, gamma_protein].
    # mRNA production: 0 -> mRNA (zeroth order)
    tx = Reaction((u, p, t, _) -> p[1], (u, _) -> (u[1] += 1))
    # mRNA degradation: mRNA -> 0
    mrna_deg = Reaction((u, p, t, _) -> p[3] * u[1], (u, _) -> (u[1] -= 1))
    # Translation: mRNA -> mRNA + protein
    tl = Reaction((u, p, t, _) -> p[2] * u[1], (u, _) -> (u[2] += 1))
    # Protein degradation: protein -> 0
    prot_deg = Reaction((u, p, t, _) -> p[4] * u[2], (u, _) -> (u[2] -= 1))

    return [tx, mrna_deg, tl, prot_deg]
end
