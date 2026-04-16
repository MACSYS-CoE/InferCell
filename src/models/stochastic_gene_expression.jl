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
    # mRNA production: 0 -> mRNA (zeroth order)
    tx_rate(u, p, t) = p[1]         # k_tx
    tx_affect!(integrator) = (integrator.u[1] += 1)
    tx = ConstantRateJump(tx_rate, tx_affect!)

    # mRNA degradation: mRNA -> 0
    mrna_deg_rate(u, p, t) = p[3] * u[1]  # gamma_mRNA * mRNA
    mrna_deg_affect!(integrator) = (integrator.u[1] -= 1)
    mrna_deg = ConstantRateJump(mrna_deg_rate, mrna_deg_affect!)

    # Translation: mRNA -> mRNA + protein
    tl_rate(u, p, t) = p[2] * u[1]  # k_tl * mRNA
    tl_affect!(integrator) = (integrator.u[2] += 1)
    tl = ConstantRateJump(tl_rate, tl_affect!)

    # Protein degradation: protein -> 0
    prot_deg_rate(u, p, t) = p[4] * u[2]  # gamma_protein * protein
    prot_deg_affect!(integrator) = (integrator.u[2] -= 1)
    prot_deg = ConstantRateJump(prot_deg_rate, prot_deg_affect!)

    return [tx, mrna_deg, tl, prot_deg]
end
