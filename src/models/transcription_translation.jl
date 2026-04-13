struct TranscriptionTranslation <: AbstractSubModel
    params::Vector{InferParameter}
end

function TranscriptionTranslation(;
        k_tx=1.0, k_tl=2.0, gamma_mRNA=0.5, gamma_protein=0.1,
        mRNA0=0.0, protein0=0.0, sigma_obs=0.3)
    params = [
        InferParameter(k_tx,          LogNormal(0, 1),                     false, :k_tx,          :txl, :rate),
        InferParameter(k_tl,          LogNormal(0, 1),                     false, :k_tl,          :txl, :rate),
        InferParameter(gamma_mRNA,    LogNormal(0, 1),                     false, :gamma_mRNA,    :txl, :rate),
        InferParameter(gamma_protein, LogNormal(0, 1),                     false, :gamma_protein, :txl, :rate),
        InferParameter(mRNA0,         Normal(0, 1),                        true,  :mRNA0,         :txl, :initial_condition),
        InferParameter(protein0,      Normal(0, 1),                        true,  :protein0,      :txl, :initial_condition),
        InferParameter(sigma_obs,     truncated(Normal(0, 1); lower=0.0),  false, :sigma_obs,     :txl, :observation),
    ]
    return TranscriptionTranslation(params)
end

states(::TranscriptionTranslation) = [:mRNA, :protein]
parameters(m::TranscriptionTranslation) = m.params
formalism(::TranscriptionTranslation) = :ode
inference_mode(::TranscriptionTranslation) = :differentiable

function dynamics(u, p, t, ::TranscriptionTranslation)
    mRNA, protein = u
    k_tx, k_tl, gamma_mRNA, gamma_protein = p
    SA[k_tx - gamma_mRNA * mRNA, k_tl * mRNA - gamma_protein * protein]
end
