"""
    LightMetabolism(; k_atp, k_ntp, k_aa, k_tx, k_tl, ATP_max, …, mRNA_source=:mRNA, enzyme_source=:protein)

Lightweight ODE metabolism block: ATP production (gated by enzyme via
Michaelis-Menten with constant `K_M_enzyme`), NTP and amino-acid synthesis,
and TX/TL-driven sinks. States `[:ATP, :NTP, :AA]`. Couples to a gene-
expression block by reading `mRNA_source` (drives TX consumption) and
`enzyme_source` (modulates ATP production). Closes the v0.0.1 autocatalytic
loop with [`TranscriptionTranslation`](@ref).
"""
struct LightMetabolism <: AbstractSubModel
    params::Vector{InferParameter}
    ATP_max::Float64
    e_tx::Float64
    e_tl::Float64
    n_tx::Float64
    n_tl::Float64
    gamma_ntp::Float64
    gamma_aa::Float64
    K_M_enzyme::Float64
    mRNA_source::Symbol
    enzyme_source::Symbol
end

function LightMetabolism(;
        k_atp=1.0, k_ntp=0.2, k_aa=0.4,
        k_tx=1.0, k_tl=2.0,
        ATP_max=10.0, e_tx=0.1, e_tl=0.05, n_tx=0.5, n_tl=0.5,
        gamma_ntp=0.1, gamma_aa=0.1, K_M_enzyme=5.0,
        ATP0=5.0, NTP0=2.0, AA0=2.0,
        sigma_obs=0.3,
        mRNA_source::Symbol=:mRNA,
        enzyme_source::Symbol=:protein)
    params = [
        InferParameter(k_atp, LogNormal(0, 1), false, :k_atp, :metab, :rate),
        InferParameter(k_ntp, LogNormal(0, 1), false, :k_ntp, :metab, :rate),
        InferParameter(k_aa,  LogNormal(0, 1), false, :k_aa,  :metab, :rate),
        InferParameter(k_tx,  LogNormal(0, 1), false, :k_tx,  :metab, :rate),
        InferParameter(k_tl,  LogNormal(0, 1), false, :k_tl,  :metab, :rate),
        InferParameter(ATP0,  Normal(5, 1),    true,  :ATP0,   :metab, :initial_condition),
        InferParameter(NTP0,  Normal(2, 1),    true,  :NTP0,   :metab, :initial_condition),
        InferParameter(AA0,   Normal(2, 1),    true,  :AA0,    :metab, :initial_condition),
        InferParameter(sigma_obs, truncated(Normal(0, 1); lower=0.0), false, :sigma_obs, :metab, :observation),
    ]
    return LightMetabolism(params, ATP_max, e_tx, e_tl, n_tx, n_tl, gamma_ntp, gamma_aa,
                           K_M_enzyme, mRNA_source, enzyme_source)
end

states(::LightMetabolism) = [:ATP, :NTP, :AA]
parameters(m::LightMetabolism) = m.params
formalism(::LightMetabolism) = :ode
inference_mode(::LightMetabolism) = :differentiable
inputs(m::LightMetabolism) = [m.mRNA_source, m.enzyme_source]

function dynamics(u_local, p_local, t, m::LightMetabolism, u_inputs)
    ATP, NTP, AA = u_local
    k_atp, k_ntp, k_aa, k_tx, k_tl = p_local
    mRNA, enzyme = u_inputs

    # Michaelis-Menten enzyme modulation: factor in [0, 1), bounded.
    # enzyme=0 preserves the Step-3 baseline rate; enzyme→∞ doubles it.
    enzyme_factor = enzyme / (m.K_M_enzyme + enzyme)
    v_atp_prod = k_atp * (1 + enzyme_factor) * (m.ATP_max - ATP)
    v_ntp_synth = k_ntp * ATP
    v_aa_synth = k_aa * ATP

    dATP = v_atp_prod - m.e_tx * k_tx - m.e_tl * k_tl * mRNA - v_ntp_synth - v_aa_synth
    dNTP = v_ntp_synth - m.n_tx * k_tx - m.gamma_ntp * NTP
    dAA  = v_aa_synth - m.n_tl * k_tl * mRNA - m.gamma_aa * AA

    SA[dATP, dNTP, dAA]
end
