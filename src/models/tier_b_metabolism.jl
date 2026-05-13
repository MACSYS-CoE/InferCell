# Tier B: a "v0.0.1 + ε" metabolism with two reactions LightMetabolism does not
# model. Used to generate synthetic data that the v0.0.1 fitter must absorb under
# model mismatch; not part of the inference pipeline itself.

struct TierBMetabolism <: AbstractSubModel
    params::Vector{InferParameter}
    ATP_max::Float64
    e_tx::Float64
    e_tl::Float64
    n_tx::Float64
    n_tl::Float64
    gamma_ntp::Float64
    gamma_aa::Float64
    K_M_enzyme::Float64
    # Tier B extras: a constant ATP leak and a Hill-shaped enzyme response (n_hill > 1)
    k_atp_leak::Float64
    n_hill::Float64
    mRNA_source::Symbol
    enzyme_source::Symbol
end

function TierBMetabolism(;
        k_atp=1.0, k_ntp=0.2, k_aa=0.4,
        k_tx=1.0, k_tl=2.0,
        ATP_max=10.0, e_tx=0.1, e_tl=0.05, n_tx=0.5, n_tl=0.5,
        gamma_ntp=0.1, gamma_aa=0.1, K_M_enzyme=5.0,
        k_atp_leak=0.2, n_hill=2.0,
        ATP0=5.0, NTP0=2.0, AA0=2.0,
        sigma_obs=0.3,
        mRNA_source::Symbol=:mRNA,
        enzyme_source::Symbol=:protein)
    params = [
        InferParameter(k_atp, LogNormal(0, 1), false, :k_atp, :tierb, :rate),
        InferParameter(k_ntp, LogNormal(0, 1), false, :k_ntp, :tierb, :rate),
        InferParameter(k_aa,  LogNormal(0, 1), false, :k_aa,  :tierb, :rate),
        InferParameter(k_tx,  LogNormal(0, 1), false, :k_tx,  :tierb, :rate),
        InferParameter(k_tl,  LogNormal(0, 1), false, :k_tl,  :tierb, :rate),
        InferParameter(ATP0,  Normal(5, 1),    true,  :ATP0,   :tierb, :initial_condition),
        InferParameter(NTP0,  Normal(2, 1),    true,  :NTP0,   :tierb, :initial_condition),
        InferParameter(AA0,   Normal(2, 1),    true,  :AA0,    :tierb, :initial_condition),
        InferParameter(sigma_obs, truncated(Normal(0, 1); lower=0.0), false, :sigma_obs, :tierb, :observation),
    ]
    return TierBMetabolism(params, ATP_max, e_tx, e_tl, n_tx, n_tl, gamma_ntp, gamma_aa,
                           K_M_enzyme, k_atp_leak, n_hill, mRNA_source, enzyme_source)
end

states(::TierBMetabolism) = [:ATP, :NTP, :AA]
parameters(m::TierBMetabolism) = m.params
formalism(::TierBMetabolism) = :ode
inference_mode(::TierBMetabolism) = :differentiable
inputs(m::TierBMetabolism) = [m.mRNA_source, m.enzyme_source]

function dynamics(u_local, p_local, t, m::TierBMetabolism, u_inputs)
    ATP, NTP, AA = u_local
    k_atp, k_ntp, k_aa, k_tx, k_tl = p_local
    mRNA, enzyme = u_inputs

    # Hill response (n_hill > 1 makes the enzyme→rate response sharper than MM).
    # LightMetabolism uses Hill with n=1, so this is one source of mismatch.
    enzyme_n = enzyme^m.n_hill
    K_n = m.K_M_enzyme^m.n_hill
    enzyme_factor = enzyme_n / (K_n + enzyme_n)
    v_atp_prod = k_atp * (1 + enzyme_factor) * (m.ATP_max - ATP)
    v_ntp_synth = k_ntp * ATP
    v_aa_synth = k_aa * ATP

    # Second source of mismatch: a constant ATP leak LightMetabolism doesn't model.
    dATP = v_atp_prod - m.e_tx * k_tx - m.e_tl * k_tl * mRNA -
           v_ntp_synth - v_aa_synth - m.k_atp_leak
    dNTP = v_ntp_synth - m.n_tx * k_tx - m.gamma_ntp * NTP
    dAA  = v_aa_synth  - m.n_tl * k_tl * mRNA - m.gamma_aa * AA

    SA[dATP, dNTP, dAA]
end
