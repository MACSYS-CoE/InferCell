struct LightMetabolism <: AbstractSubModel
    params::Vector{InferParameter}
    ATP_max::Float64
    e_tx::Float64
    e_tl::Float64
    n_tx::Float64
    n_tl::Float64
    gamma_ntp::Float64
    gamma_aa::Float64
end

function LightMetabolism(;
        k_atp=1.0, k_ntp=0.2, k_aa=0.4,
        k_tx=1.0, k_tl=2.0,
        ATP_max=10.0, e_tx=0.1, e_tl=0.05, n_tx=0.5, n_tl=0.5,
        gamma_ntp=0.1, gamma_aa=0.1,
        ATP0=5.0, NTP0=2.0, AA0=2.0,
        sigma_obs=0.3)
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
    return LightMetabolism(params, ATP_max, e_tx, e_tl, n_tx, n_tl, gamma_ntp, gamma_aa)
end

states(::LightMetabolism) = [:ATP, :NTP, :AA]
parameters(m::LightMetabolism) = m.params
formalism(::LightMetabolism) = :ode
inference_mode(::LightMetabolism) = :differentiable
inputs(::LightMetabolism) = [:mRNA]

function dynamics(u_local, p_local, t, m::LightMetabolism, u_inputs)
    ATP, NTP, AA = u_local
    k_atp, k_ntp, k_aa, k_tx, k_tl = p_local
    mRNA = u_inputs[1]

    v_atp_prod = k_atp * (m.ATP_max - ATP)
    v_ntp_synth = k_ntp * ATP
    v_aa_synth = k_aa * ATP

    dATP = v_atp_prod - m.e_tx * k_tx - m.e_tl * k_tl * mRNA - v_ntp_synth - v_aa_synth
    dNTP = v_ntp_synth - m.n_tx * k_tx - m.gamma_ntp * NTP
    dAA  = v_aa_synth - m.n_tl * k_tl * mRNA - m.gamma_aa * AA

    SA[dATP, dNTP, dAA]
end
