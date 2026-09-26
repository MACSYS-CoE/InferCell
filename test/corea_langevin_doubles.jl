using InferCell
using LinearAlgebra
using Random
using StaticArrays: SVector
using ForwardDiff
using Statistics: quantile

# Check 1b's chemical-Langevin ensemble (spec §3, §11 task 14.2). A test-local
# double: the `:sde` formalism stays refused in `src` and §7's non-goal stands
# (§12, 2026-09-25 B).
#
# **What it simulates.** The ODE block's reactions as a chemical Langevin
# equation, with the jump path frozen at a reference run's. "Frozen" means
# every handshake of the reference is replayed exactly as it acted on the ODE
# block: the dilution its growth applied, the enzyme concentrations and radius
# it wrote into the rate laws, and the debits and credits it paid. Only the
# metabolic reactions fluctuate.
#
# **How the path is frozen.** At handshake k the reference driver diluted the
# ODE state by `dilₖ`, wrote `θₖ`, integrated one interval and paid its
# counters. The recorder reruns the interval deterministically with this
# file's own integrator, from the diluted previous state, and stores
#
#     aₖ = u_postₖ − Φ₁(dilₖ ⊙ u_postₖ₋₁; θₖ)
#
# the part of the handshake no rate law explains. Replay is `u ← dilₖ ⊙ u`,
# one interval of the Langevin integrator at `θₖ`, then `u += aₖ`. With the
# noise off, the replay is the reference, by construction. That puts the
# integrator's own error into `aₖ`. So the self-test is that `aₖ` is small on
# every species no counter touches, the glycolytic intermediates and redox,
# where it measures how well this integrator reproduces the pinned Rodas5P
# over one interval.
#
# **The integrator: local linearization, not Euler–Maruyama** (§12,
# 2026-09-26). The block's fastest mode is about −2.1e4 /s, so explicit
# Euler–Maruyama needs steps under 5e-5 s, and a drift-implicit step damps the
# noise on exactly the fast near-equilibrium modes that set the small pools.
# Each step here is instead the exact solution of the SDE linearised at the
# step's start:
#
#     u ← u + h·φ₁(hJ)·f(u) + η,   η ~ N(0, ∫₀ʰ e^{Js} D e^{Jᵀs} ds)
#     D = Σ_c N_c N_cᵀ |v_c| / Ω
#
# The sum is over **channels**, each reaction's forward and reverse separately,
# so a near-equilibrium reaction carries its gross noise and not its small net.
# `Ω` is the particles per mM at the handshake. The covariance integral is
# built by scaling and squaring: a third-order series over a step short enough
# that `‖J‖·h₀ ≤ 10⁻³`, then doubled with `Σ(2t) = F(t)Σ(t)F(t)ᵀ + Σ(t)`. Van
# Loan's block exponential would need `e^{−Jh}`, which overflows at `|λ|h ≈
# 1000`.
#
# **Small pools take no noise: the channels are partitioned** (§12,
# 2026-09-26 G). At each step's start, a channel keeps its noise only if
# every intracellular species it touches holds at least `min_particles`,
# `CONTINUUM_MEDIAN` by default. The others stay in the drift, noiselessly
# (Salis & Kaznessis 2005's partition, by population alone). Without it, a
# pool of a few particles under a large gross flux, GMP under GK1 or phospho-EI
# in the cascade, draws below zero on about half its steps. The clamp that
# follows creates mass: in job 17558677 the guanylate total rose from 1.98 to
# 22 mM over 1,200 handshakes, GK1 spent the injected GMP's phosphate from
# ATP, the energy charge collapsed, and a negative ATP gave the linearisation
# a +1,625 /s mode that overflowed `exp`. The clamp is kept as a last resort,
# and what each clamp injects is booked and reported, as is each moiety's
# drift from the reference, which the injection bounds. Most of what is booked
# is not drift: where the reference's own debit clipped a pool to zero, the
# clamp after `aₖ` returns the replay to that zero.
#
# **Channels are read off the right-hand side.** Each channel is one parameter
# the composed `f` is linear in: a `kcatF` or `kcatR`, a PTS `kf` or `kr`,
# `p_lact2r` or `k_chg`. `g_c = p_c·∂f/∂p_c = N_c·v_c`. The integer column
# `N_c` is recovered once and its integrality asserted, so no stoichiometry is
# restated here.

# The ODE block composed without a jump block. The enzymes are the nominal
# ones on the structs, and each handshake's enzyme concentration is applied by
# scaling both catalytic constants by `E/E_nominal`, which is the same rate
# law. The PTS scalars are freed as channels. `r_cell_nm` is always free.
const _LGLYC = Tuple(r.id for r in GLYCOLYTIC_REACTIONS)
const _LPTS = (:glcpts0, :glcpts1, :glcpts2, :glcpts3, :glcpts4)
_lkcats(ids) = Symbol[s for r in ids for s in (Symbol("kcatF_", r), Symbol("kcatR_", r))]

"""
    langevin_ode_models() -> Vector{AbstractSubModel}

The four ODE modules as the assembled model composes them, at nominal enzymes
with every channel parameter freed, plus phase 6's held protein counts.
"""
langevin_ode_models() = AbstractSubModel[
    CentralGlycolysis(free = _lkcats(_LGLYC)),
    PtsTransport(free = vcat([Symbol(k, "_", s) for s in _LPTS for k in (:kf, :kr)], [:p_lact2r])),
    NucleotideRecycling(free = _lkcats(RECYCLING_REACTIONS)),
    TrnaCharging(free_k_chg = true),
    HeldProteins(),
]

"""
    LangevinBlock

The ODE-only composition and what the integrator needs from it:
- `f(u, p)`, the composed right-hand side;
- the state names and the map from the hybrid's ODE layout into this one;
- the parameter vector at nominal values;
- the channel slots and their integer stoichiometries;
- the enzyme slots each handshake rescales.
"""
struct LangevinBlock{F}
    f::F
    names::Vector{Symbol}
    n::Int                                  # states the integrator moves
    u_held::Vector{Float64}                 # the held proteins, appended
    p0::Vector{Float64}
    channels::Vector{Int}                   # slots of p, one per channel
    channel_names::Vector{Symbol}
    N::Matrix{Float64}                      # n × channels, integer except external rows
    ref::Vector{Int}                        # per channel, the row v_c is read on
    touches::Vector{Vector{Int}}            # per channel, the intracellular rows it moves
    enzyme_slots::Vector{Tuple{Symbol, Vector{Int}, Float64}}  # (enz slot in hybrid p, kcat slots here, E nominal)
    radius_slot::Int
end

function LangevinBlock(hybrid_models::Vector{<:AbstractSubModel}, hybrid_pnames::Vector{Symbol};
                       models = langevin_ode_models())
    prob = build_problem(models; tspan = (0.0, 1.0))
    pnames = unique(Symbol[q.name for m in models for q in model_free_params(parameters(m))])
    length(pnames) == length(prob.p) || error("parameter layout mismatch")
    names = Symbol[s for m in models for s in states(m)]
    hnames = Symbol[s for m in hybrid_models if formalism(m) === :ode for s in states(m)]
    n = length(hnames)
    names[1:n] == hnames || throw(ArgumentError(
        "The ODE-only composition must lay the hybrid's $n ODE states out first, in order"))
    slot(s) = (i = findfirst(==(s), pnames); i === nothing && throw(ArgumentError("no slot :$s")); i)

    chn = vcat(_lkcats(_LGLYC), [Symbol(k, "_", s) for s in _LPTS for k in (:kf, :kr)],
               [:p_lact2r], _lkcats(RECYCLING_REACTIONS), [:k_chg])
    ch = [slot(s) for s in chn]
    f(u, p) = prob.f(u, p, 0.0)
    u0 = collect(prob.u0)

    # Stoichiometry from the channel matrix at the initial state.
    G = _channels(f, u0, prob.p, ch)[1:n, :]
    ext = Set(i for (i, s) in enumerate(hnames) if s === :M_lac__L_e)
    N = similar(G)
    ref = zeros(Int, length(ch))
    for c in axes(G, 2)
        col = view(G, :, c)
        r = argmax([i in ext ? 0.0 : abs(col[i]) for i in eachindex(col)])
        col[r] == 0 && throw(ArgumentError("channel :$(chn[c]) has zero rate at the initial state"))
        # The coefficient on r is the smallest integer m making every
        # intracellular entry of m·col/col[r] integral.
        m = findfirst(1:6) do m
            all(i -> i in ext || abs(m * col[i] / col[r] - round(m * col[i] / col[r])) < 1e-8,
                eachindex(col))
        end
        m === nothing && throw(ArgumentError(
            "channel :$(chn[c]) has no integer stoichiometry up to 6: the parameter does not " *
            "scale one reaction direction"))
        N[:, c] = m .* col ./ col[r]
        for i in eachindex(col)
            i in ext || (N[i, c] = round(N[i, c]))
        end
        ref[c] = r
    end
    touches = [[i for i in 1:n if !(i in ext) && N[i, c] != 0] for c in axes(N, 2)]

    # The enzyme each kcat pair is scaled by, named as the hybrid's slot.
    enz = Tuple{Symbol, Vector{Int}, Float64}[]
    enom = InferCell.nominal_enzyme_concentrations()
    for (k, r) in enumerate(_LGLYC)
        push!(enz, (Symbol("enz_", r), [slot(Symbol("kcatF_", r)), slot(Symbol("kcatR_", r))], enom[k]))
    end
    rec = only(m for m in models if m isa NucleotideRecycling)
    for r in RECYCLING_REACTIONS
        e = only(q for q in parameters(rec) if q.name === Symbol("enz_", r)).value
        push!(enz, (Symbol("enz_", r), [slot(Symbol("kcatF_", r)), slot(Symbol("kcatR_", r))], e))
    end
    for (s, _, _) in enz
        s in hybrid_pnames || throw(ArgumentError("the hybrid has no enzyme slot :$s"))
    end
    return LangevinBlock(f, names, n, u0[n+1:end], collect(prob.p), ch, chn, N, ref, touches,
                         enz, slot(:r_cell_nm))
end

# g_c = p_c · ∂f/∂p_c at (u, p), for each channel slot.
function _channels(f, u, p, ch)
    U = SVector{length(u)}(u)
    J = ForwardDiff.jacobian(a -> f(U, _with(p, ch, a)), p[ch])
    return J .* p[ch]'
end
_with(p, ch, a::AbstractVector{T}) where {T} = (q = T.(p); q[ch] .= a; q)

"""
    handshake_parameters(lb, hybrid_p, hybrid_pnames) -> Vector{Float64}

The ODE-only parameter vector for one handshake: every kcat pair scaled by
the enzyme concentration the hybrid wrote, over its nominal, and the radius.
"""
function handshake_parameters(lb::LangevinBlock, hp::AbstractVector, hpn::Vector{Symbol})
    p = copy(lb.p0)
    for (s, slots, e0) in lb.enzyme_slots
        e = hp[findfirst(==(s), hpn)]
        for i in slots
            p[i] = lb.p0[i] * e / e0
        end
    end
    p[lb.radius_slot] = hp[findfirst(==(:r_cell_nm), hpn)]
    return p
end

# The full ODE-only state from the n integrated states.
_ufull(lb::LangevinBlock, u) = SVector{lb.n + length(lb.u_held)}(vcat(u, lb.u_held))

"""
    ll_step(lb, u, p, h, Ω, rng; min_particles = CONTINUUM_MEDIAN) -> (u, injected, frozen)

One local-linearization step of length `h`. `rng === nothing` takes the
deterministic step. Only channels whose every intracellular species holds at
least `min_particles` at the step's start are noisy; `frozen` counts the
others. A state still taken below zero is set to zero, and `injected` is what
that added, per state in mM (zero if nothing clamped).
"""
function ll_step(lb::LangevinBlock, u::AbstractVector, p::AbstractVector, h::Real,
                 Ω::Real, rng; min_particles::Real = CONTINUUM_MEDIAN)
    n = lb.n
    u = collect(Float64, u)
    U = _ufull(lb, u)
    f0 = lb.f(U, p)[1:n]
    J = ForwardDiff.jacobian(y -> lb.f(_ufull(lb, y), p)[1:n], u)
    # φ₁(hJ)·h·f from the exponential of [J f; 0 0]·h.
    A = zeros(n + 1, n + 1)
    A[1:n, 1:n] .= J .* h
    A[1:n, n+1] .= f0 .* h
    drift = exp(A)[1:n, n+1]
    unew = u .+ drift
    frozen = 0
    if rng !== nothing
        noisy = [all(i -> u[i] * Ω >= min_particles, lb.touches[c]) for c in eachindex(lb.touches)]
        frozen = count(!, noisy)
        G = _channels(lb.f, U, p, lb.channels)[1:n, :]
        v = [noisy[c] ? G[lb.ref[c], c] / lb.N[lb.ref[c], c] : 0.0 for c in axes(G, 2)]
        D = (lb.N .* (abs.(v) ./ Ω)') * lb.N'
        λ, V = eigen(Symmetric(noise_covariance(J, D, h)))
        unew .+= V * (sqrt.(max.(λ, 0.0)) .* randn(rng, n))
    end
    injected = max.(.-unew, 0.0)
    unew .+= injected
    return unew, injected, frozen
end

"""
    langevin_moieties(lb) -> Vector{Pair{Symbol, Vector{Float64}}}

The moieties every channel conserves, as weights on the integrated states:
14a's eight without carbon (redox, adenylate, guanylate, phosphate and the
four carriers, from [`corea_moieties`](@ref)), and the tRNA pair. A noisy
replay can move them off the reference only through what its clamps inject.
"""
function langevin_moieties(lb::LangevinBlock)
    w(ps) = (x = zeros(lb.n); for (s, v) in ps; x[findfirst(==(s), lb.names)] = v; end; x)
    ms = vcat([m.name => w(m.ode) for m in corea_moieties(carbon = false)],
              [:trna => w((:M_trna_c => 1.0, :M_trna_chg_c => 1.0))])
    for (m, x) in ms, c in axes(lb.N, 2)
        abs(x' * lb.N[:, c]) < 1e-12 || throw(ArgumentError(
            "channel :$(lb.channel_names[c]) does not conserve $m"))
    end
    return ms
end

"""
    noise_covariance(J, D, h) -> Matrix

`∫₀ʰ e^{Js} D e^{Jᵀs} ds`, by scaling and squaring (see the file header).
Every intermediate is bounded by the result, so it neither overflows nor
cancels when `J` is stiff.
"""
function noise_covariance(J::AbstractMatrix, D::AbstractMatrix, h::Real)
    m = max(0, ceil(Int, log2(max(opnorm(J, Inf) * h, eps()) / 1e-3)))
    h0 = h / 2^m
    A = J .* h0
    F = I + A + A * A / 2 + A * A * A / 6
    JD = J * D
    Σ = D .* h0 .+ (JD .+ JD') .* (h0^2 / 2) .+
        (J * JD .+ 2 .* JD * J' .+ (J * JD)') .* (h0^3 / 6)
    for _ in 1:m
        Σ = F * Σ * F' .+ Σ
        F = F * F
    end
    return Σ
end

"""
    FrozenPath

A reference run's handshakes as they acted on the ODE block: per handshake,
the dilution ratio, the ODE-only parameters `θₖ`, the particles per mM, the
unexplained part `aₖ`, and the reference state after the handshake.
"""
struct FrozenPath
    t::Vector{Float64}
    dil::Vector{Float64}
    θ::Vector{Vector{Float64}}
    Ω::Vector{Float64}
    a::Vector{Vector{Float64}}
    u_ref::Vector{Vector{Float64}}
    u0::Vector{Float64}
    dilute::Vector{Int}
    h::Float64
end

"""
    record_frozen_path(lb, hybrid_models, d, n; h = 0.05) -> FrozenPath

Run the hybrid driver `d` for `n` handshakes and record its path as it acted
on the ODE block (see the file header). `h` is the integrator step the path is
recorded for. `aₖ` depends on it, so a replay must use the same `h`.
"""
function record_frozen_path(lb::LangevinBlock, hybrid_models, d, n::Integer; h::Real = 0.05)
    hpn = unique(Symbol[q.name for m in hybrid_models if formalism(m) === :ode
                        for q in model_free_params(parameters(m))])
    length(hpn) == length(d.ode.p) || error("hybrid parameter layout mismatch")
    steps = round(Int, d.interval / h)
    steps * h ≈ d.interval || throw(ArgumentError("h must divide the handshake interval"))
    u0 = collect(Float64, d.ode.u)
    prev = copy(u0)
    fp = FrozenPath(Float64[], Float64[], Vector{Float64}[], Float64[], Vector{Float64}[],
                    Vector{Float64}[], u0, copy(d.dilute_idxs), Float64(h))
    for k in 1:n
        f_old = d.factor
        handshake_step!(d)
        dil = f_old / d.factor
        θ = handshake_parameters(lb, d.ode.p, hpn)
        u = copy(prev)
        u[fp.dilute] .*= dil
        for _ in 1:steps
            u, _ = ll_step(lb, u, θ, h, d.factor, nothing)
        end
        post = collect(Float64, d.ode.u)
        push!(fp.t, d.ode.t); push!(fp.dil, dil); push!(fp.θ, θ); push!(fp.Ω, d.factor)
        push!(fp.a, post .- u); push!(fp.u_ref, post)
        prev = post
    end
    return fp
end

"""
    replay(lb, fp, rng; noise = true, min_particles = CONTINUUM_MEDIAN) -> (traj, stats)

One trajectory along the frozen path: the ODE state after every handshake.
With `noise = false` it reproduces the reference. Adding `aₖ` can take a pool
the noise has left below the reference negative; that is clamped too, and
booked apart. `stats`:
- `clamped`, the steps whose noise clamped, and `clamped_a`, the handshakes
  whose `aₖ` did;
- `injected` and `injected_a`, per state, what those clamps added, in mM;
- `frozen`, the fraction of channel-steps the partition held noiseless;
- `drift` and `booked`, per moiety of [`langevin_moieties`](@ref), the
  largest `|w·(u − u_ref)|` over the handshakes and `w·` all the injection,
  in mM. Every channel conserves the moieties, so only a clamp moves one off
  the reference, and the drift is at most what was booked.
"""
function replay(lb::LangevinBlock, fp::FrozenPath, rng; noise::Bool = true,
                min_particles::Real = CONTINUUM_MEDIAN)
    steps = round(Int, 1.0 / fp.h)
    u = copy(fp.u0)
    traj = Vector{Vector{Float64}}(undef, length(fp.t))
    ms = langevin_moieties(lb)
    clamped = clamped_a = frozen = 0
    injected = zeros(lb.n)
    injected_a = zeros(lb.n)
    drift = zeros(length(ms))
    for k in eachindex(fp.t)
        u[fp.dilute] .*= fp.dil[k]
        for _ in 1:steps
            u, inj, fr = ll_step(lb, u, fp.θ[k], fp.h, fp.Ω[k], noise ? rng : nothing;
                                 min_particles)
            clamped += any(>(0), inj)
            injected .+= inj
            frozen += fr
        end
        u .+= fp.a[k]
        inj = max.(.-u, 0.0)
        u .+= inj
        clamped_a += any(>(0), inj)
        injected_a .+= inj
        for (j, (_, w)) in enumerate(ms)
            drift[j] = max(drift[j], abs(w' * (u .- fp.u_ref[k])))
        end
        traj[k] = copy(u)
    end
    booked = [w' * (injected .+ injected_a) for (_, w) in ms]
    return traj, (clamped = clamped, clamped_a = clamped_a, injected = injected,
                  injected_a = injected_a, frozen = frozen / (length(fp.t) * steps * length(lb.touches)),
                  moieties = first.(ms), drift = drift, booked = booked)
end

"""
    langevin_band(lb, fp, species; n = 100, seed = 1, q = (0.05, 0.95)) -> NamedTuple

An ensemble of `n` replays and, per species and handshake, the `q` quantiles
in mM. Also reported: the reference, the largest relative `excursion` of the
reference outside the band (zero if it never leaves), the time it occurs, the
steps any trajectory clamped at zero, and each replay's `stats`.
"""
function langevin_band(lb::LangevinBlock, fp::FrozenPath, species::Vector{Symbol};
                       n::Integer = 100, seed::Integer = 1, q = (0.05, 0.95),
                       reference = fp.u_ref, min_particles::Real = CONTINUUM_MEDIAN)
    idx = [findfirst(==(s), lb.names) for s in species]
    any(isnothing, idx) && throw(ArgumentError("unknown species in $species"))
    K = length(fp.t)
    vals = zeros(n, K, length(idx))
    clamped = 0
    stats = []
    for r in 1:n
        tr, st = replay(lb, fp, MersenneTwister(seed + r); min_particles)
        clamped += st.clamped
        push!(stats, st)
        for k in 1:K, (j, i) in enumerate(idx)
            vals[r, k, j] = tr[k][i]
        end
    end
    lo = [quantile(view(vals, :, k, j), q[1]) for k in 1:K, j in eachindex(idx)]
    hi = [quantile(view(vals, :, k, j), q[2]) for k in 1:K, j in eachindex(idx)]
    ref = [reference[k][i] for k in 1:K, i in idx]
    rel = @. max(lo - ref, ref - hi, 0.0) / max(abs(ref), eps())
    exc = [maximum(view(rel, :, j)) for j in eachindex(idx)]
    at = [fp.t[argmax(view(rel, :, j))] for j in eachindex(idx)]
    return (species = species, t = fp.t, lo = lo, hi = hi, ref = ref,
            excursion = exc, at = at, clamped = clamped, n = n, stats = stats)
end
