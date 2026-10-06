"""
Block 3 of the phase 16 sampler: one cell's stochastic path given θ and its
data, by conditional SMC over its 60 s observation windows
(spec/phases/16-recovery.md D16.3, task 16a.7).

Window `m` covers `(60(m − 1), 60m]`. Every particle's proposal for a window is
built outside the SSA, and the hybrid is advanced by replaying it, which is
checked bitwise (V2). So a particle is a driver snapshot, and its randomness is
its own generator, seeded from the sweep's.

**The proposal.** On the jump clock the rate constants change at 60m − 1 s, since
`handshake_step!` rebuilds before the jump step.
1. Over the window's first 59 s, at the window's starting constants `R_{m−1}`,
   each gene's transcript path is the first 59 s of a 60 s bridge from its
   observed count `y_{m−1}` to `y_m`. Every other reaction, translation and
   translocation here, is simulated exactly given those transcripts, by a
   Gillespie on a working copy of the jump state with the composition's own
   propensities.
2. The hybrid replays those 59 s, then runs step 60m's exchange, which rebuilds
   `R_m` from this particle's pools.
3. Over the last second, each gene is bridged from its count at 60m − 1 to `y_m`
   at `R_m`, the other reactions are simulated again, and the events are applied.

**The incremental weight** is target over proposal. The translation path is
proposed from its own conditional, so it cancels, and

```
w_m = g_m(y_met) · Π_g P_{R_{m−1}}(y_m | y_{m−1}; 60) · P_{R_m}(y_m | c_g; 1) / P_{R_{m−1}}(y_m | c_g; 1)
```

where `c_g` is gene g's count at 60m − 1 and `g_m` is window m's metabolite
likelihood at 60m. The transition probabilities differ between particles,
because each particle's constants are rebuilt from its own pools
(spec/phases/16-recovery.md §12, the #80 review).
"""

"""
    TranscriptMap(d, states, births, deaths)

Which jump states are the observed transcripts, and which composed reactions
make and remove each one, by index into the driver's jump state and its event
log. [`transcript_map`](@ref) builds Core A′'s and M0's.
"""
struct TranscriptMap
    states::Vector{Int}
    births::Vector{Int}
    deaths::Vector{Int}
end

"""
    transcript_map(d, models) -> TranscriptMap

The transcripts of a Core A′ composition (Core A′ or M0), in transcription's gene
order, with transcription's and decay's reactions for each.
"""
function transcript_map(d::HandshakeDriver, models::AbstractVector{<:AbstractSubModel})
    tx = only(m for m in models if m isa CoreATranscription)
    dec = only(m for m in models if m isa CoreATranscriptDecay)
    names = _block_names(models, :jump)
    labels = d.events.labels
    states = [findfirst(==(transcript_state(g.locus)), names) for g in tx.genes]
    births = [findfirst(==(Symbol(:CoreATranscription_, i)), labels) for i in eachindex(tx.genes)]
    jdec = [findfirst(h -> h.locus === g.locus, dec.genes) for g in tx.genes]
    deaths = [findfirst(==(Symbol(:CoreATranscriptDecay_, j)), labels) for j in jdec]
    any(isnothing, vcat(states, births, deaths)) && error("a transcript's state or reaction is missing")
    return TranscriptMap(states, births, deaths)
end

"""
    CellData(transcripts, metabolites)

One cell's observations by window: `transcripts` is genes × windows, the exact
counts at 60m for m = 1..T, in the [`TranscriptMap`](@ref)'s order, and
`metabolites` is panel × windows, the observed values at the same saves.
"""
struct CellData
    transcripts::Matrix{Int}
    metabolites::Matrix{Float64}
end

"""
    CSMCSpec(tmap, panel_rows; sigma, floor = 1.0, window = 60, cap_floor)

What a sweep holds fixed: the transcript map, the ODE rows of the observed panel,
the metabolite noise σ and its floor, the window in handshakes, and the floor
under every bridge's count cap (`bridge_cap` of the data's largest count).
"""
struct CSMCSpec
    tmap::TranscriptMap
    panel_rows::Vector{Int}
    sigma::Float64
    floor::Float64
    window::Int
    cap_floor::Int
end
CSMCSpec(tmap::TranscriptMap, panel_rows::AbstractVector{Int}; sigma::Real, floor::Real = 1.0,
         window::Integer = 60, cap_floor::Integer) =
    CSMCSpec(tmap, collect(panel_rows), Float64(sigma), Float64(floor), Int(window), Int(cap_floor))

"""
    BulkObservation(z, others, ncells)

One save's bulk metabolite term, as one cell's update sees it
(spec/phases/16-recovery.md §3 and D16.3, §12 2026-10-05). `z` is the observed bulk
value per panel pool, `others` the other cells' summed latent at that save, and
`ncells` the population size C. A particle with latent `x` is scored at the
population mean `(others + x)/C`.
"""
struct BulkObservation
    z::Vector{Float64}
    others::Vector{Float64}
    ncells::Int
end

"""
    bulk_window_observations(z, others, ncells) -> Vector{BulkObservation}

One [`BulkObservation`](@ref) per window from `z` and `others`, both panel × windows.
"""
bulk_window_observations(z::AbstractMatrix, others::AbstractMatrix, ncells::Integer) =
    [BulkObservation(z[:, m], others[:, m], Int(ncells)) for m in axes(z, 2)]

"""
    WindowEvents(times, reactions)

One window's events for one particle, in order, by composed reaction index.
"""
struct WindowEvents
    times::Vector{Float64}
    reactions::Vector{Int}
end

# A working jump state that the composed affects can mutate.
mutable struct _WorkU{U}
    u::U
end

# Gene g's birth rate and per-transcript death rate at the driver's live state.
# The death propensity is first order, so its rate at a count of one is μ.
function _birth_death(d::HandshakeDriver, tm::TranscriptMap, g::Int, u = d.jump.u)
    rates = d.events.rates
    k = rates[tm.births[g]](u, d.jump.p, d.jump.t)
    u1 = copy(u)
    u1[tm.states[g]] = 1
    μ = rates[tm.deaths[g]](u1, d.jump.p, d.jump.t)
    return k, μ
end

_chain(k, μ, τ, a, cap_floor) = TranscriptChain(k, μ, adequate_cap(k, μ, τ, a; floor = cap_floor))

"""
    propose_events(rng, d, spec, imposed, t0, t1; u0 = d.jump.u) -> WindowEvents

Every event over `[t0, t1)` given the imposed transcript events: those are
applied at their times, and every other reaction is simulated exactly from the
jump state `u0` (by default the driver's live state) and the driver's constants,
with the composition's own propensities, on a working copy. The driver is not
changed.
"""
function propose_events(rng::AbstractRNG, d::HandshakeDriver, spec::CSMCSpec,
                        imposed::WindowEvents, t0::Real, t1::Real; u0 = d.jump.u)
    tm = spec.tmap
    fixed = Set(vcat(tm.births, tm.deaths))
    free = [j for j in eachindex(d.events.rates) if !(j in fixed)]
    rates, affects, p = d.events.rates, d.events.affects, d.jump.p
    w = _WorkU(copy(u0))
    times, rxns = Float64[], Int[]
    t = Float64(t0)
    k = 1
    λ = zeros(length(free))
    while true
        for (i, j) in enumerate(free)
            λ[i] = rates[j](w.u, p, t)
        end
        total = sum(λ)
        tnext = total > 0 ? t + randexp(rng) / total : Inf
        timp = k <= length(imposed.times) ? imposed.times[k] : Inf
        if timp <= tnext && timp < t1
            affects[imposed.reactions[k]](w)
            push!(times, timp)
            push!(rxns, imposed.reactions[k])
            t = timp
            k += 1
        elseif tnext < t1
            c = rand(rng) * total
            i = findfirst(>=(c), cumsum(λ))
            # Rounding can leave c above the cumulative sum: take the last
            # reaction that can fire, never one with zero propensity.
            i === nothing && (i = findlast(>(0), λ))
            affects[free[i]](w)
            push!(times, tnext)
            push!(rxns, free[i])
            t = tnext
        else
            break
        end
    end
    return WindowEvents(times, rxns)
end

# Bridge every gene over `τ` from its count in `u` to `yend`, keeping events
# before `keep`, as absolute times from `t0`. Returns the imposed events and
# each gene's count at `t0 + keep`.
function _bridge_genes(rng, d, spec, u, yend, t0, τ, keep)
    tm = spec.tmap
    ts, rs = Float64[], Int[]
    counts = Int[]
    for g in eachindex(tm.states)
        a = u[tm.states[g]]
        k, μ = _birth_death(d, tm, g, u)
        times, deltas = sample_bridge(rng, _chain(k, μ, τ, a, spec.cap_floor), τ, a, yend[g])
        c = a
        for (x, δ) in zip(times, deltas)
            x < keep || break
            push!(ts, t0 + x)
            push!(rs, δ > 0 ? tm.births[g] : tm.deaths[g])
            c += δ
        end
        push!(counts, c)
    end
    o = sortperm(ts)
    return WindowEvents(ts[o], rs[o]), counts
end

# The log transition probability of every gene from `ystart` to `yend` over τ
# at the live constants, and each chain, at the counts in `u`.
function _log_transitions(d, spec, u, ystart, yend, τ)
    tm = spec.tmap
    s = 0.0
    for g in eachindex(tm.states)
        k, μ = _birth_death(d, tm, g, u)
        s += log(transition_probability(_chain(k, μ, τ, ystart[g], spec.cap_floor), τ,
                                        ystart[g], yend[g]))
    end
    return s
end

"""
    window_log_weight(lp60_prev, lp1_new, lp1_prev, loglik_met) -> Float64

The incremental log weight: `log P_{R_{m−1}}(y_m | y_{m−1}; 60) + log P_{R_m}(y_m | c; 1)
− log P_{R_{m−1}}(y_m | c; 1)`, summed over genes, plus the metabolite
log-likelihood.
"""
window_log_weight(lp60_prev, lp1_new, lp1_prev, loglik_met) =
    lp60_prev + lp1_new - lp1_prev + loglik_met

_panel_latent(d::HandshakeDriver, spec::CSMCSpec) =
    [d.ode.u[i] * d.factor + d.rounding.remainders[i] for i in spec.panel_rows]

# The bulk log-likelihood with this particle's latent in the population mean.
function _metabolite_loglik(d::HandshakeDriver, spec::CSMCSpec, y::BulkObservation)
    s = 0.0
    for (j, x) in enumerate(_panel_latent(d, spec))
        xbar = (y.others[j] + x) / y.ncells
        s += logpdf(Normal(log(max(xbar, spec.floor)), spec.sigma), log(y.z[j]))
    end
    return s
end

# The panel's lognormal log-likelihood at the driver's live state.
function _metabolite_loglik(d::HandshakeDriver, spec::CSMCSpec, y::AbstractVector)
    s = 0.0
    for (j, i) in enumerate(spec.panel_rows)
        x = d.ode.u[i] * d.factor + d.rounding.remainders[i]
        s += logpdf(Normal(log(max(x, spec.floor)), spec.sigma), log(y[j]))
    end
    return s
end

# Apply `events` from the driver's current time through the end of the
# handshake that contains `t_end`: whole handshakes by replay, with the events
# each one owns.
function _replay_events!(d::HandshakeDriver, events::WindowEvents, nsteps::Integer)
    path = JumpPath(events.times, events.reactions, d.events.labels)
    r = PathReplay(path, d)
    for _ in 1:nsteps
        replay_step!(d, r)
    end
    r.next == length(events.times) + 1 || error("events beyond the replayed handshakes")
    return d
end

# Step 60m: the exchange, which rebuilds R_m, then the jump interval with
# `events`. The tail of `replay_step!`, with the events supplied after the
# exchange because they are proposed at R_m.
function _last_step!(d::HandshakeDriver, events_after_exchange)
    _, drained, clipped = _exchange!(d)
    target = d.jump.t + d.interval
    ev = events_after_exchange(d)
    affects = d.events.affects
    for (x, k) in zip(ev.times, ev.reactions)
        d.jump.t <= x < target || error("an event at $x lies outside [$(d.jump.t), $target)")
        affects[k](d.jump)
    end
    d.jump.tprev = d.jump.t
    d.jump.t = target
    _close_handshake!(d, drained, clipped)
    return ev
end

"""
    advance_window!(rng, d, spec, y_prev, y_end, ymet; fixed = nothing, head = nothing)
        -> (log_weight, events)

Advance one particle across one window and return its incremental log weight and
the window's events. With `fixed` (a reference particle's [`WindowEvents`](@ref))
the events are replayed rather than proposed, and the same weight is computed for
them. With `head = (events, s)`, the events before `s` are kept and the rest of
the window is drawn from the proposal's own conditional given them: each gene's
bridge from its count at `s` to `y_end`, the other reactions forward from the
state at `s`, then the last second as usual. `s` must lie in the window's first
59 s. This is [`annealed_window`](@ref)'s move. The driver ends at 60m, and its
transcripts are checked against `y_end`.
"""
function advance_window!(rng::AbstractRNG, d::HandshakeDriver, spec::CSMCSpec,
                         y_prev::AbstractVector{<:Integer}, y_end::AbstractVector{<:Integer},
                         ymet; fixed::Union{Nothing, WindowEvents} = nothing,
                         head::Union{Nothing, Tuple{WindowEvents, Float64}} = nothing)
    fixed === nothing || head === nothing ||
        throw(ArgumentError("a window is either replayed or proposed from a head, not both"))
    tm = spec.tmap
    W = spec.window
    t0 = d.jump.t
    tlast = t0 + (W - 1) * d.interval            # 60m − 1
    u0 = copy(d.jump.u)
    [u0[s] for s in tm.states] == y_prev || error("the particle's transcripts are not y_{m−1}")
    lp60 = _log_transitions(d, spec, u0, y_prev, y_end, W * d.interval)
    # A proposal cannot reach y_end at these constants: its target density is
    # zero, so it gets weight zero and is not advanced. It is never resampled.
    fixed === nothing && lp60 == -Inf && return -Inf, WindowEvents(Float64[], Int[])

    if fixed === nothing && head === nothing
        imposed, _ = _bridge_genes(rng, d, spec, u0, y_end, t0, W * d.interval,
                                   (W - 1) * d.interval)
        first_part = propose_events(rng, d, spec, imposed, t0, tlast)
    elseif fixed === nothing
        hev, s = head
        t0 <= s < tlast || throw(ArgumentError("the head's cut $s lies outside [$t0, $tlast)"))
        k = searchsortedfirst(hev.times, s)
        # The jump state at s: the window's start with the head's events applied.
        w = _WorkU(copy(u0))
        for j in view(hev.reactions, 1:(k - 1))
            d.events.affects[j](w)
        end
        imposed, _ = _bridge_genes(rng, d, spec, w.u, y_end, s, t0 + W * d.interval - s, tlast - s)
        tail = propose_events(rng, d, spec, imposed, s, tlast; u0 = w.u)
        first_part = WindowEvents(vcat(hev.times[1:(k - 1)], tail.times),
                                  vcat(hev.reactions[1:(k - 1)], tail.reactions))
    else
        cut = searchsortedfirst(fixed.times, tlast)
        first_part = WindowEvents(fixed.times[1:cut - 1], fixed.reactions[1:cut - 1])
    end
    _replay_events!(d, first_part, W - 1)
    c = [d.jump.u[s] for s in tm.states]
    lp1_prev = _log_transitions(d, spec, d.jump.u, c, y_end, d.interval)  # still R_{m−1}

    lp1_new = Ref(0.0)
    last = _last_step!(d, e -> begin
        lp1_new[] = _log_transitions(e, spec, e.jump.u, c, y_end, e.interval)   # R_m
        if fixed === nothing && lp1_new[] == -Inf
            WindowEvents(Float64[], Int[])      # weight zero, as above
        elseif fixed === nothing
            imp, _ = _bridge_genes(rng, e, spec, e.jump.u, y_end, tlast, e.interval, e.interval)
            propose_events(rng, e, spec, imp, tlast, tlast + e.interval)
        else
            cut = searchsortedfirst(fixed.times, tlast)
            WindowEvents(fixed.times[cut:end], fixed.reactions[cut:end])
        end
    end)
    lp1_new[] == -Inf && fixed === nothing && return -Inf, WindowEvents(Float64[], Int[])
    [d.jump.u[s] for s in tm.states] == y_end || error(
        "the window ended at transcripts $([d.jump.u[s] for s in tm.states]), not $y_end")
    lw = window_log_weight(lp60, lp1_new[], lp1_prev, _metabolite_loglik(d, spec, ymet))
    return lw, WindowEvents(vcat(first_part.times, last.times),
                            vcat(first_part.reactions, last.reactions))
end

_normalise(lw) = (w = exp.(lw .- maximum(lw)); w ./ sum(w))

function _categorical(rng::AbstractRNG, w)
    u = rand(rng)
    c = 0.0
    for i in eachindex(w)
        c += w[i]
        u < c && return i
    end
    # Rounding can leave u above the sum: the last particle with weight, never
    # one with weight zero.
    return findlast(>(0), w)
end

# The log density of `ref`'s windows m..m+L−1 given a particle's state at the
# start of window m: path density plus metabolite likelihood, by replay.
function _future_logdensity(d0::HandshakeDriver, spec::CSMCSpec, data::CellData,
                            ref::Vector{WindowEvents}, m::Int, L::Int)
    d = restore(d0)
    T = size(data.transcripts, 2)
    s = 0.0
    for j in m:min(T, m + L - 1)
        path = JumpPath(ref[j].times, ref[j].reactions, d.events.labels)
        r = PathReplay(path, d; density = true)
        for _ in 1:spec.window
            replay_step!(d, r)
        end
        s += path_logdensity(r) + _metabolite_loglik(d, spec, view(data.metabolites, :, j))
    end
    return s
end

"""
    geometric_schedule(K; β1 = 1e-4) -> Vector{Float64}

`K` inverse temperatures rising geometrically from `β1` to exactly 1, for
[`annealed_window`](@ref) (spec/phases/16-recovery.md D16.3, §12 2026-10-02).
"""
function geometric_schedule(K::Integer; β1::Real = 1e-4)
    K >= 1 || throw(ArgumentError("a schedule needs at least one stage"))
    0 < β1 <= 1 || throw(ArgumentError("β1 must lie in (0, 1]"))
    K == 1 && return [1.0]
    β = [β1 * (1 / β1)^((k - 1) / (K - 1)) for k in 1:K]
    β[end] = 1.0
    return β
end

function _check_schedule(schedule)
    isempty(schedule) && throw(ArgumentError("a schedule needs at least one stage"))
    first(schedule) > 0 && issorted(schedule; lt = <=) && last(schedule) == 1 ||
        throw(ArgumentError("a schedule must rise strictly from above 0 to exactly 1"))
    return nothing
end

# One stage move at β, from the window-start driver `d`: keep `x` before a cut s
# drawn uniformly in the window's first 59 s, and redraw the rest from the
# proposal's conditional. The proposal cancels against the target π_β ∝ q · w^β,
# so the acceptance is min(1, (w′/w)^β), and the move is π_β-reversible. Returns
# the path, its log weight, and its driver at the window's end, or `nothing` for
# the driver when the move is rejected.
function _stage_move(rng, d, spec, y_prev, y_end, ymet, x::WindowEvents, lwx, β)
    s = d.jump.t + rand(rng) * (spec.window - 1) * d.interval
    e = restore(d)
    lw, ev = advance_window!(rng, e, spec, y_prev, y_end, ymet; head = (x, s))
    isfinite(lw) && log(rand(rng)) < β * (lw - lwx) && return ev, lw, e
    return x, lwx, nothing
end

"""
    annealed_window(rng, d, spec, y_prev, y_end, ymet, schedule; fixed = nothing)
        -> (log_weight, events, driver)

One particle's window by annealed importance sampling (spec/phases/16-recovery.md
D16.3, §12 2026-10-02). `d` is the particle at the window's start, and it is not
changed. The tempered targets are π_β ∝ q · w^β, where q is
[`advance_window!`](@ref)'s proposal and w its incremental weight, at each β of
`schedule`, which rises to 1. Each stage is one move: keep the path before a
uniform cut and redraw the rest from q's conditional, accepted at min(1, (w′/w)^β).

- **A fresh particle** draws x₀ from q and moves it through the stages. Its log
  weight is the AIS log weight, Σₖ (βₖ − βₖ₋₁) log w(xₖ₋₁), whose exponential is
  an unbiased estimate of ∫ q·w. Its events are the chain's final state.
- **The reference** (`fixed`) keeps its path. Its auxiliary chain is drawn
  backwards, xₖ₋₁ by the stage-k move from xₖ, for k = K down to 1, which is the
  extended target's conditional because each move is its own reversal. Its weight
  is computed the same way.

Returns the log weight, the window's events, and the particle's driver at 60m. A
proposal that cannot reach `y_end` has log weight −Inf.
"""
function annealed_window(rng::AbstractRNG, d::HandshakeDriver, spec::CSMCSpec,
                         y_prev::AbstractVector{<:Integer}, y_end::AbstractVector{<:Integer},
                         ymet, schedule::AbstractVector{<:Real};
                         fixed::Union{Nothing, WindowEvents} = nothing)
    _check_schedule(schedule)
    K = length(schedule)
    if fixed === nothing
        final = restore(d)
        lwx, x = advance_window!(rng, final, spec, y_prev, y_end, ymet)
        isfinite(lwx) || return -Inf, x, final
        logw = 0.0
        βprev = 0.0
        for β in schedule
            logw += (β - βprev) * lwx
            x, lwx, e = _stage_move(rng, d, spec, y_prev, y_end, ymet, x, lwx, β)
            e === nothing || (final = e)
            βprev = β
        end
        return logw, x, final
    end
    final = restore(d)
    lwx, _ = advance_window!(rng, final, spec, y_prev, y_end, ymet; fixed)
    x = fixed
    lws = zeros(K)                       # lws[k] = log w(x_{k−1})
    for k in K:-1:1
        x, lwx, _ = _stage_move(rng, d, spec, y_prev, y_end, ymet, x, lwx, schedule[k])
        lws[k] = lwx
    end
    logw = sum((schedule[k] - (k == 1 ? 0.0 : schedule[k - 1])) * lws[k] for k in 1:K)
    return logw, fixed, final
end

"""
    csmc_sweep(rng, base, spec, data, ref; N, lag = 0, schedule = nothing, bulk = nothing)
        -> (path, changed, latent)

One conditional SMC sweep for one cell (spec/phases/16-recovery.md D16.3).
`base` is the cell's driver at t = 0 at the current θ, `ref` the current path as
one [`WindowEvents`](@ref) per window, and `N` the particle count, the reference
included. Returns the new path and, per window, whether its events changed.

- `schedule` (a rising vector ending at 1, as [`geometric_schedule`](@ref) gives)
  anneals each particle's window by [`annealed_window`](@ref). It needs `lag = 0`.
- `lag = 0` is particle Gibbs: the reference keeps its own ancestor.
- `lag = L > 0` is ancestor sampling truncated at `L` windows. The reference's
  ancestor at window m is drawn in proportion to each candidate's weight times
  the density of the reference's next `L` windows from that candidate's state.
  At `L ≥ T` that is full PGAS, which is exact. Below it, the dropped windows
  bias the kernel (D16.5).

Each particle draws from its own generator, seeded from `rng`, so a sweep is
reproducible from one seed and uses no global stream.

`bulk` (one [`BulkObservation`](@ref) per window) replaces the per-cell panel in
`data.metabolites` with the bulk likelihood at the population mean, this cell's
particle substituted (§12 2026-10-05). `latent` is the returned path's panel
latent at each window's end, panel × windows, so a caller updating cells in turn
can refresh the others' sums without a replay.
"""
function csmc_sweep(rng::AbstractRNG, base::HandshakeDriver, spec::CSMCSpec, data::CellData,
                    ref::Vector{WindowEvents}; N::Integer, lag::Integer = 0,
                    schedule::Union{Nothing, AbstractVector{<:Real}} = nothing,
                    bulk::Union{Nothing, AbstractVector{BulkObservation}} = nothing)
    T = size(data.transcripts, 2)
    bulk === nothing || length(bulk) == T ||
        throw(ArgumentError("$(length(bulk)) bulk observations for $T windows"))
    lag > 0 && bulk !== nothing &&
        throw(ArgumentError("ancestor sampling's future density does not score bulk data"))
    obs(m) = bulk === nothing ? view(data.metabolites, :, m) : bulk[m]
    if schedule !== nothing
        lag == 0 || throw(ArgumentError("an annealed window runs under particle Gibbs only (lag = 0)"))
        _check_schedule(schedule)
    end
    length(ref) == T || throw(ArgumentError("the reference has $(length(ref)) windows, not $T"))
    N >= 2 || throw(ArgumentError("need at least two particles"))
    y0 = [base.jump.u[s] for s in spec.tmap.states]
    particles = [restore(base) for _ in 1:N]
    events = Matrix{WindowEvents}(undef, N, T)
    ancestors = zeros(Int, N, T)
    lw = zeros(N)
    lat = Matrix{Vector{Float64}}(undef, N, T)
    for m in 1:T
        yprev = m == 1 ? y0 : data.transcripts[:, m - 1]
        if m > 1
            w = _normalise(lw)
            a = [_categorical(rng, w) for _ in 1:N]
            if lag > 0
                la = [log(w[i]) + _future_logdensity(particles[i], spec, data, ref, m, lag)
                      for i in 1:N]
                a[1] = _categorical(rng, _normalise(la))
            else
                a[1] = 1
            end
            particles = [restore(particles[a[i]]) for i in 1:N]
            ancestors[:, m] = a
        else
            ancestors[:, 1] = 1:N
        end
        seeds = rand(rng, UInt64, N)
        for i in 1:N
            prng = Xoshiro(seeds[i])
            fixed = i == 1 ? ref[m] : nothing
            if schedule === nothing
                lw[i], events[i, m] = advance_window!(prng, particles[i], spec, yprev,
                                                      data.transcripts[:, m], obs(m); fixed)
            else
                lw[i], events[i, m], particles[i] =
                    annealed_window(prng, particles[i], spec, yprev, data.transcripts[:, m],
                                    obs(m), schedule; fixed)
            end
            lat[i, m] = _panel_latent(particles[i], spec)
        end
    end
    k = _categorical(rng, _normalise(lw))
    path = Vector{WindowEvents}(undef, T)
    latent = zeros(length(spec.panel_rows), T)
    for m in T:-1:1
        path[m] = events[k, m]
        latent[:, m] = lat[k, m]
        k = ancestors[k, m]
    end
    changed = [path[m].times != ref[m].times || path[m].reactions != ref[m].reactions for m in 1:T]
    return path, changed, latent
end

"""
    window_events(path, T; window = 60.0) -> Vector{WindowEvents}

A recorded [`JumpPath`](@ref) split into its `T` windows, `[60(m − 1), 60m)`.
"""
function window_events(path::JumpPath, T::Integer; window::Real = 60.0)
    out = WindowEvents[]
    for m in 1:T
        idx = findall(t -> (m - 1) * window <= t < m * window, path.times)
        push!(out, WindowEvents(path.times[idx], path.reactions[idx]))
    end
    return out
end

"The windows' events joined back into one [`JumpPath`](@ref)."
join_path(ws::AbstractVector{WindowEvents}, labels) =
    JumpPath(reduce(vcat, [w.times for w in ws]), reduce(vcat, [w.reactions for w in ws]), labels)

export TranscriptMap, transcript_map, CellData, CSMCSpec, WindowEvents, propose_events,
       window_log_weight, advance_window!, csmc_sweep, window_events, join_path,
       geometric_schedule, annealed_window, BulkObservation, bulk_window_observations
