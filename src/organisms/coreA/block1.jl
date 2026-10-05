"""
Block 1 of the phase 16 sampler: the ODE constants and the metabolite noise
scale, given every cell's stochastic path (spec/phases/16-recovery.md D16.1,
task 16a.5).

The target, for the free forward constants `θF` and σ, is

```
log p(θF) + log p(σ) + Σ_c [ log p(X_c | θ) + Σ_k Σ_j log N(log y_cjk; log max(x_cjk, 1), σ) ]
```

- `x` is the carry-inclusive particle count of the replayed published model at
  each save.
- `log p(X_c | θ)` is the replay's path density ([`path_logdensity`](@ref)), which
  sees θF through the rebuilt rate constants.
- **The prior is on the forward constants only.** Each reverse constant is derived
  through the nominal equilibrium constant, as the truth was
  ([`derived_ode_values`](@ref)), and carries no prior term and no Jacobian.
- Every step is on `u = ln θ`, so each prior enters as `Normal(u; μ, s)`.

**Amended 2026-10-05 (spec §12): bulk metabolites.** With `bulk` (panel × saves),
the metabolite term is `Σ_k Σ_j log N(log z_jk; log max(x̄_jk, 1), σ_b)`, where
`x̄` is the cells' mean latent. The squared residuals and their count feed the
same σ conditional, so σ_b is updated as σ was.
"""

"""
    derived_ode_values(models, forwards) -> Vector{Pair{Symbol, Float64}}

The parameter writes for forward constants `forwards` (`name => value`), each with
its reverse constant derived through the nominal equilibrium constant. That is
the rule [`draw_truth`](@ref) draws under (spec §11 task 14c.1).
"""
function derived_ode_values(models::AbstractVector{<:AbstractSubModel}, forwards)
    nominal = nominal_parameter_values(models)
    out = Pair{Symbol, Float64}[]
    for (name, v) in forwards
        append!(out, perturbed_values(models, name, v / nominal[name]))
    end
    return out
end

"""
    Block1Cell(base, path, observed)

One cell for block 1: a snapshot of the driver at t = 0 with the stochastic
parameters written, its recorded path, and its observed metabolites (panel ×
save times). Block 1 restores the snapshot, writes θF, and replays the path.
"""
struct Block1Cell
    base::HandshakeDriver
    path::JumpPath
    observed::Matrix{Float64}
end

"""
    Block1(models, cells; forwards, panel, save_every = 60, floor = 1.0,
           sigma_prior = SIGMA_MET_PRIOR, bulk = nothing)

Block 1's fixed ingredients: the composition, the cells, the free forward
constants' names in order, the observed panel and how often it is saved, in
handshakes, and the one-particle floor of the lognormal observation. With `bulk`
(panel × saves, §12 2026-10-05) the panel is scored at the cells' mean, and each
cell's own `observed` is not read.
"""
struct Block1
    models::Vector{AbstractSubModel}
    cells::Vector{Block1Cell}
    forwards::Vector{Symbol}
    priors::Vector{Normal{Float64}}
    rows::Vector{Int}
    save_every::Int
    floor::Float64
    sigma_prior::Normal{Float64}
    bulk::Union{Nothing, Matrix{Float64}}
end

function Block1(models::AbstractVector{<:AbstractSubModel}, cells::AbstractVector{Block1Cell};
                forwards::AbstractVector{Symbol}, panel::AbstractVector{Symbol},
                save_every::Integer = 60, floor::Real = 1.0,
                sigma_prior::LogNormal = SIGMA_MET_PRIOR,
                bulk::Union{Nothing, AbstractMatrix} = nothing)
    reverses = Set(r.reverse for r in haldane_relations(models))
    bad = [f for f in forwards if f in reverses]
    isempty(bad) || throw(ArgumentError(
        "$bad are reverse constants, which block 1 derives rather than samples"))
    priors = map(forwards) do f
        p = _parameter_by_name(models, f).prior
        p isa LogNormal || throw(ArgumentError(":$f has a $(typeof(p)) prior, not LogNormal"))
        Normal(params(p)...)
    end
    rows = state_rows(collect(panel), _block_names(models, :ode))
    if bulk === nothing
        for c in cells
            size(c.observed, 1) == length(rows) || throw(DimensionMismatch(
                "a cell observes $(size(c.observed, 1)) rows for a $(length(rows))-pool panel"))
        end
    else
        size(bulk, 1) == length(rows) || throw(DimensionMismatch(
            "the bulk record has $(size(bulk, 1)) rows for a $(length(rows))-pool panel"))
    end
    return Block1(collect(models), collect(cells), collect(forwards), priors, rows,
                  Int(save_every), Float64(floor), Normal(params(sigma_prior)...),
                  bulk === nothing ? nothing : Matrix{Float64}(bulk))
end

"""
    block1_log_prior(b, u) -> Float64

The log prior of `u = ln θF`: a `Normal` per forward constant and nothing for
any reverse constant (spec/phases/16-recovery.md D16.1).
"""
block1_log_prior(b::Block1, u::AbstractVector) =
    sum(logpdf(p, x) for (p, x) in zip(b.priors, u))

"""
    block1_replay(b, u) -> (logpath, sumsq, n)

Replay every cell at `θF = exp.(u)`, with its reverse constants derived. Returns
the summed path log-density, the summed squared log residuals of the observed
panel against `max(x, floor)`, and their count. The last two are all the σ step
needs.
"""
function block1_replay(b::Block1, u::AbstractVector)
    writes = derived_ode_values(b.models, [f => exp(x) for (f, x) in zip(b.forwards, u)])
    logpath = 0.0
    sumsq = 0.0
    n = 0
    xsum = b.bulk === nothing ? nothing : zeros(size(b.bulk))
    for c in b.cells
        d = restore(c.base)
        set_parameters!(d, b.models, writes)
        r = PathReplay(c.path, d; density = true)
        nsave = b.bulk === nothing ? size(c.observed, 2) : size(b.bulk, 2)
        for k in 1:nsave
            for _ in 1:b.save_every
                replay_step!(d, r)
            end
            for (j, i) in enumerate(b.rows)
                x = d.ode.u[i] * d.factor + d.rounding.remainders[i]
                if b.bulk === nothing
                    sumsq += (log(c.observed[j, k]) - log(max(x, b.floor)))^2
                    n += 1
                else
                    xsum[j, k] += x
                end
            end
        end
        logpath += path_logdensity(r)
    end
    if b.bulk !== nothing
        C = length(b.cells)
        for k in axes(b.bulk, 2), j in axes(b.bulk, 1)
            sumsq += (log(b.bulk[j, k]) - log(max(xsum[j, k] / C, b.floor)))^2
            n += 1
        end
    end
    return logpath, sumsq, n
end

# The Gaussian log-likelihood of `n` log residuals with squared sum `sumsq`.
_log_lognormal(sumsq, n, σ) = -n * log(σ) - sumsq / (2σ^2) - n * log(2π) / 2

"""
    block1_logtarget(b, u, σ) -> Float64

Block 1's log target at `u = ln θF` and σ, up to a constant: the prior on `u`,
the path densities and the panel's lognormal likelihood.
"""
function block1_logtarget(b::Block1, u::AbstractVector, σ::Real)
    lp, ss, n = block1_replay(b, u)
    return block1_log_prior(b, u) + lp + _log_lognormal(ss, n, σ)
end

"""
    sigma_log_conditional(b, sumsq, n) -> Function

σ's conditional on `v = ln σ`, given the squared log residuals: the lognormal
likelihood plus `σ`'s prior carried to `v`. It needs no replay.
"""
sigma_log_conditional(b::Block1, sumsq::Real, n::Integer) =
    v -> _log_lognormal(sumsq, n, exp(v)) + logpdf(b.sigma_prior, v)

"""
    block1_update!(rng, b, state; w = 0.5, update_sigma = true) -> state

One block 1 sweep. It takes a slice step on each `u` coordinate in turn, each
evaluation replaying every cell, then a slice step on `ln σ` from the accepted
residuals. `update_sigma = false` holds σ, which V4 does. `state` is a mutable NamedTuple-like `(u, σ, sumsq, n)`. Every
accepted θF keeps each equilibrium constant, as `assert_haldane` confirms.
"""
function block1_update!(rng::AbstractRNG, b::Block1, state; w::Real = 0.5,
                        update_sigma::Bool = true)
    for i in eachindex(state.u)
        cache = Dict{Float64, Tuple{Float64, Float64, Int}}()
        function logf(x)
            u = copy(state.u)
            u[i] = x
            lp, ss, n = get!(() -> block1_replay(b, u), cache, x)
            return block1_log_prior(b, u) + lp + _log_lognormal(ss, n, state.σ)
        end
        x = slice_step(rng, logf, state.u[i]; w)
        state.u[i] = x
        _, state.sumsq, state.n = cache[x]
    end
    assert_haldane(b.models, derived_ode_values(b.models,
                                                [f => exp(x) for (f, x) in zip(b.forwards, state.u)]))
    update_sigma &&
        (state.σ = exp(slice_step(rng, sigma_log_conditional(b, state.sumsq, state.n),
                                  log(state.σ))))
    return state
end

"""
    Block1State(u, σ, sumsq, n)

The chain's block 1 state: `u = ln θF`, σ, and the accepted replay's squared log
residuals and their count.
"""
mutable struct Block1State
    u::Vector{Float64}
    σ::Float64
    sumsq::Float64
    n::Int
end

"""
    Block1State(b, u, σ) -> Block1State

A state at `u` and σ, with the residuals from one replay.
"""
function Block1State(b::Block1, u::AbstractVector, σ::Real)
    _, ss, n = block1_replay(b, u)
    return Block1State(collect(Float64, u), Float64(σ), ss, n)
end

export derived_ode_values, Block1Cell, Block1, Block1State, block1_log_prior,
       block1_replay, block1_logtarget, sigma_log_conditional, block1_update!
