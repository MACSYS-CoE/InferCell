"""
Track B's approximate likelihood, and NUTS on it (spec/phases/approximate-method.md
§3, task B5; §12 2026-10-23, "B5 as planned").

For θ, with σ_b sampled beside it,

```
p̃(θ | D) ∝ p(θ) · Π_cells Π_g L_tx,g(θ; k̄_tx(θ)) · N(log z; m(θ) + δ, σ_b² I + S)
```

- `m(θ)` is `log max(x̄, 1)` on the panel at the 60 s saves, with x̄ in particles,
  from route (ii) ([`build_mean_field_problem`](@ref)).
- `k̄_tx(θ)` are the transcription constants route (ii) rebuilt, read at t = 0 and
  at each switch, `60m − 1`, from the route's own pools. They are scored by B1's
  exact factor ([`transcript_loglik`](@ref)).
- δ and S come from a reference ensemble at a reference θ and are held fixed
  ([`reference_ensemble`](@ref), [`delta_and_S`](@ref)). At δ = 0 and S = 0 the
  bulk factor is block 1's lognormal.

The density is over x = (log θ, log σ_b), with the exact sampler's priors: a
`Normal` on each log, and each reverse constant derived through Haldane with no
prior term. A θ the route or the factor cannot score has density −∞, not an error.
"""

"""
    MixedSetup

What a fit of the mixed likelihood holds fixed: the original composition, its route
(ii) build, one [`TranscriptFactor`](@ref) per cell, where each gene's rebuilt
constant sits, the panel rows, the save and switch times, the log bulk record, δ
and S, and the parameters with their priors on the log scale. Built by
[`mixed_setup`](@ref).
"""
struct MixedSetup{R <: MeanFieldRoute}
    models::Vector{AbstractSubModel}
    route::R
    factors::Vector{TranscriptFactor}
    loci::Vector{Symbol}
    lengths::Vector{Int}
    tx_slots::Vector{Tuple{Int, Int}}
    rows::Vector{Int}
    saves::Vector{Float64}
    switches::Vector{Float64}
    logz::Vector{Float64}
    delta::Vector{Float64}
    S::Matrix{Float64}
    floor::Float64
    names::Vector{Symbol}
    forwards::Vector{Symbol}
    prior_mu::Vector{Float64}
    prior_sd::Vector{Float64}
    sigma_prior::Normal{Float64}
    abstol::Float64
    reltol::Float64
    maxiters::Int
    box_sd::Float64
end

"""
    mixed_setup(ds; models = m0_models(), targets = M0_TARGETS, panel = M0_PANEL,
                delta = nothing, S = nothing, floor = 1.0,
                sigma_prior = SIGMA_MET_PRIOR, window = 60, abstol = 1e-10,
                reltol = 1e-8, maxiters = 10^4, box_sd = 6) -> MixedSetup

A fit's fixed ingredients for a dataset `ds` with a bulk record, as
[`m0_dataset`](@ref) returns. `delta` (panel × saves, or its `vec`) and `S`
(square, in the same order, pools fastest) default to zero, which is pass 1 of
the reference θ (§12 2026-10-23, "B5 as planned").

Each cell's transcripts start from the composition's fixed initial condition,
since the record drops t = 0. `abstol` and `reltol` are route (ii)'s solve
tolerances, [`solve_mean_field`](@ref)'s by default. `maxiters` caps its steps: an
ordinary M0 solve takes about 1,200 out to ±4 prior SD. x outside the box
`prior mean ± box_sd · prior SD`, per coordinate, has density −∞ and is not solved
(§12 2026-10-23, "B5's support box"). Each gene's factor is sized for the support:
`kmax = TURNOVER_CEILING / (len − 1)`, which bounds the transcription constant, and
`μmin` = krnadeg's prior quantile at 1e-6 over `len`.
"""
function mixed_setup(ds; models = m0_models(), targets = M0_TARGETS,
                     panel::AbstractVector{Symbol} = M0_PANEL, delta = nothing, S = nothing,
                     floor::Real = 1.0, sigma_prior::LogNormal = SIGMA_MET_PRIOR,
                     window::Real = 60, abstol::Real = 1e-10, reltol::Real = 1e-8,
                     maxiters::Integer = 10^4, box_sd::Real = 6)
    ms = collect(AbstractSubModel, models)
    times = Float64.(ds.times)
    M = length(times)
    times == window .* (1:M) || throw(ArgumentError(
        "the record's saves must be every $window s from $window s, as the rebuild's are"))
    rt = build_mean_field_problem(ms; tspan = (0.0, times[end]), complete = true)
    rt.rebuild_times == times .- rt.interval || throw(ArgumentError(
        "route (ii)'s switches $(rt.rebuild_times) are not one handshake before the saves"))

    tx = only(m for m in ms if m isa CoreATranscription)
    loci = [g.locus for g in tx.genes]
    lens = [g.length for g in tx.genes]
    dec = only(m for m in ms if m isa CoreATranscriptDecay)
    [g.length for g in dec.genes] == lens || error(
        "transcription and decay disagree on the genes or their lengths")
    tx_slots = map(loci) do l
        name = Symbol(:k_tx_, l)
        hits = [(i, j) for (i, r) in enumerate(rt.rebuilds) for (j, n) in enumerate(r.names)
                if n === name]
        length(hits) == 1 || error("route (ii) rebuilds :$name $(length(hits)) times, not once")
        only(hits)
    end

    gi = [findfirst(==(l), ds.genes) for l in loci]
    any(isnothing, gi) && error("the dataset does not record every transcribed gene")
    y0 = map(loci) do l
        v = rt.prob.u0[findfirst(==(transcript_state(l)), rt.state_names)]
        isinteger(v) || error("the initial transcript count of $l is $v, not whole")
        Int(v)
    end
    krnadeg_prior = _parameter_by_name(ms, :krnadeg).prior
    kmax = TURNOVER_CEILING ./ (lens .- 1)
    μmin = quantile(krnadeg_prior, 1e-6) ./ lens
    factors = [TranscriptFactor(y0, Matrix{Int}(ds.transcripts[c, gi, :]); kmax, μmin)
               for c in axes(ds.transcripts, 1)]

    rows = [findfirst(==(o), rt.state_names) for o in panel]
    any(isnothing, rows) && error("route (ii) does not hold every panel pool")
    size(ds.bulk) == (length(panel), M) || throw(DimensionMismatch(
        "the bulk record is $(size(ds.bulk)), not $(length(panel)) pools × $M saves"))
    n = length(panel) * M
    δ = delta === nothing ? zeros(n) : Float64.(vec(delta))
    Σ = S === nothing ? zeros(n, n) : Matrix{Float64}(S)
    length(δ) == n && size(Σ) == (n, n) || throw(DimensionMismatch(
        "δ and S must have $n entries per side, pools fastest"))

    priors = map(n -> _parameter_by_name(ms, n).prior, targets)
    all(p -> p isa LogNormal, priors) || throw(ArgumentError("every target needs a LogNormal prior"))
    forwards = [r.forward for r in haldane_relations(ms) if r.forward in targets]
    return MixedSetup(ms, rt, factors, loci, lens, tx_slots, rows, times,
                      vcat(0.0, rt.rebuild_times), vec(log.(ds.bulk)), δ, Σ, Float64(floor),
                      collect(Symbol, targets), forwards,
                      [Float64(params(p)[1]) for p in priors],
                      [Float64(params(p)[2]) for p in priors],
                      Normal(Float64.(params(sigma_prior))...), Float64(abstol), Float64(reltol),
                      Int(maxiters), Float64(box_sd))
end

"""
    with_reference(s, delta, S) -> MixedSetup

`s` with δ and S replaced: pass 2 of the reference θ.
"""
function with_reference(s::MixedSetup, delta, S)
    δ, Σ = Float64.(vec(delta)), Matrix{Float64}(S)
    length(δ) == length(s.delta) && size(Σ) == size(s.S) ||
        throw(DimensionMismatch("δ or S does not match the setup's panel × saves"))
    return MixedSetup(s.models, s.route, s.factors, s.loci, s.lengths, s.tx_slots, s.rows,
                      s.saves, s.switches, s.logz, δ, Σ, s.floor, s.names, s.forwards,
                      s.prior_mu, s.prior_sd, s.sigma_prior, s.abstol, s.reltol,
                      s.maxiters, s.box_sd)
end

"""
    mixed_parameter_names(s) -> Vector{Symbol}

The coordinates of x: `log_<target>` for each target, then `log_sigma_b`.
"""
mixed_parameter_names(s::MixedSetup) = _log_names(s.names)

# The chain coordinates for a list of targets, shared with `infer`'s exact route.
_log_names(targets) = vcat([Symbol(:log_, n) for n in targets], :log_sigma_b)

"""
    mixed_prior_median(s) -> Vector{Float64}

x at the prior median: each log prior's mean.
"""
mixed_prior_median(s::MixedSetup) = vcat(s.prior_mu, s.sigma_prior.μ)

"""
    mixed_box(s) -> (lower, upper)

The support box in x: each coordinate's prior mean ± `box_sd` prior SDs.
"""
function mixed_box(s::MixedSetup)
    sd = vcat(s.prior_sd, s.sigma_prior.σ)
    return (mixed_prior_median(s) .- s.box_sd .* sd, mixed_prior_median(s) .+ s.box_sd .* sd)
end

# Every parameter write for θ, each reverse constant derived through Haldane.
function _mixed_writes(s::MixedSetup, θ::AbstractVector)
    plain = [n => v for (n, v) in zip(s.names, θ) if !(n in s.forwards)]
    derived = derived_ode_values(s.models, [n => v for (n, v) in zip(s.names, θ) if n in s.forwards])
    return vcat(plain, derived)
end

# Route (ii) at θ, saved at t = 0, the switches and the saves; `nothing` when the
# route cannot run there (an overdraw, or a solve that does not succeed).
function _mixed_solve(s::MixedSetup, θ::AbstractVector)
    p = mean_field_parameters(s.route, _mixed_writes(s, θ))
    saveat = sort(unique(vcat(s.switches, s.saves)))
    sol = try
        solve_mean_field(s.route; p, saveat, abstol = s.abstol, reltol = s.reltol,
                         maxiters = s.maxiters)
    catch e
        e isa ErrorException && startswith(e.msg, "Route (ii): the pool") && return nothing
        rethrow()
    end
    successful_retcode(sol) || return nothing
    return (; sol, p)
end

_state_at(sol, t) = sol.u[findfirst(τ -> isapprox(τ, t; atol = 1e-9), sol.t)]

"""
    route_transcription_rates(s, sol, p) -> Matrix

The transcription constants route (ii) ran, genes × (windows + 1), in the order
[`transcript_loglik`](@ref) reads them. Column 1 is rebuilt from the pools at
t = 0, and column m + 1 from those at the switch `60m − 1`. `sol` must be saved
at both. Generic in the solution's number type.
"""
function route_transcription_rates(s::MixedSetup, sol, p)
    rt = s.route
    cols = map(s.switches) do t
        u = _state_at(sol, t)
        vals = Dict(i => _route_rebuilt_values(rt.rebuilds[i], p, u, t)
                    for i in unique(first.(s.tx_slots)))
        [vals[i][j] for (i, j) in s.tx_slots]
    end
    return reduce(hcat, cols)
end

"""
    mean_field_log_bulk(s, sol) -> Vector

`log max(x̄, floor)` on the panel at the saves, pools fastest, with x̄ in
particles at the route's live volume.
"""
function mean_field_log_bulk(s::MixedSetup, sol)
    rt = s.route
    out = map(s.saves) do t
        u = _state_at(sol, t)
        f = _route_geometry(rt, rt.n_idx == 0 ? zero(eltype(u)) : u[rt.n_idx]).factor
        [log(max(u[i] * f, s.floor)) for i in s.rows]
    end
    return reduce(vcat, out)
end

"""
    bulk_loglik(s, m, σ) -> Real

`log N(log z; m + δ, σ² I + S)`, one Cholesky per call.
"""
function bulk_loglik(s::MixedSetup, m::AbstractVector, σ::Real)
    r = s.logz .- m .- s.delta
    Σ = s.S + σ^2 * I
    C = cholesky(Symmetric(Σ))
    return -(dot(r, C \ r) + logdet(C) + length(r) * log(2π)) / 2
end

"""
    mixed_loglik(s, x) -> Real

The transcript and bulk factors at x = (log θ, log σ_b), without the prior. −Inf
outside the support box ([`mixed_box`](@ref)), where route (ii) cannot run or
reaches its step cap, or where a rate leaves the support the factors were sized
for (§12 2026-10-23, "B5 as planned" and "B5's support box").
"""
function mixed_loglik(s::MixedSetup, x::AbstractVector)
    n = length(s.names)
    lo, hi = mixed_box(s)
    all(i -> lo[i] <= _primal(x[i]) <= hi[i], eachindex(lo)) || return oftype(exp(x[n + 1]), -Inf)
    θ = exp.(x[1:n])
    σ = exp(x[n + 1])
    r = _mixed_solve(s, θ)
    r === nothing && return oftype(σ, -Inf)
    k = route_transcription_rates(s, r.sol, r.p)
    μ = θ[findfirst(==(:krnadeg), s.names)] ./ s.lengths
    f1 = first(s.factors)
    for g in eachindex(μ)
        _primal(μ[g]) >= f1.μmin[g] || return oftype(σ, -Inf)
        all(j -> _primal(k[g, j]) <= f1.kmax[g], axes(k, 2)) || return oftype(σ, -Inf)
    end
    return transcript_loglik(s.factors, k, μ) + bulk_loglik(s, mean_field_log_bulk(s, r.sol), σ)
end

"""
    mixed_log_prior(s, x) -> Real

The prior on x: a `Normal` on each log target, and σ_b's on log σ_b.
"""
mixed_log_prior(s::MixedSetup, x::AbstractVector) =
    sum(logpdf(Normal(μ, sd), xi) for (μ, sd, xi) in zip(s.prior_mu, s.prior_sd, x)) +
    logpdf(s.sigma_prior, x[end])

"""
    mixed_logdensity(s, x) -> Real

The approximate posterior's log density at x = (log θ, log σ_b), up to a constant.
"""
function mixed_logdensity(s::MixedSetup, x::AbstractVector)
    ll = mixed_loglik(s, x)
    isfinite(ll) || return ll
    return mixed_log_prior(s, x) + ll
end

@model function _mixed_model(s::MixedSetup)
    x ~ MvNormal(mixed_prior_median(s), Diagonal(vcat(s.prior_sd, s.sigma_prior.σ) .^ 2))
    Turing.@addlogprob! mixed_loglik(s, x)
end

"""
    mixed_model(s) -> DynamicPPL.Model

The approximate posterior as a Turing model over the vector `x`, ready for
`sample(mixed_model(s), NUTS(), n)`. The prior is [`mixed_log_prior`](@ref)'s, and
`x[i]` is [`mixed_parameter_names`](@ref)`(s)[i]`.
"""
mixed_model(s::MixedSetup) = _mixed_model(s)

"""
    reference_theta(s; x0 = mixed_prior_median(s), maxiters = 200) -> NamedTuple

Pass 1 of the reference θ (§12 2026-10-23, "B5 as planned"): the MAP of
`mixed_model(s)` from `x0`, bounded by the support box ([`mixed_box`](@ref)),
where `s` normally has δ = 0 and S = 0. If the optimiser fails, stops at
`maxiters` without converging, or ends at a non-finite density, the prior median
is used. Returns `(x, method, logdensity)`, with `method` `:map` or `:prior_median`.
"""
function reference_theta(s::MixedSetup; x0::AbstractVector = mixed_prior_median(s),
                         maxiters::Integer = 200)
    res = try
        lb, ub = mixed_box(s)
        maximum_a_posteriori(mixed_model(s); initial_params = collect(Float64, x0), lb, ub,
                             maxiters)
    catch e
        e isa InterruptException && rethrow()
        @warn "reference_theta: the MAP failed, so the prior median is used" exception = e
        nothing
    end
    if res !== nothing && successful_retcode(res.optim_result) && isfinite(res.lp)
        return (x = collect(Float64, res.values.array), method = :map, logdensity = res.lp)
    end
    x = mixed_prior_median(s)
    return (x = x, method = :prior_median, logdensity = mixed_logdensity(s, x))
end

"""
    mixed_theta_values(s, x) -> Vector{Pair{Symbol, Float64}}

Every parameter write for x's θ, reverse constants derived: what the reference
composition runs at.
"""
mixed_theta_values(s::MixedSetup, x::AbstractVector) =
    [n => Float64(v) for (n, v) in _mixed_writes(s, exp.(x[1:length(s.names)]))]

"""
    reference_ensemble(values, seeds; models = m0_models(), horizon = M0_HORIZON_S,
                       panel = M0_PANEL) -> Array{Float64, 3}

The reference composition at `values`, one cell per seed, as
[`generate_dataset`](@ref) runs a dataset's cells. Returns the panel in particles,
cells × pools × saves, without t = 0.
"""
function reference_ensemble(values, seeds::AbstractVector{<:Integer};
                            models = m0_models(), horizon::Real = M0_HORIZON_S,
                            panel::AbstractVector{Symbol} = M0_PANEL)
    ms = collect(AbstractSubModel, models)
    assert_haldane(ms, values)
    function build()
        d = build_problem(ms; tspan = (0.0, Float64(horizon)), complete = true)
        set_parameters!(d, ms, values)
        return ms, d
    end
    e = emit_ensemble(build, seeds)
    keep = findall(>(0), e.times)
    return e.counts[:, state_rows(collect(panel), e.species), keep]
end

const REFERENCE_SEED_OFFSET = 89_000_000_000

"""
    reference_seeds(ds, n) -> Vector{Int}

The `n` cell seeds of `ds`'s reference ensemble (§12 2026-10-23, "B5 as planned"):
drawn from a generator seeded at the dataset's replicate seed plus
`REFERENCE_SEED_OFFSET`, so each dataset has its own range. They are refused if
one repeats or meets one of the dataset's cell seeds.
"""
function reference_seeds(ds, n::Integer)
    seeds = Int.(rand(Xoshiro(REFERENCE_SEED_OFFSET + ds.replicate_seed), UInt32, n))
    allunique(seeds) && isdisjoint(seeds, ds.seeds) || error(
        "replicate $(ds.replicate_seed)'s reference seeds repeat or meet its cells'")
    return seeds
end

"""
    delta_and_S(cells, m_ref; C) -> (delta, S)

δ and S from a reference ensemble `cells` (cells × pools × saves, in particles)
and the mean field's `m_ref` = [`mean_field_log_bulk`](@ref) at the same θ, pools
fastest (§12 2026-10-23, "B5 as planned"):
- δ = `log max(x̄, 1) − m_ref`, with x̄ the mean over every cell;
- S = Cov(per-cell counts) / (x̄ x̄ᵀ C), the delta method's covariance of
  `log x̄_C` for a C-cell mean.
"""
function delta_and_S(cells::AbstractArray{<:Real, 3}, m_ref::AbstractVector; C::Integer)
    X = reshape(cells, size(cells, 1), :)
    size(X, 2) == length(m_ref) || throw(DimensionMismatch(
        "the ensemble has $(size(X, 2)) pool × save entries and m_ref $(length(m_ref))"))
    size(X, 1) > 1 || throw(ArgumentError("S needs at least two cells"))
    x̄ = vec(mean(X; dims = 1))
    δ = log.(max.(x̄, 1.0)) .- m_ref
    S = cov(X) ./ (x̄ * x̄') ./ C
    return (delta = δ, S = Matrix(Symmetric(S)))
end

export MixedSetup, mixed_setup, with_reference, mixed_parameter_names, mixed_prior_median,
       mixed_box,
       route_transcription_rates, mean_field_log_bulk, bulk_loglik, mixed_loglik,
       mixed_log_prior, mixed_logdensity, mixed_model, reference_theta, mixed_theta_values,
       reference_ensemble, reference_seeds, delta_and_S
