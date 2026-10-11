"""
The phase 16 sampler composed: one Gibbs chain over the three blocks
(spec/phases/16-recovery.md §3, D16.1 to D16.3, task 16a.9a).

One sweep is a systematic scan:
1. **Block 3,** cells in turn. Each cell's path is updated by particle Gibbs
   ([`csmc_sweep`](@ref)) given θ, σ_b and the other cells' current latents,
   which enter through the bulk likelihood of the population mean. The others'
   sum is refreshed after each cell.
2. **Block 1,** the free forward ODE constants and σ_b given every path
   ([`block1_update!`](@ref)), with the reverse constants derived.
3. **Block 2,** each promoter and `krnadeg` given the paths and the pools, from
   block 1's accepted replay: firings and exposures per unit θ
   ([`rate_update`](@ref)). Given the path the pools do not depend on these, so
   the replay at the old value serves, and its latents are the cells' latents
   at the new θ, which the next sweep's block 3 reads through the bulk term.
4. Every cell's t = 0 driver is rebuilt at the new θ by writing all of θ onto a
   pristine snapshot, so transcription's starting constants are re-derived.

A chain starts from a path drawn by plain SMC with transcript weights only, at a
starting θ. Any path that meets every observed count has positive density, so
that is a valid start.
"""

"""
    GibbsSetup

What a chain holds fixed: the composition, each cell's pristine t = 0 driver,
the transcript map, the panel's ODE rows, each cell's transcripts (genes ×
windows), the bulk record (panel × windows), the particle count, and the free
parameters with their priors. Built by [`gibbs_setup`](@ref).
"""
struct GibbsSetup
    models::Vector{AbstractSubModel}
    pristine::Vector{HandshakeDriver}
    tmap::TranscriptMap
    panel::Vector{Symbol}
    rows::Vector{Int}
    transcripts::Vector{Matrix{Int}}
    bulk::Matrix{Float64}
    window::Int
    floor::Float64
    cap_floor::Int
    N::Int
    cme::Vector{Symbol}
    cme_priors::Vector{LogNormal{Float64}}
    forwards::Vector{Symbol}
    sigma_prior::LogNormal{Float64}
end

"""
    gibbs_setup(ds; models = m0_models(), build = build_m0, N = 20, horizon,
                panel = M0_PANEL, forwards = [:kcatF_R_ENO, :kcatF_R_FBA],
                window = 60, floor = 1.0, sigma_prior = SIGMA_MET_PRIOR) -> GibbsSetup

A chain's fixed ingredients for a dataset `ds` with a bulk record, as
[`m0_dataset`](@ref) returns. Each cell's driver is built as the dataset's was,
under its seed, at its nominal parameters.
"""
function gibbs_setup(ds; models = m0_models(), build = build_m0, N::Integer = 20,
                     horizon::Real = M0_HORIZON_S, panel::AbstractVector{Symbol} = M0_PANEL,
                     forwards::AbstractVector{Symbol} = [:kcatF_R_ENO, :kcatF_R_FBA],
                     window::Integer = 60, floor::Real = 1.0,
                     sigma_prior::LogNormal = SIGMA_MET_PRIOR)
    pristine = map(ds.seeds) do seed
        Random.seed!(seed)
        build(; tspan = (0.0, Float64(horizon)))
    end
    tmap = transcript_map(first(pristine), models)
    tx = only(m for m in models if m isa CoreATranscription)
    loci = [g.locus for g in tx.genes]
    gi = [findfirst(==(l), ds.genes) for l in loci]
    any(isnothing, gi) && error("the dataset does not record every transcribed gene")
    transcripts = [Matrix{Int}(ds.transcripts[c, gi, :]) for c in eachindex(ds.seeds)]
    T = size(first(transcripts), 2)
    size(ds.bulk) == (length(panel), T) || throw(DimensionMismatch(
        "the bulk record is $(size(ds.bulk)), not $(length(panel)) pools × $T windows"))
    rows = state_rows(collect(panel), _block_names(models, :ode))
    cme = vcat(promoter_param.(loci), :krnadeg)
    cme_priors = map(n -> _parameter_by_name(models, n).prior, cme)
    all(p -> p isa LogNormal, cme_priors) || throw(ArgumentError("block 2 needs LogNormal priors"))
    cap_floor = bridge_cap(maximum(maximum, transcripts))
    return GibbsSetup(collect(models), pristine, tmap, collect(panel), rows, transcripts,
                      Matrix{Float64}(ds.bulk), Int(window), Float64(floor), cap_floor, Int(N),
                      cme, [LogNormal(Float64.(params(p))...) for p in cme_priors],
                      collect(forwards), LogNormal(Float64.(params(sigma_prior))...))
end

"""
    GibbsState

A chain's state: the stochastic block's free values `cme` (in
[`GibbsSetup`](@ref)'s `cme` order), `u = ln θF`, σ_b, each cell's path by
window and its panel latent at each window's end, each cell's t = 0 driver at
the current θ, the sweep count, the chain's generator, and the trace of past
sweeps.
"""
mutable struct GibbsState
    cme::Vector{Float64}
    u::Vector{Float64}
    σ::Float64
    paths::Vector{Vector{WindowEvents}}
    latents::Vector{Matrix{Float64}}
    bases::Vector{HandshakeDriver}
    sweep::Int
    rng::Xoshiro
    trace::Vector{NamedTuple}
end

"""
    gibbs_writes(g, cme, u) -> Vector{Pair{Symbol, Float64}}

Every parameter write for a θ: the stochastic block's values, and the forward ODE
constants `exp.(u)` with their reverse constants derived.
"""
gibbs_writes(g::GibbsSetup, cme::AbstractVector, u::AbstractVector) =
    vcat([n => Float64(v) for (n, v) in zip(g.cme, cme)],
         derived_ode_values(g.models, [f => exp(x) for (f, x) in zip(g.forwards, u)]))

# Each cell's t = 0 driver at θ: a pristine snapshot with all of θ written.
function _gibbs_bases(g::GibbsSetup, cme, u; threads::Bool)
    writes = gibbs_writes(g, cme, u)
    bases = Vector{HandshakeDriver}(undef, length(g.pristine))
    _each(length(bases), threads) do c
        d = restore(g.pristine[c])
        set_parameters!(d, g.models, writes)
        bases[c] = d
    end
    return bases
end

_gibbs_spec(g::GibbsSetup, σ) =
    CSMCSpec(g.tmap, g.rows; sigma = σ, floor = g.floor, window = g.window,
             cap_floor = g.cap_floor)

"""
    gibbs_init(rng, g; cme, u, σ, threads = Threads.nthreads() > 1) -> GibbsState

A chain at `(cme, u, σ)`, each cell's path drawn by plain SMC at that θ with
transcript weights only ([`csmc_sweep`](@ref) with `ref = nothing`). The chain's
generator is seeded from `rng`.
"""
function gibbs_init(rng::AbstractRNG, g::GibbsSetup; cme::AbstractVector, u::AbstractVector,
                    σ::Real, threads::Bool = Threads.nthreads() > 1)
    length(cme) == length(g.cme) && length(u) == length(g.forwards) ||
        throw(DimensionMismatch("θ does not match the setup's free parameters"))
    chain = Xoshiro(rand(rng, UInt64))
    bases = _gibbs_bases(g, cme, u; threads)
    spec = _gibbs_spec(g, σ)
    C = length(bases)
    paths = Vector{Vector{WindowEvents}}(undef, C)
    latents = Vector{Matrix{Float64}}(undef, C)
    for c in 1:C
        data = CellData(g.transcripts[c], g.bulk)
        paths[c], _, latents[c] = csmc_sweep(chain, bases[c], spec, data, nothing; N = g.N,
                                            threads, metabolites = false)
    end
    return GibbsState(collect(Float64, cme), collect(Float64, u), Float64(σ), paths, latents,
                      bases, 0, chain, NamedTuple[])
end

"""
    gibbs_sweep!(s, g; threads = Threads.nthreads() > 1) -> s

One sweep of the chain: block 3 cell by cell, block 1, block 2, then the cells'
drivers rebuilt at the new θ. Appends to `s.trace` the new values, which windows
of which cells changed, and each block's seconds.
"""
function gibbs_sweep!(s::GibbsState, g::GibbsSetup; threads::Bool = Threads.nthreads() > 1)
    rng = s.rng
    C = length(s.paths)
    T = size(g.bulk, 2)
    changed = falses(T, C)
    t3 = @elapsed begin
        spec = _gibbs_spec(g, s.σ)
        total = sum(s.latents)
        for c in 1:C
            others = total .- s.latents[c]
            bulk = bulk_window_observations(g.bulk, others, C)
            data = CellData(g.transcripts[c], g.bulk)
            s.paths[c], ch, lat = csmc_sweep(rng, s.bases[c], spec, data, s.paths[c];
                                             N = g.N, bulk, threads)
            total .+= lat .- s.latents[c]
            s.latents[c] = lat
            changed[:, c] = ch
        end
    end
    labels = first(s.bases).events.labels
    st = Block1State(s.u, s.σ, 0.0, 0)
    t1 = @elapsed begin
        cells = [Block1Cell(s.bases[c], join_path(s.paths[c], labels), zeros(0, 0)) for c in 1:C]
        b = Block1(g.models, cells; forwards = g.forwards, panel = g.panel,
                   save_every = g.window, floor = g.floor, sigma_prior = g.sigma_prior,
                   bulk = g.bulk)
        block1_update!(rng, b, st; threads)
    end
    t2 = @elapsed begin
        # Block 2 reads block 1's accepted replay, which ran at the current θ_CME.
        cme = copy(s.cme)
        bound = TURNOVER_CEILING / RNAPOL_KCAT
        for (k, j) in enumerate(g.tmap.births)
            cme[k] = rate_update(rng, cme[k], st.counts[j], st.exposure[j], g.cme_priors[k];
                                 bound)
        end
        jd = g.tmap.deaths
        cme[end] = rate_update(rng, cme[end], sum(st.counts[jd]), sum(st.exposure[jd]),
                               g.cme_priors[end])
    end
    tb = @elapsed begin
        s.cme = cme
        s.u = copy(st.u)
        s.σ = st.σ
        # The latents block 3 left were at the old θF. Block 1's accepted replay
        # has them at the new one, and block 2's values do not move the pools
        # given the path, so these are the latents at the chain's new θ.
        s.latents = st.latents
        s.bases = _gibbs_bases(g, s.cme, s.u; threads)
    end
    s.sweep += 1
    push!(s.trace, (sweep = s.sweep, cme = copy(s.cme), u = copy(s.u), σ = s.σ,
                    changed = changed, seconds = (b3 = t3, b1 = t1, b2 = t2, rebuild = tb)))
    return s
end

"""
    gibbs_checkpoint(s) -> NamedTuple

The chain's state as plain data, everything but the drivers, which
[`gibbs_resume`](@ref) rebuilds from θ. A script serializes it.
"""
gibbs_checkpoint(s::GibbsState) =
    (; cme = copy(s.cme), u = copy(s.u), s.σ, paths = deepcopy(s.paths),
     latents = deepcopy(s.latents), s.sweep, rng = copy(s.rng), trace = copy(s.trace))

"""
    gibbs_resume(saved, g; threads = Threads.nthreads() > 1) -> GibbsState

A chain from a [`gibbs_checkpoint`](@ref), its drivers rebuilt at the saved θ.
Continuing it gives what the uninterrupted chain would have.
"""
function gibbs_resume(r::NamedTuple, g::GibbsSetup; threads::Bool = Threads.nthreads() > 1)
    return GibbsState(copy(r.cme), copy(r.u), r.σ, deepcopy(r.paths), deepcopy(r.latents),
                      _gibbs_bases(g, r.cme, r.u; threads), r.sweep, copy(r.rng), copy(r.trace))
end

export GibbsSetup, gibbs_setup, GibbsState, gibbs_writes, gibbs_init, gibbs_sweep!,
       gibbs_checkpoint, gibbs_resume
