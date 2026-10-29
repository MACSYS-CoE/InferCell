# Posterior draws for `predict` (KIM-28): one type, whatever produced the posterior.
#
# Four sources write a posterior over named model parameters, each in its own
# layout:
# - the exact route's V7 replicate files, `(; state, meta)`, whose `state.trace`
#   holds one Gibbs sweep per entry (`dev/scripts/sbc_v7_16a9b.jl`);
# - the mixed route's B7 chain files, a `Chains` over `x[1]`…`x[6]` on the log
#   scale beside the run's record (`dev/scripts/av6_89_b7.jl`);
# - `infer(models, ds; route)`, a `Chains` over `log_<name>` columns;
# - the toy models' `infer(models, ::ObservedData)`, a `Chains` over
#   natural-scale names, or an `ABCPosterior`.
#
# `PosteriorDraws` holds all of them the same way: natural-scale values by
# parameter name, ready to write into a model, with the bulk noise σ_b kept
# apart. σ_b is observation noise, not a cell parameter, and `predict` never
# writes it to a model (dev/notes/2026-10-29-cut-predictive-distribution.md §1).

"""
    PosteriorDraws

A posterior's draws over named model parameters, as [`predict`](@ref) reads
them, whatever produced them. Build one with [`posterior_draws`](@ref).

Fields:
- `names`: the model parameters, as `set_parameters!` names them;
- `values`: draws × `names`, on the natural scale;
- `sigma_b`: the bulk observation noise of each draw, on the natural scale, or
  `nothing` when the posterior has none. It is not a cell parameter, and is never
  written to a model;
- `chain` and `iteration`: where each draw sits in the run that produced it;
- `weights`: normalised importance weights (an ABC posterior), or `nothing` for
  equally weighted draws;
- `truth`: the dataset's true values by name on the natural scale, σ_b under
  `:sigma_b`, when the source records them; otherwise `nothing`;
- `source`: what the draws were read from. `kind` is one of `:v7_replicate`,
  `:b7_chain`, `:chains`, `:abc` or `:merged`, beside what the source records.
"""
struct PosteriorDraws
    names::Vector{Symbol}
    values::Matrix{Float64}
    sigma_b::Union{Nothing, Vector{Float64}}
    chain::Vector{Int}
    iteration::Vector{Int}
    weights::Union{Nothing, Vector{Float64}}
    truth::Union{Nothing, Dict{Symbol, Float64}}
    source::NamedTuple
    function PosteriorDraws(names, values, sigma_b, chain, iteration, weights, truth, source)
        n = size(values, 1)
        size(values, 2) == length(names) || throw(ArgumentError(
            "$(length(names)) names for $(size(values, 2)) value columns"))
        allunique(names) || throw(ArgumentError("a parameter name repeats: $names"))
        :sigma_b in names && throw(ArgumentError(
            "σ_b is observation noise and goes in `sigma_b`, not among the parameters"))
        n > 0 || throw(ArgumentError("a posterior with no draws"))
        all(isfinite, values) || throw(ArgumentError("a draw is not finite"))
        length(chain) == n && length(iteration) == n || throw(ArgumentError(
            "`chain` and `iteration` must give one entry per draw"))
        sigma_b === nothing || length(sigma_b) == n || throw(ArgumentError(
            "`sigma_b` must give one entry per draw"))
        if weights !== nothing
            length(weights) == n && all(>=(0), weights) && sum(weights) > 0 ||
                throw(ArgumentError("`weights` must be one non-negative entry per draw"))
            weights = weights ./ sum(weights)
        end
        return new(collect(Symbol, names), Matrix{Float64}(values), sigma_b,
                   collect(Int, chain), collect(Int, iteration), weights, truth, source)
    end
end

Base.length(pd::PosteriorDraws) = size(pd.values, 1)

function Base.show(io::IO, pd::PosteriorDraws)
    print(io, "PosteriorDraws($(length(pd)) draws over $(pd.names)",
          pd.sigma_b === nothing ? "" : ", σ_b carried",
          pd.weights === nothing ? "" : ", weighted", "; ", pd.source.kind, ")")
end

"""
    draw_values(pd, i) -> Vector{Pair{Symbol, Float64}}

Draw `i` as `name => value` pairs on the natural scale. σ_b is not among them.
"""
draw_values(pd::PosteriorDraws, i::Integer) = [n => pd.values[i, j] for (j, n) in enumerate(pd.names)]

# A column name's parameter and scale: `log_<name>` is on the log scale.
function _decode_column(c::Symbol)
    s = String(c)
    return startswith(s, "log_") ? (Symbol(s[5:end]), :log) : (c, :natural)
end

"""
    posterior_draws(chain::Chains; names = nothing, truth = nothing,
                    source = (kind = :chains,)) -> PosteriorDraws
    posterior_draws(result::ABCPosterior) -> PosteriorDraws
    posterior_draws(path::AbstractString; burn = 100, thin = 1) -> PosteriorDraws
    posterior_draws(paths::AbstractVector{<:AbstractString}; kwargs...) -> PosteriorDraws

Read a posterior into a [`PosteriorDraws`](@ref), mapping every column to its
parameter by name.

**A `Chains`** reads its `:parameters` section and leaves the sampler's internals.
A column `log_<name>` is `<name>` on the log scale, as `infer(models, ds; route)`
writes it, and any other column is on the natural scale, as the toy models' NUTS
path writes it. `sigma_b`, on either scale, becomes the `sigma_b` field. A chain
whose columns are `x[1]`, `x[2]`, … (Turing's names for a vector parameter) needs
`names`, the column names in index order. Every chain of the `Chains` is read,
not only the first.

**A file** is deserialized and read by what it holds:
- a V7 replicate file, `(; state, meta)`, from `dev/scripts/sbc_v7_16a9b.jl`. Its
  draws are the sweeps after the first `burn`, every `thin`-th. The default
  `burn = 100` is V7's. The sweeps come from one correlated chain (23 of the 50
  replicates have a smallest bulk ESS under 100), so `thin` is there to thin it.
  The names are `meta.names`, and the truth is `meta.truth`;
- a B7 chain file, from `dev/scripts/av6_89_b7.jl`. Its `x[i]` columns are the
  mixed route's names in order, `log_` and then `M0_TARGETS` plus σ_b. The truth
  comes from the fit's `prep.jls` beside it, when that file holds the dataset;
- a `Chains` saved on its own.

A `PosteriorDraws` is returned as it is.

**Several files** are read in turn and merged, as a fit's two B7 chain files are.
Their parameters must match, and so must their truths when they record one: a
posterior merged over datasets is refused.
"""
function posterior_draws(ch::Chains; names::Union{Nothing, AbstractVector{Symbol}} = nothing,
                         truth::Union{Nothing, AbstractDict} = nothing,
                         source::NamedTuple = (kind = :chains,))
    cols = if names === nothing
        Base.names(ch, :parameters)
    else
        xs = [Symbol("x[$i]") for i in eachindex(names)]
        have = Set(Base.names(ch, :parameters))
        all(in(have), xs) && length(have) == length(xs) || throw(ArgumentError(
            "`names` gives $(length(names)) columns, but the chain's parameters are " *
            "$(sort(collect(have)))"))
        xs
    end
    decoded = names === nothing ? _decode_column.(cols) : _decode_column.(collect(Symbol, names))
    keep = findall(d -> first(d) !== :sigma_b, decoded)
    sig = findfirst(d -> first(d) === :sigma_b, decoded)
    nit, _, nch = size(ch)
    its = collect(Int, range(ch))
    # The raw iterations × variables × chains array, chains stacked draw-wise.
    idx = [findfirst(==(c), Base.names(ch)) for c in cols]
    X = reduce(vcat, (Float64.(ch.value.data[:, idx, c]) for c in 1:nch))
    natural(j) = last(decoded[j]) === :log ? exp.(X[:, j]) : X[:, j]
    values = reduce(hcat, (natural(j) for j in keep); init = zeros(nit * nch, 0))
    return PosteriorDraws(first.(decoded[keep]), values,
                          sig === nothing ? nothing : natural(sig),
                          repeat(1:nch; inner = nit), repeat(its, nch), nothing,
                          truth === nothing ? nothing : Dict{Symbol, Float64}(truth),
                          merge(source, (route = get(ch.info, :route, nothing),)))
end

posterior_draws(pd::PosteriorDraws) = pd

function posterior_draws(r::ABCPosterior)
    n = size(r.particles, 2)
    return PosteriorDraws(r.param_names, permutedims(r.particles), nothing, fill(1, n),
                          collect(1:n), r.weights, nothing,
                          (kind = :abc, tolerance = r.tolerance,
                           n_populations = r.n_populations))
end

const V7_BURN = 100

function posterior_draws(path::AbstractString; burn::Integer = V7_BURN, thin::Integer = 1)
    x = deserialize(path)
    x isa Chains && return posterior_draws(x; source = (kind = :chains, path = String(path)))
    x isa NamedTuple && haskey(x, :state) && haskey(x, :meta) &&
        return _v7_draws(x, String(path); burn, thin)
    x isa NamedTuple && haskey(x, :chain) && x.chain isa Chains &&
        return _b7_draws(x, String(path))
    throw(ArgumentError(
        "$path holds a $(typeof(x)), not a V7 replicate file, a B7 chain file or a Chains"))
end

function posterior_draws(paths::AbstractVector{<:AbstractString}; kwargs...)
    isempty(paths) && throw(ArgumentError("no files to read"))
    return reduce(vcat, [posterior_draws(p; kwargs...) for p in paths])
end

# A V7 replicate: each sweep's draw is `vcat(cme, exp.(u), σ)` in `meta.names`
# order, the CME parameters on the natural scale and the forward constants as
# u = ln θ (`_infer_exact`, gibbs.jl).
function _v7_draws(x, path; burn, thin)
    names = collect(Symbol, x.meta.names)
    last(names) === :sigma_b || throw(ArgumentError(
        "$path: meta.names should end with :sigma_b, got $names"))
    trace = x.state.trace
    0 <= burn < length(trace) || throw(ArgumentError(
        "$path: burn = $burn leaves no draws of $(length(trace)) sweeps"))
    thin >= 1 || throw(ArgumentError("thin must be at least 1, got $thin"))
    kept = trace[(burn + 1):thin:end]
    rows = [vcat(t.cme, exp.(t.u)) for t in kept]
    all(r -> length(r) == length(names) - 1, rows) || throw(ArgumentError(
        "$path: a sweep's cme and u do not cover meta.names"))
    truth = Dict(n => exp(v) for (n, v) in zip(names, x.meta.truth))
    return PosteriorDraws(names[1:(end - 1)], permutedims(reduce(hcat, rows)),
                          [Float64(t.σ) for t in kept], fill(1, length(kept)),
                          [t.sweep for t in kept], nothing, truth,
                          (kind = :v7_replicate, path, replicate = x.meta.replicate,
                           seed = x.meta.seed, burn, thin))
end

# A B7 chain file: x[i] is `mixed_parameter_names(mixed_setup(ds))[i]`, which is
# `_log_names(M0_TARGETS)`. The truth is the fit's `prep.jls`, when it holds `ds`.
function _b7_draws(x, path)
    names = _log_names(M0_TARGETS)
    truth = nothing
    prep = joinpath(dirname(path), "prep.jls")
    if isfile(prep)
        p = deserialize(prep)
        if haskey(p, :ds)
            tv = Dict(p.ds.truth.values)
            truth = Dict{Symbol, Float64}(n => tv[n] for n in M0_TARGETS)
            truth[:sigma_b] = p.ds.sigma
        end
    end
    pd = posterior_draws(x.chain; names, truth,
                         source = (kind = :b7_chain, path, fit = x.fit, seed = x.seed,
                                   commit = x.commit))
    return PosteriorDraws(pd.names, pd.values, pd.sigma_b, fill(x.c, length(pd)),
                          pd.iteration, pd.weights, pd.truth, pd.source)
end

function Base.vcat(a::PosteriorDraws, bs::PosteriorDraws...)
    parts = (a, bs...)
    for b in bs
        b.names == a.names || throw(ArgumentError(
            "the posteriors name different parameters: $(a.names) and $(b.names)"))
        (b.sigma_b === nothing) == (a.sigma_b === nothing) || throw(ArgumentError(
            "one posterior carries σ_b and another does not"))
        (b.weights === nothing && a.weights === nothing) || throw(ArgumentError(
            "weighted posteriors cannot be merged by concatenation"))
        a.truth === nothing || b.truth === nothing || a.truth == b.truth || throw(ArgumentError(
            "the posteriors record different truths, so they are not one dataset's"))
    end
    k = findfirst(p -> p.truth !== nothing, parts)
    truth = k === nothing ? nothing : parts[k].truth
    return PosteriorDraws(a.names, reduce(vcat, (p.values for p in parts)),
                          a.sigma_b === nothing ? nothing : reduce(vcat, (p.sigma_b for p in parts)),
                          reduce(vcat, (p.chain for p in parts)),
                          reduce(vcat, (p.iteration for p in parts)), nothing, truth,
                          (kind = :merged, parts = [p.source for p in parts]))
end

export PosteriorDraws, posterior_draws, draw_values
