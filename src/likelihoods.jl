"""
    ObservedData(times, observations, species[, spread])

Time-series observations for inference: `times` is a length-T vector of sampling
times, `observations` is a `species × T` matrix of measured values, and
`species` names the rows. Returned by [`observe`](@ref) and consumed by
[`infer`](@ref).

Rows are matched to model states **by name**, not by position (spec §11 task
15.9): any subset of the states, in any order, may be observed, and a name that
is not a state throws.

`spread`, when not `nothing`, is a `species × T` matrix of the ensemble's
standard deviation at each time. [`observe`](@ref) fills it from replicate
trajectories, so the ABC summary can keep the spread rather than only the
per-time mean (§12, 2026-09-23 G).
"""
struct ObservedData
    times::Vector{Float64}
    observations::Matrix{Float64}  # species x timepoints
    species::Vector{Symbol}
    spread::Union{Nothing, Matrix{Float64}}  # species x timepoints, or nothing

    function ObservedData(times, observations, species, spread = nothing)
        size(observations, 1) == length(species) || throw(DimensionMismatch(
            "$(size(observations, 1)) observation rows for $(length(species)) species"))
        size(observations, 2) == length(times) || throw(DimensionMismatch(
            "$(size(observations, 2)) observation columns for $(length(times)) times"))
        spread === nothing || size(spread) == size(observations) || throw(
            DimensionMismatch("spread is $(size(spread)), observations $(size(observations))"))
        new(collect(Float64, times), Matrix{Float64}(observations), collect(Symbol, species),
            spread === nothing ? nothing : Matrix{Float64}(spread))
    end
end

"""
    state_rows(observed, names) -> Vector{Int}

The index in `names`, the model's state order, of each observed species. Throws
naming any species that is not a state.
"""
function state_rows(observed::AbstractVector{Symbol}, names::AbstractVector{Symbol})
    allunique(observed) || throw(ArgumentError("observed species repeat: $observed"))
    rows = indexin(observed, names)
    unknown = [o for (o, r) in zip(observed, rows) if r === nothing]
    isempty(unknown) || throw(ArgumentError(
        "observed species $(unknown) are not states of the model, whose states are $(names)"))
    return Int[r for r in rows]
end

"""
    PosteriorPredictive(solutions, times, species)

Container for forward simulations drawn from a posterior. `solutions` is a
vector of `DiffEq` / `JumpProcesses` solution objects evaluated at `times` for
the named `species`. Built by [`posterior_predictive`](@ref).
"""
struct PosteriorPredictive
    solutions::Vector
    times::Vector{Float64}
    species::Vector{Symbol}
end

# ---------------------------------------------------------------------------
# Per-modality observation noise (spec §11 task 15.4, §4 D11)
# ---------------------------------------------------------------------------

"The observation kinds a [`Modality`](@ref) may take."
const MODALITY_KINDS = (:lognormal, :poisson)

"""
    Modality(name, kind, species; prior = LogNormal(log(0.2), 1.0), floor = 1e-12)

One observation stream with its own noise model (spec §4 D11: per-modality noise
on a log scale, replacing the single additive term).

- `:lognormal`: `log y ~ Normal(log max(x, floor), σ)`. The scale `σ` is a free
  parameter named `sigma_<name>`, with `prior`. Because the noise is relative, a
  0.01 mM pool and a 17.8 mM pool are fitted to the same precision. One additive
  σ lets the small pools contribute nothing. The density is on `log y`, so it
  omits the `−log y` Jacobian: that term is constant in the parameters and σ,
  but it matters when comparing likelihoods across noise models.
- `:poisson`: `y ~ Poisson(max(x, floor))`, a count model with no scale, for
  transcript counts.
"""
struct Modality
    name::Symbol
    kind::Symbol
    species::Vector{Symbol}
    prior::Distribution
    floor::Float64

    function Modality(name::Symbol, kind::Symbol, species::AbstractVector{Symbol};
                      prior::Distribution = LogNormal(log(0.2), 1.0),
                      floor::Real = 1e-12)
        kind in MODALITY_KINDS || throw(ArgumentError(
            "modality kind must be one of $MODALITY_KINDS, not :$kind"))
        isempty(species) && throw(ArgumentError("modality :$name observes no species"))
        floor > 0 || throw(ArgumentError("the floor must be positive, got $floor"))
        new(name, kind, collect(Symbol, species), prior, Float64(floor))
    end
end

"The free scale parameter's name for a `:lognormal` modality."
scale_name(m::Modality) = Symbol(:sigma_, m.name)

"""
    NoiseModel(modalities...)

The observation model: a partition of the observed species into
[`Modality`](@ref) streams. Pass it to [`build_turing_model`](@ref) or
[`observe`](@ref) as `noise`.
"""
struct NoiseModel
    modalities::Vector{Modality}

    function NoiseModel(ms::AbstractVector{Modality})
        sp = reduce(vcat, [m.species for m in ms])
        allunique(sp) || throw(ArgumentError(
            "a species is in more than one modality: $(unique(filter(s -> count(==(s), sp) > 1, sp)))"))
        allunique([m.name for m in ms]) || throw(ArgumentError("modality names repeat"))
        new(collect(ms))
    end
end
NoiseModel(ms::Modality...) = NoiseModel(collect(ms))

"The free scale parameters, as `(name, prior)` pairs, in modality order."
noise_scales(nm::NoiseModel) =
    [(scale_name(m), m.prior) for m in nm.modalities if m.kind === :lognormal]

"""
    observation_loglik(nm, data, pred_at, scales) -> Real

The log-likelihood of `data` under `nm`. `pred_at(i)` returns the model's
predicted value of every observed species at `data.times[i]`, in `data.species`
order, and `scales` holds the `:lognormal` modalities' σ in
[`noise_scales`](@ref) order. Every species in `data` must be in exactly one
modality. Called by the Turing model, and directly by the tests that check a
known scale is recovered (spec §3 check 9).
"""
observation_loglik(nm::NoiseModel, data::ObservedData, pred_at, scales) =
    observation_loglik(nm, data, pred_at, scales, _modality_rows(nm, data))

# With the rows checked once, as `build_turing_model` does, so a sampler does not
# re-validate the data at every log-density evaluation.
function observation_loglik(nm::NoiseModel, data::ObservedData, pred_at, scales,
                            idx::AbstractVector{<:AbstractVector{Int}})
    ll = zero(eltype(scales))
    for i in eachindex(data.times)
        pred = pred_at(i)
        k = 0
        for (m, rows) in zip(nm.modalities, idx)
            if m.kind === :lognormal
                k += 1
                σ = scales[k]
                for r in rows
                    y = data.observations[r, i]
                    ll += logpdf(Normal(log(max(pred[r], m.floor)), σ), log(y))
                end
            else
                for r in rows
                    ll += logpdf(Poisson(max(pred[r], m.floor)),
                                 round(Int, data.observations[r, i]))
                end
            end
        end
    end
    return ll
end

# For each modality, the rows of `data` it observes. Checks the partition and
# that each row's values suit its kind.
function _modality_rows(nm::NoiseModel, data::ObservedData)
    covered = reduce(vcat, [m.species for m in nm.modalities])
    extra = setdiff(data.species, covered)
    isempty(extra) || throw(ArgumentError(
        "observed species $(extra) are in no modality of the noise model"))
    idx = Vector{Vector{Int}}()
    for m in nm.modalities
        rows = Int[r for r in indexin(m.species, data.species) if r !== nothing]
        for r in rows
            y = view(data.observations, r, :)
            m.kind === :lognormal && any(<=(0), y) && throw(ArgumentError(
                "lognormal modality :$(m.name) needs positive observations; " *
                "$(data.species[r]) has $(minimum(y))"))
            m.kind === :poisson && !all(v -> v >= 0 && isinteger(v), y) && throw(
                ArgumentError("poisson modality :$(m.name) needs counts; " *
                              "$(data.species[r]) is not integer-valued"))
        end
        push!(idx, rows)
    end
    return idx
end
