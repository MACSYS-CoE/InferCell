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
