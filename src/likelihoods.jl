"""
    ObservedData(times, observations, species)

Time-series observations for inference: `times` is a length-T vector of sampling
times, `observations` is a `species × T` matrix of measured values, and
`species` names the rows. Returned by [`observe`](@ref) and consumed by
[`infer`](@ref).
"""
struct ObservedData
    times::Vector{Float64}
    observations::Matrix{Float64}  # species x timepoints
    species::Vector{Symbol}
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
