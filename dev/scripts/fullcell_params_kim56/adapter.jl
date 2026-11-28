# KIM-56: a parameter vector per run, for syn3a-ensemble's SpeedyCell (`Params`, lane 10a).
#
# `params_from_vector(model, v)` takes one multiplier per named kinetic constant, in
# `param_names()` order (805 ODE constants, then 33 CME constants; see param_names.csv),
# and returns the `Params` that `run_cell(...; params=P)` takes. v = ones(838) is the
# published model on the parameterised path. Validation (finite, ≥ 0) is Params's.

using SpeedyCell.Cell

function params_from_vector(m::CellModel, v::AbstractVector{<:Real}; names::AbstractVector{String}=param_names())
    length(v) == length(names) ||
        throw(DimensionMismatch("$(length(v)) multipliers for $(length(names)) named constants"))
    return Params(m; multipliers=Dict{String,Float64}(n => Float64(x) for (n, x) in zip(names, v)))
end
