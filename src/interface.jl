abstract type AbstractSubModel end

function states end
function parameters end
function dynamics end

inputs(::AbstractSubModel) = Symbol[]
formalism(::AbstractSubModel) = :ode
inference_mode(::AbstractSubModel) = :differentiable

reactions(m::AbstractSubModel) = error("reactions() not implemented for $(typeof(m))")

# Fallback: models with no inputs ignore u_inputs
dynamics(u, p, t, m::AbstractSubModel, u_inputs) = dynamics(u, p, t, m)

struct SubModelContext
    state_idxs::UnitRange{Int}
    param_idxs::Vector{Int}
    input_map::Dict{Symbol, Int}
end

const VALID_FORMALISMS = (:ode, :sde, :jump)

function validate_formalism(f::Symbol)
    f in VALID_FORMALISMS || throw(ArgumentError(
        "Invalid formalism :$f. Must be one of $VALID_FORMALISMS"))
end
