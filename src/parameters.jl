"""
    InferParameter(value, prior, fixed, name, module_id, role)

A single parameter declared by a sub-model: its nominal value, prior, identity,
and role. `role` is one of `:rate`, `:initial_condition`, `:observation`; the
inference layer uses `role` and `fixed` to decide which parameters are sampled.
"""
struct InferParameter
    value::Float64
    prior::Distribution
    fixed::Bool
    name::Symbol
    module_id::Symbol
    role::Symbol  # :rate, :initial_condition, :observation
end

"""
    free_params(params)

Return the subset of `params` that are not fixed (i.e. inferred).
"""
free_params(params::Vector{InferParameter}) = filter(p -> !p.fixed, params)

"""
    rate_params(params)

Return parameters whose `role` is `:rate` (kinetic rate constants).
"""
rate_params(params::Vector{InferParameter}) = filter(p -> p.role == :rate, params)

"""
    ic_params(params)

Return parameters whose `role` is `:initial_condition` (initial state values).
"""
ic_params(params::Vector{InferParameter}) = filter(p -> p.role == :initial_condition, params)

"""
    obs_params(params)

Return parameters whose `role` is `:observation` (e.g. measurement noise σ).
"""
obs_params(params::Vector{InferParameter}) = filter(p -> p.role == :observation, params)

"""
    model_free_params(params)

Free parameters of the dynamical model — everything except observation
parameters. These are the entries that appear in the ODE/SSA parameter vector.
"""
model_free_params(params::Vector{InferParameter}) =
    filter(p -> p.role != :observation, free_params(params))

"""
    obs_free_params(params)

Free observation parameters (e.g. measurement noise σ). These appear in the
likelihood but not in the dynamics RHS.
"""
obs_free_params(params::Vector{InferParameter}) =
    filter(p -> p.role == :observation, free_params(params))

"""
    ode_free_params(params)

Backward-compatibility alias for [`model_free_params`](@ref).
"""
const ode_free_params = model_free_params

"""
    unique_params(params)

Deduplicate `params` by `name`, preserving first-occurrence order. Used when
merging parameter lists from several sub-models that share named parameters.
"""
function unique_params(params::Vector{InferParameter})
    seen = Set{Symbol}()
    result = InferParameter[]
    for p in params
        if !(p.name in seen)
            push!(seen, p.name)
            push!(result, p)
        end
    end
    return result
end
