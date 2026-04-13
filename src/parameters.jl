struct InferParameter
    value::Float64
    prior::Distribution
    fixed::Bool
    name::Symbol
    module_id::Symbol
    role::Symbol  # :rate, :initial_condition, :observation
end

free_params(params::Vector{InferParameter}) = filter(p -> !p.fixed, params)
rate_params(params::Vector{InferParameter}) = filter(p -> p.role == :rate, params)
ic_params(params::Vector{InferParameter}) = filter(p -> p.role == :initial_condition, params)
obs_params(params::Vector{InferParameter}) = filter(p -> p.role == :observation, params)

model_free_params(params::Vector{InferParameter}) =
    filter(p -> p.role != :observation, free_params(params))

obs_free_params(params::Vector{InferParameter}) =
    filter(p -> p.role == :observation, free_params(params))

# Backward compatibility alias
const ode_free_params = model_free_params
