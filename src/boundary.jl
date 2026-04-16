import KernelDensity

"""
    KDEPrior <: ContinuousUnivariateDistribution

Univariate distribution fitted to posterior samples via kernel density estimation.
Satisfies the `rand(rng, d)` and `pdf(d, x)` interface required by `abc_smc()`.
"""
struct KDEPrior <: ContinuousUnivariateDistribution
    samples::Vector{Float64}
    kde_result::KernelDensity.UnivariateKDE
    bandwidth::Float64
end

function KDEPrior(samples::Vector{Float64})
    k = KernelDensity.kde(samples)
    # Silverman's rule of thumb bandwidth
    h = 1.06 * std(samples) * length(samples)^(-1/5)
    return KDEPrior(samples, k, h)
end

function Distributions.rand(rng::Random.AbstractRNG, d::KDEPrior)
    i = rand(rng, 1:length(d.samples))
    return d.samples[i] + d.bandwidth * randn(rng)
end

function Distributions.pdf(d::KDEPrior, x::Real)
    # Evaluate the KDE on a single point via interpolation from the grid
    grid = d.kde_result.x
    dens = d.kde_result.density
    # Clamp to grid range
    x < first(grid) && return 0.0
    x > last(grid) && return 0.0
    # Linear interpolation
    step = grid[2] - grid[1]
    idx_f = (x - first(grid)) / step + 1.0
    idx_lo = floor(Int, idx_f)
    idx_lo = clamp(idx_lo, 1, length(dens) - 1)
    frac = idx_f - idx_lo
    return dens[idx_lo] + frac * (dens[idx_lo + 1] - dens[idx_lo])
end

function Distributions.logpdf(d::KDEPrior, x::Real)
    p = pdf(d, x)
    return p > 0 ? log(p) : -Inf
end

"""
    boundary_condition(chain, ssa_models; method=:kde)

Build conditioned priors for the SSA block from an ODE posterior chain.

For each free parameter in `ssa_models`, if a matching column exists in `chain`
(by name), its prior is replaced with a distribution fitted to the posterior
samples. Parameters not found in the chain keep their original prior.

Returns a named tuple `(priors, param_names, conditioned)`.
"""
function boundary_condition(chain, ssa_models::Vector{<:AbstractSubModel};
                            method::Symbol=:kde)
    all_free = unique_params(model_free_params(reduce(vcat, parameters.(ssa_models))))
    chain_names = Set(String.(names(chain, :parameters)))

    priors = Distribution[]
    param_names = Symbol[]
    conditioned = Bool[]

    for p in all_free
        pname = string(p.name)
        if pname in chain_names
            samples = vec(chain[pname].data)
            if method == :kde
                push!(priors, KDEPrior(samples))
            elseif method == :normal
                push!(priors, Normal(mean(samples), std(samples)))
            else
                error("Unknown boundary method: $method. Use :kde or :normal.")
            end
            push!(conditioned, true)
        else
            push!(priors, p.prior)
            push!(conditioned, false)
        end
        push!(param_names, p.name)
    end

    return (priors=priors, param_names=param_names, conditioned=conditioned)
end

boundary_condition(chain, model::AbstractSubModel; kwargs...) =
    boundary_condition(chain, [model]; kwargs...)

"""
    sequential_infer(ode_models, ode_data, ssa_models, ssa_data; kwargs...)

Two-stage sequential inference (cut posterior):
1. Run NUTS on the ODE block
2. Fit boundary priors from ODE posterior for shared parameters
3. Run ABC-SMC on the SSA block with conditioned priors

Returns `(ode_chain, ssa_posterior, boundary)`.
"""
function sequential_infer(ode_models::Vector{<:AbstractSubModel}, ode_data::ObservedData,
                          ssa_models::Vector{<:AbstractSubModel}, ssa_data::ObservedData;
                          ode_kwargs=(;), ssa_kwargs=(;),
                          boundary_method::Symbol=:kde, verbose::Bool=false)
    # Stage 1: ODE inference via NUTS
    verbose && println("Stage 1: ODE inference (NUTS)...")
    ode_chain = infer(ode_models, ode_data; ode_kwargs...)

    # Stage 2: Boundary conditioning
    verbose && println("Stage 2: Boundary conditioning (method=$boundary_method)...")
    bc = boundary_condition(ode_chain, ssa_models; method=boundary_method)

    n_cond = sum(bc.conditioned)
    verbose && println("  Conditioned $n_cond / $(length(bc.param_names)) SSA parameters")

    # Stage 3: SSA inference via ABC-SMC with conditioned priors
    verbose && println("Stage 3: SSA inference (ABC-SMC with conditioned priors)...")
    merged_ssa = merge((verbose=verbose,), ssa_kwargs,
                       (priors=bc.priors, param_names=bc.param_names))
    ssa_posterior = _infer_abc(ssa_models, ssa_data; merged_ssa...)

    return (ode_chain=ode_chain, ssa_posterior=ssa_posterior, boundary=bc)
end

sequential_infer(ode_model::AbstractSubModel, ode_data::ObservedData,
                 ssa_model::AbstractSubModel, ssa_data::ObservedData; kwargs...) =
    sequential_infer([ode_model], ode_data, [ssa_model], ssa_data; kwargs...)
