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

# Bijectors / Turing need a declared support to derive a transformation.
# The KDE grid is finite; bound the support to the grid range.
Distributions.minimum(d::KDEPrior) = first(d.kde_result.x)
Distributions.maximum(d::KDEPrior) = last(d.kde_result.x)
Distributions.insupport(d::KDEPrior, x::Real) =
    first(d.kde_result.x) <= x <= last(d.kde_result.x)

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
    boundary_condition(posterior::ABCPosterior, ode_models; method=:kde, n_samples=1000, rng=...)

Build conditioned priors for the ODE block from an SSA `ABCPosterior` — the
return direction of the boundary protocol. Weighted particles are resampled to
unweighted samples; the rest mirrors the chain-based overload.
"""
function boundary_condition(posterior::ABCPosterior,
                            ode_models::Vector{<:AbstractSubModel};
                            method::Symbol=:kde, n_samples::Int=1000,
                            rng::Random.AbstractRNG=Random.default_rng())
    all_free = unique_params(model_free_params(reduce(vcat, parameters.(ode_models))))
    abc_names = posterior.param_names

    priors = Distribution[]
    param_names = Symbol[]
    conditioned = Bool[]

    for p in all_free
        j = findfirst(==(p.name), abc_names)
        if j !== nothing
            idxs = [_weighted_sample(posterior.weights, rng) for _ in 1:n_samples]
            samples = Float64[posterior.particles[j, k] for k in idxs]
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

boundary_condition(posterior::ABCPosterior, model::AbstractSubModel; kwargs...) =
    boundary_condition(posterior, [model]; kwargs...)

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


# --- Level-2 boundary protocol: iterative message-passing ---

"""
    kl_divergence(samples_p, samples_q; n_grid=200)

Estimate KL(p || q) where p, q are univariate distributions represented by
samples. Both are KDE-smoothed on a shared grid covering the union of supports,
then KL is computed by numerical integration. Returns 0.0 when one or both
KDEs are degenerate; returns Inf if any p > 0 region has q ≈ 0.
"""
function kl_divergence(samples_p::AbstractVector{<:Real},
                       samples_q::AbstractVector{<:Real};
                       n_grid::Int=200, q_floor::Float64=1e-12)
    p_kde = KDEPrior(collect(Float64, samples_p))
    q_kde = KDEPrior(collect(Float64, samples_q))
    lo = min(minimum(samples_p), minimum(samples_q))
    hi = max(maximum(samples_p), maximum(samples_q))
    pad = 0.05 * (hi - lo + 1e-12)
    grid = range(lo - pad, hi + pad; length=n_grid)
    dx = step(grid)
    kl = 0.0
    for x in grid
        p = pdf(p_kde, x)
        q = max(pdf(q_kde, x), q_floor)  # floor q to keep KL finite
        if p > q_floor
            kl += p * (log(p) - log(q)) * dx
        end
    end
    return max(kl, 0.0)
end

"""
    chain_kl(chain_new, chain_old, param_names)

Per-parameter KL divergences between two MCMCChains, summed over `param_names`.
"""
function chain_kl(chain_new, chain_old, param_names::Vector{Symbol})
    new_names = Set(String.(names(chain_new, :parameters)))
    old_names = Set(String.(names(chain_old, :parameters)))
    total = 0.0
    for name in param_names
        s = string(name)
        (s in new_names && s in old_names) || continue
        sp = vec(chain_new[s].data)
        sq = vec(chain_old[s].data)
        total += kl_divergence(sp, sq)
    end
    return total
end

"""
    iterative_infer(ode_models, ode_data, ssa_models, ssa_data; kwargs...)

Level-2 boundary protocol: iterate ODE↔SSA message-passing until the ODE
posterior stops changing (by summed KL across shared parameters) or `max_iters`
is reached.

At each iteration k ≥ 1:
1. NUTS on the ODE block, with priors built from the previous SSA posterior
   (default priors at k=1).
2. ABC-SMC on the SSA block, with priors built from this iteration's ODE chain.
3. Compute KL(ode_chain_k, ode_chain_{k-1}); stop if < `kl_tol`.

Returns `(ode_chain, ssa_posterior, n_iters, kl_trace, ode_history, ssa_history, boundary)`.
The history vectors are empty when `keep_history=false` (default; set true for
the integration test).
"""
function iterative_infer(ode_models::Vector{<:AbstractSubModel}, ode_data::ObservedData,
                         ssa_models::Vector{<:AbstractSubModel}, ssa_data::ObservedData;
                         ode_kwargs=(;), ssa_kwargs=(;),
                         boundary_method::Symbol=:kde,
                         max_iters::Int=5, kl_tol::Float64=0.05,
                         n_backward_samples::Int=1000,
                         keep_history::Bool=false,
                         verbose::Bool=false,
                         rng::Random.AbstractRNG=Random.default_rng())
    kl_trace = Float64[]
    ode_history = []
    ssa_history = []

    shared_names = let
        ode_free = unique_params(model_free_params(reduce(vcat, parameters.(ode_models))))
        ssa_free = unique_params(model_free_params(reduce(vcat, parameters.(ssa_models))))
        ssa_set = Set(p.name for p in ssa_free)
        Symbol[p.name for p in ode_free if p.name in ssa_set]
    end
    verbose && println("Shared parameters across blocks: $shared_names")

    ssa_posterior = nothing
    ode_chain = nothing
    prev_ode_chain = nothing
    bc = nothing

    for iter in 1:max_iters
        verbose && println("\n=== Iteration $iter/$max_iters ===")

        # Forward direction: SSA posterior -> ODE priors
        ode_priors_override = nothing
        if ssa_posterior !== nothing
            ode_bc = boundary_condition(ssa_posterior, ode_models;
                                        method=boundary_method,
                                        n_samples=n_backward_samples, rng=rng)
            ode_priors_override = Dict{Symbol, Distribution}()
            for (name, prior, cond) in zip(ode_bc.param_names, ode_bc.priors, ode_bc.conditioned)
                cond && (ode_priors_override[name] = prior)
            end
            verbose && println("ODE conditioned on $(length(ode_priors_override)) params from SSA posterior")
        end

        verbose && println("Stage 1: ODE inference (NUTS)...")
        ode_chain = infer(ode_models, ode_data;
                          priors_override=ode_priors_override, ode_kwargs...)

        # Backward direction: ODE chain -> SSA priors
        verbose && println("Stage 2: Boundary conditioning (ODE -> SSA)...")
        bc = boundary_condition(ode_chain, ssa_models; method=boundary_method)
        verbose && println("  Conditioned $(sum(bc.conditioned)) / $(length(bc.param_names)) SSA params")

        verbose && println("Stage 3: SSA inference (ABC-SMC)...")
        merged_ssa = merge((verbose=verbose,), ssa_kwargs,
                           (priors=bc.priors, param_names=bc.param_names))
        ssa_posterior = _infer_abc(ssa_models, ssa_data; merged_ssa...)

        # Convergence check
        if prev_ode_chain !== nothing
            kl = chain_kl(ode_chain, prev_ode_chain, shared_names)
            push!(kl_trace, kl)
            verbose && println("KL(ode_iter$iter || ode_iter$(iter-1)) = $(round(kl; digits=4))")
            if kl < kl_tol
                verbose && println("Converged (kl=$kl < tol=$kl_tol).")
                keep_history && (push!(ode_history, ode_chain); push!(ssa_history, ssa_posterior))
                return (ode_chain=ode_chain, ssa_posterior=ssa_posterior,
                        n_iters=iter, kl_trace=kl_trace,
                        ode_history=ode_history, ssa_history=ssa_history, boundary=bc)
            end
        end

        prev_ode_chain = ode_chain
        if keep_history
            push!(ode_history, ode_chain)
            push!(ssa_history, ssa_posterior)
        end
    end

    verbose && println("Reached max_iters=$max_iters without converging.")
    return (ode_chain=ode_chain, ssa_posterior=ssa_posterior,
            n_iters=max_iters, kl_trace=kl_trace,
            ode_history=ode_history, ssa_history=ssa_history, boundary=bc)
end

iterative_infer(ode_model::AbstractSubModel, ode_data::ObservedData,
                ssa_model::AbstractSubModel, ssa_data::ObservedData; kwargs...) =
    iterative_infer([ode_model], ode_data, [ssa_model], ssa_data; kwargs...)
