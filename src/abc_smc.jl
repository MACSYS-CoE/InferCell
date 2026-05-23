"""
    ABCPosterior(particles, weights, param_names, tolerance, n_populations)

Result of [`abc_smc`](@ref): a weighted particle approximation of the
posterior. `particles` is `n_params × n_particles`, `weights` sums to one, and
`param_names` orders the rows of `particles`. `tolerance` and `n_populations`
record the final acceptance threshold and the number of SMC populations run.
"""
struct ABCPosterior
    particles::Matrix{Float64}   # n_params x n_particles
    weights::Vector{Float64}
    param_names::Vector{Symbol}
    tolerance::Float64
    n_populations::Int
end

"""
    abc_smc(simulate, observed_stats, priors, param_names; n_particles=1000, n_populations=10, alpha=0.5, verbose=false, rng=...)

Approximate Bayesian Computation with Sequential Monte Carlo (Beaumont et al.).
`simulate(theta)` must return a summary-statistic vector comparable to
`observed_stats`. The tolerance schedule is adaptive: at each population the
new tolerance is the `alpha`-quantile of the previous distances.

Returns an [`ABCPosterior`](@ref).
"""
function abc_smc(simulate, observed_stats::Vector{Float64},
                 priors::Vector{<:Distribution}, param_names::Vector{Symbol};
                 n_particles=1000, n_populations=10, alpha=0.5,
                 verbose=false, rng=Random.default_rng())
    n_params = length(priors)
    particles = zeros(n_params, n_particles)
    weights = fill(1.0 / n_particles, n_particles)
    distances = zeros(n_particles)

    verbose && (print("ABC-SMC population 1/$n_populations (prior)..."); flush(stdout))
    for i in 1:n_particles
        for j in 1:n_params
            particles[j, i] = rand(rng, priors[j])
        end
        sim_stats = simulate(particles[:, i])
        distances[i] = summary_distance(sim_stats, observed_stats)
    end

    tolerance = quantile(distances, alpha)
    if verbose
        println(" tolerance=$(round(tolerance; digits=3))")
        flush(stdout)
    end

    for pop in 2:n_populations
        verbose && (print("ABC-SMC population $pop/$n_populations..."); flush(stdout))
        t_pop = time()

        prev_particles = copy(particles)
        prev_weights = copy(weights)
        kernel_widths = _compute_kernel_widths(prev_particles, prev_weights)

        new_particles = zeros(n_params, n_particles)
        new_distances = zeros(n_particles)
        n_simulations = 0

        candidate = zeros(n_params)
        for i in 1:n_particles
            accepted = false
            while !accepted
                idx = _weighted_sample(prev_weights, rng)
                _perturb!(candidate, @view(prev_particles[:, idx]), kernel_widths, rng)
                all(pdf(priors[j], candidate[j]) > 0 for j in 1:n_params) || continue
                sim_stats = simulate(candidate)
                n_simulations += 1
                d = summary_distance(sim_stats, observed_stats)
                if d <= tolerance
                    new_particles[:, i] = candidate
                    new_distances[i] = d
                    accepted = true
                end
            end
        end

        for i in 1:n_particles
            prior_prob = prod(pdf(priors[j], new_particles[j, i]) for j in 1:n_params)
            kernel_sum = sum(
                prev_weights[k] * _kernel_density(
                    @view(new_particles[:, i]),
                    @view(prev_particles[:, k]), kernel_widths)
                for k in 1:n_particles
            )
            weights[i] = prior_prob / kernel_sum
        end
        weights ./= sum(weights)

        particles = new_particles
        tolerance = quantile(new_distances, alpha)

        if verbose
            elapsed = round(time() - t_pop; digits=1)
            acc_rate = round(n_particles / n_simulations * 100; digits=1)
            println(" tolerance=$(round(tolerance; digits=3)), acc_rate=$(acc_rate)%, sims=$n_simulations, $(elapsed)s")
            flush(stdout)
        end
    end

    return ABCPosterior(particles, weights, param_names, tolerance, n_populations)
end

function _compute_kernel_widths(particles::Matrix{Float64}, weights::Vector{Float64})
    n_params = size(particles, 1)
    n_particles = size(particles, 2)
    widths = zeros(n_params)
    for j in 1:n_params
        wmean = 0.0
        for i in 1:n_particles
            wmean += weights[i] * particles[j, i]
        end
        wvar = 0.0
        for i in 1:n_particles
            wvar += weights[i] * (particles[j, i] - wmean)^2
        end
        widths[j] = 2.0 * sqrt(wvar)
    end
    return widths
end

function _weighted_sample(weights::Vector{Float64}, rng)
    r = rand(rng)
    cumw = 0.0
    for i in eachindex(weights)
        cumw += weights[i]
        cumw >= r && return i
    end
    return length(weights)
end

function _perturb!(out::Vector{Float64}, particle::AbstractVector{Float64},
                    widths::Vector{Float64}, rng)
    for j in eachindex(out)
        out[j] = particle[j] + widths[j] * randn(rng)
    end
    return out
end

function _kernel_density(x::AbstractVector, mu::AbstractVector, widths::Vector{Float64})
    d = length(x)
    exponent = -0.5 * sum(((x[j] - mu[j]) / widths[j])^2 for j in 1:d)
    norm_const = prod((sqrt(2π) * widths[j]) for j in 1:d)
    return exp(exponent) / norm_const
end

