struct ABCPosterior
    particles::Matrix{Float64}   # n_params x n_particles
    weights::Vector{Float64}
    param_names::Vector{Symbol}
    tolerance::Float64
    n_populations::Int
end

function abc_smc(simulate, observed_stats::Vector{Float64},
                 priors::Vector{<:Distribution}, param_names::Vector{Symbol};
                 n_particles=1000, n_populations=10, alpha=0.5,
                 verbose=false, rng=Random.default_rng())
    n_params = length(priors)
    particles = zeros(n_params, n_particles)
    weights = fill(1.0 / n_particles, n_particles)
    distances = zeros(n_particles)

    # Population 1: sample from prior, accept all
    verbose && (print("ABC-SMC population 1/$n_populations (prior)..."); flush(stdout))
    for i in 1:n_particles
        for j in 1:n_params
            particles[j, i] = rand(rng, priors[j])
        end
        sim_stats = simulate(particles[:, i])
        distances[i] = summary_distance(sim_stats, observed_stats)
    end

    # Initial tolerance: alpha-quantile of distances
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

        # Compute perturbation kernel widths (twice the weighted std of each parameter)
        kernel_widths = _compute_kernel_widths(prev_particles, prev_weights)

        new_particles = zeros(n_params, n_particles)
        new_distances = zeros(n_particles)
        n_simulations = 0

        for i in 1:n_particles
            accepted = false
            while !accepted
                # Pick a particle from previous population (weighted)
                idx = _weighted_sample(prev_weights, rng)
                # Perturb
                candidate = _perturb(prev_particles[:, idx], kernel_widths, rng)
                # Check prior support
                all(pdf(priors[j], candidate[j]) > 0 for j in 1:n_params) || continue
                # Simulate and check distance
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

        # Update weights
        for i in 1:n_particles
            prior_prob = prod(pdf(priors[j], new_particles[j, i]) for j in 1:n_params)
            kernel_sum = sum(
                prev_weights[k] * _kernel_density(new_particles[:, i],
                    prev_particles[:, k], kernel_widths)
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
    widths = zeros(n_params)
    for j in 1:n_params
        wmean = sum(weights .* particles[j, :])
        wvar = sum(weights .* (particles[j, :] .- wmean).^2)
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

function _perturb(particle::Vector{Float64}, widths::Vector{Float64}, rng)
    return particle .+ widths .* randn(rng, length(particle))
end

function _kernel_density(x::AbstractVector, mu::AbstractVector, widths::Vector{Float64})
    d = length(x)
    exponent = -0.5 * sum(((x[j] - mu[j]) / widths[j])^2 for j in 1:d)
    norm_const = prod((sqrt(2π) * widths[j]) for j in 1:d)
    return exp(exponent) / norm_const
end

function quantile(v::AbstractVector, p::Real)
    sorted = sort(v)
    idx = clamp(ceil(Int, p * length(sorted)), 1, length(sorted))
    return sorted[idx]
end
