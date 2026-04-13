function compute_summary_stats(trajectories::Vector, species::Vector{Symbol};
                                times=nothing)
    n_traj = length(trajectories)
    n_species = length(species)

    if times === nothing
        final_vals = zeros(n_species, n_traj)
        for (j, sol) in enumerate(trajectories)
            for i in 1:n_species
                final_vals[i, j] = sol[i, end]
            end
        end
        means = vec(sum(final_vals; dims=2) ./ n_traj)
        vars  = vec(sum((final_vals .- means).^2; dims=2) ./ (n_traj - 1))
        fanos = [m > 0 ? v / m : 0.0 for (m, v) in zip(means, vars)]
        return vcat(means, vars, fanos)
    else
        stats = Float64[]
        sizehint!(stats, n_species * length(times))
        vals = zeros(n_species, n_traj)
        for t_idx in eachindex(times)
            fill!(vals, 0.0)
            for (j, sol) in enumerate(trajectories)
                snapshot = sol(times[t_idx])
                for i in 1:n_species
                    vals[i, j] = snapshot[i]
                end
            end
            means = vec(sum(vals; dims=2) ./ n_traj)
            append!(stats, means)
        end
        return stats
    end
end

function summary_distance(s1::Vector{Float64}, s2::Vector{Float64})
    d = 0.0
    @inbounds for i in eachindex(s1, s2)
        diff = s1[i] - s2[i]
        d += diff * diff
    end
    return sqrt(d)
end
