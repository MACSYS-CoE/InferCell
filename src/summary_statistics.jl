function compute_summary_stats(trajectories::Vector, species::Vector{Symbol};
                                times=nothing)
    n_traj = length(trajectories)
    n_species = length(species)

    if times === nothing
        # Use the last timepoint (steady state)
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
        # Summary stats at each observed timepoint
        stats = Float64[]
        for t_idx in eachindex(times)
            vals = zeros(n_species, n_traj)
            for (j, sol) in enumerate(trajectories)
                for i in 1:n_species
                    vals[i, j] = sol(times[t_idx])[i]
                end
            end
            means = vec(sum(vals; dims=2) ./ n_traj)
            append!(stats, means)
        end
        return stats
    end
end

function summary_distance(s1::Vector{Float64}, s2::Vector{Float64})
    return norm(s1 .- s2)
end
