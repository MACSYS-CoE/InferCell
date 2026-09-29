"""
    compute_summary_stats(trajectories, species; times=nothing, spread=true,
                          state_names=species)

Compress an ensemble of SSA trajectories into a summary-statistic vector for
ABC-SMC. `species` are looked up by name in `state_names`, the trajectories'
state order, so any subset in any order may be summarised, and an unknown name
throws (spec §11 task 15.9).

When `times === nothing`, returns `[means; vars; fanos]` of the final state per
species. When `times` is provided, returns, for each time in turn, the species
means followed (when `spread`) by the species standard deviations across
replicates. The spread is what keeps the ensemble's low-copy information, which
per-time means alone discard (spec §10 R14, §12 2026-09-23 G). With one
trajectory the standard deviation is zero.
"""
function compute_summary_stats(trajectories::Vector, species::Vector{Symbol};
                               times=nothing, spread::Bool=true,
                               state_names::AbstractVector{Symbol}=species)
    n_traj = length(trajectories)
    n_species = length(species)
    rows = state_rows(species, state_names)

    if times === nothing
        final_vals = zeros(n_species, n_traj)
        for (j, sol) in enumerate(trajectories)
            for (i, r) in enumerate(rows)
                final_vals[i, j] = sol[r, end]
            end
        end
        means = vec(sum(final_vals; dims=2) ./ n_traj)
        vars  = vec(sum((final_vals .- means).^2; dims=2) ./ (n_traj - 1))
        fanos = [m > 0 ? v / m : 0.0 for (m, v) in zip(means, vars)]
        return vcat(means, vars, fanos)
    else
        stats = Float64[]
        sizehint!(stats, (spread ? 2 : 1) * n_species * length(times))
        vals = zeros(n_species, n_traj)
        for t_idx in eachindex(times)
            for (j, sol) in enumerate(trajectories)
                snapshot = sol(times[t_idx])
                for (i, r) in enumerate(rows)
                    vals[i, j] = snapshot[r]
                end
            end
            means = vec(sum(vals; dims=2) ./ n_traj)
            append!(stats, means)
            if spread
                sds = n_traj > 1 ?
                    sqrt.(vec(sum((vals .- means).^2; dims=2) ./ (n_traj - 1))) :
                    zeros(n_species)
                append!(stats, sds)
            end
        end
        return stats
    end
end

"""
    summary_distance(s1, s2) -> Float64

Euclidean distance between two summary-statistic vectors. The default distance
metric used by [`abc_smc`](@ref).
"""
function summary_distance(s1::Vector{Float64}, s2::Vector{Float64})
    d = 0.0
    @inbounds for i in eachindex(s1, s2)
        diff = s1[i] - s2[i]
        d += diff * diff
    end
    return sqrt(d)
end
