# Export V7's merged ranks (16a.9b) to TSVs that F8's plots read
# (plot_v7_16a9b.py), so the plots use the merge's own ranks and bands.
#   julia dev/scripts/v7_16a9b_export.jl <dir holding ranks.jls>
using Serialization, Printf

function export_ranks(dir)
    m = deserialize(joinpath(dir, "ranks.jls"))
    open(joinpath(@__DIR__, "v7_16a9b_ranks.tsv"), "w") do io
        println(io, join(["replicate"; "ess_min"; m.names], '\t'))
        for j in 1:m.R
            @printf(io, "%d\t%.1f\t%s\n", m.replicates[j], m.ess_min[j], join(m.ranks[j, :], '\t'))
        end
    end
    open(joinpath(@__DIR__, "v7_16a9b_bands.tsv"), "w") do io
        println(io, "# R = $(m.R), L = $(m.L), α = $(m.α) per parameter, pointwise γ = $(m.γ)")
        println(io, "z\tlo\thi")
        for i in eachindex(m.z)
            @printf(io, "%.17g\t%d\t%d\n", m.z[i], m.lo[i], m.hi[i])
        end
    end
end

export_ranks(ARGS[1])
