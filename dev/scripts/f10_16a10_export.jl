# Export 16a.10's draws to a TSV that F10's analysis reads (f10_16a10.py): V6's
# four reference chains and the production chain 5, every sweep, on the log scale.
# It refuses unless chain 5 sits on V6's dataset (spec/phases/16-recovery.md §12
# 2026-10-13 and 2026-10-22).
#   julia --project dev/scripts/f10_16a10_export.jl <V6 dir> <F10 dir>
# The V6 dir holds chain1.jls to chain4.jls; the F10 dir holds chain5.jls.
using InferCell
using Serialization
using Printf

function export_draws(v6dir, f10dir)
    files = [joinpath(v6dir, "chain$(k).jls") for k in 1:4]
    push!(files, joinpath(f10dir, "chain5.jls"))
    saved = [deserialize(f) for f in files]
    ref = first(saved).meta
    for (k, r) in enumerate(saved), f in (:replicate, :ncells, :N, :truth, :sigma_true, :names)
        getfield(r.meta, f) == getfield(ref, f) ||
            error("chain $k differs from chain 1 in $f: not the same dataset or setup")
    end
    names = string.(ref.names)
    out = joinpath(@__DIR__, "f10_16a10_draws.tsv")
    open(out, "w") do io
        println(io, join(["chain"; "sweep"; names], '\t'))
        for (k, r) in enumerate(saved), (i, t) in enumerate(r.state.trace)
            x = vcat(log.(t.cme), t.u, log(t.σ))
            println(io, k, '\t', i, '\t', join((@sprintf("%.17g", v) for v in x), '\t'))
        end
    end
    truth = Dict(ref.truth)
    tv = vcat([log(truth[Symbol(n)]) for n in names[1:end - 1]], log(ref.sigma_true))
    open(joinpath(@__DIR__, "f10_16a10_truth.tsv"), "w") do io
        println(io, join(names, '\t'))
        println(io, join((@sprintf("%.17g", v) for v in tv), '\t'))
    end
    println("wrote $out: chains of ", join([length(r.state.trace) for r in saved], ", "), " sweeps")
end

export_draws(ARGS[1], ARGS[2])
