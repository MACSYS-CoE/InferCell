# Issue #89 B4: route (ii) against the prototype extrapolated to Δt → 0 (§12 2026-10-23,
# the amended check). The prototype's error is first order in its step h, so
# 2·proto(h/2) − proto(h) removes it, and route − that = 2·Δ(h/2) − Δ(h) per pool and save.
# The prototype's limit still reads the pools at 60m, so the difference includes
# route (ii)'s declared 1 s snapshot shift.
#   julia --project dev/scripts/route_ii_89_richardson.jl <outdir of route_ii_89.sh>
using InferCell, Serialization, Printf

const SIGMA_B = 0.10
const TOL_KG_FRACTION = 0.1
const TOL_MAX_LOG = 0.02
const ORDER = ["m0_v6", "m0_low", "m0_high", "corea1", "corea2"]
const B3_MERGED = "/fred/oz022/tkimpson/infercell-89-b3-archive-2026-10-20/merged.jls"

kg(E) = sum(E .^ 2) / (2SIGMA_B^2)

function main(out)
    d = deserialize(joinpath(out, "merged.jls"))
    b3 = deserialize(B3_MERGED)
    commit = readchomp(`git -C $(@__DIR__) rev-parse --short HEAD`)
    println("Extrapolated at $commit (job $(get(ENV, "SLURM_JOB_ID", "none"))), from $(out)/merged.jls.")
    println()
    println("| Truth | pair (s) | KG | bound | max \\|Δ\\| (pool) | verdict | 3PG's Δ at 1 / 0.5 / 0.25 s |")
    println("|---|---|---|---|---|---|---|")
    n = 0
    for t in ORDER
        D = d[t].D
        # A missing step would silently move the extrapolation to a coarser pair.
        all(haskey(D, s) for s in (1.0, 0.5, 0.25)) || error("$t: a prototype step is missing")
        h = (0.5, 0.25)
        E = 2 .* D[h[2]] .- D[h[1]]
        bound = TOL_KG_FRACTION * sum(b3[t].shift)
        mx, i = findmax(vec(maximum(abs.(E); dims = 2)))
        ok = kg(E) <= bound && mx <= TOL_MAX_LOG
        n += ok
        j = findfirst(==(:M_3pg_c), M0_PANEL)
        tpg = join([haskey(D, s) ? @sprintf("%.4f", maximum(abs.(D[s][j, :]))) : "—" for s in (1.0, 0.5, 0.25)], " / ")
        @printf("| %s | %g, %g | %.3g | %.3g | %.4f (`%s`) | %s | %s |\n", t, h..., kg(E), bound, mx,
                M0_PANEL[i], ok ? "pass" : "**FAIL**", tpg)
    end
    println()
    println("**Amended verdict:** pass at $n of $(length(ORDER)) truths.")
end

main(ARGS[1])
