# Analyse the 14c.5 decomposition (corea_smoothing_14c5_decomp.jl). For each
# seed and each variant against the published model: the handshake at which the
# pair decouples, and the GTP gap at the seed's worst coupled save point, both
# plain and net of the carried deficits. Then the published pool and the gaps
# every 10 s over the two minutes before that save point.
# Usage: julia --project dev/scripts/corea_smoothing_14c5_decomp_merge.jl <dir of .jls>
using Serialization, Printf
worst = Dict(18227 => 4500, 16184 => 5880, 19958 => 5460, 16391 => 3300)
for f in sort(filter(x -> startswith(x, "decomp_") && endswith(x, ".jls"), readdir(ARGS[1])))
    r = deserialize(joinpath(ARGS[1], f)); A = r.res[:A]; tw = worst[r.seed]
    println("=== seed $(r.seed), worst save t=$tw, GTP counters $(A.labels)")
    for v in (:B, :C, :D)
        x = r.res[v]
        dec = findfirst(k -> x.jumps[k] != A.jumps[k], 1:r.T)
        g = x.st[:M_gtp_c] .- A.st[:M_gtp_c]
        # pool minus total carried deficit on GTP: identical drains would keep this equal
        net = (x.st[:M_gtp_c] .- vec(sum(x.def; dims = 2))) .- (A.st[:M_gtp_c] .- vec(sum(A.def; dims = 2)))
        @printf("%s: decouples at %s; GTP gap at t=%d: %.2f (A pool %.1f); net(pool-deficit) gap %.2f; max|gap| before decoupling %.2f\n",
                v, dec === nothing ? "never" : string(dec), tw, g[tw], A.st[:M_gtp_c][tw], net[tw],
                maximum(abs, g[1:(dec === nothing ? r.T : dec - 1)]))
    end
    D = r.res[:D]
    println("  t    A.gtp   D-A gtp  B-A gtp  C-A gtp  D-A net  A.def(sum) D.def(sum) A.clip D.clip  D-A gdp  D-A gmp  D-A atp")
    for k in (tw - 120):10:tw
        dA = sum(A.def[k, :]); dD = sum(D.def[k, :])
        @printf("%5d %8.2f %8.3f %8.3f %8.3f %8.3f %9.3g %9.3g %6s %6s %8.2f %8.2f %8.2f\n", k, A.st[:M_gtp_c][k],
                D.st[:M_gtp_c][k] - A.st[:M_gtp_c][k], r.res[:B].st[:M_gtp_c][k] - A.st[:M_gtp_c][k],
                r.res[:C].st[:M_gtp_c][k] - A.st[:M_gtp_c][k],
                (D.st[:M_gtp_c][k] - dD) - (A.st[:M_gtp_c][k] - dA), dA, dD, any(A.clip[k, :]), any(D.clip[k, :]),
                D.st[:M_gdp_c][k] - A.st[:M_gdp_c][k], D.st[:M_gmp_c][k] - A.st[:M_gmp_c][k],
                D.st[:M_atp_c][k] - A.st[:M_atp_c][k])
    end
end
