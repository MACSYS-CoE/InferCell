# Issue #89 B7: do reference_theta's added stops move θ_ref? It calls
# reference_theta with its defaults at the commit it runs from: the relative-change
# stop at 6526b45 (§12 2026-10-25), and that plus the x-step stop from 8da6045
# (2026-10-26). Reruns one AV6 prep's MAP under it and compares x with the θ_ref the
# production run stored (when that prep finished under Optim's tests alone), in
# prior SDs. For a prep that never finished, reports only the time and method.
#
#   julia --project dev/scripts/av6_89_b7_mapcheck.jl <replicate r> <run dir>

using InferCell
using Serialization
using Printf

function main(r, run)
    ds = m0_dataset(16_200_000 + r; ncells = 50, horizon = M0_HORIZON_S)   # as the prep
    s = mixed_setup(ds)
    t = @elapsed ref = reference_theta(s; maxiters = 200)
    @printf("r%d: %s in %.2f h, log density %.6f\n", r, ref.method, t / 3600, ref.logdensity)
    file = joinpath(run, "r$r", "prep.jls")
    if isfile(file)
        p = deserialize(file)
        sd = mixed_prior_sd(s)
        z = (ref.x .- p.ref.x) ./ sd
        @printf("r%d: stored θ_ref (%s, %.2f h, log density %.6f); max |Δx| %.3g prior SD, Δ log density %.3g\n",
                r, p.ref.method, p.seconds.map / 3600, p.ref.logdensity, maximum(abs.(z)),
                ref.logdensity - p.ref.logdensity)
        println("r$r: Δx / prior SD = ", join((@sprintf("%+.2e", v) for v in z), " "))
    end
    println("r$r: x = ", join((@sprintf("%.6f", v) for v in ref.x), " "))
end

main(parse(Int, ARGS[1]), ARGS[2])
