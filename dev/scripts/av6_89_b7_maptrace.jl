# Issue #89 B7: why 11 of AV6's 51 bounded MAPs ran past 12 h. Replays one prep's
# MAP (reference_theta's call: maximum_a_posteriori from the prior median, box
# bounds, maxiters = 200) on a model identical to mixed_model but for a wrapper that
# times and logs every likelihood evaluation, values unchanged. Also logs each
# iteration the optimiser reports. Stops after MAX_HOURS (default 6). RELTOL (default
# 0, the old call) is reference_theta's relative-change stop (§12 2026-10-25). It
# replays the call before xtol (§12 2026-10-26): reference_theta now also passes
# x_abstol and outer_x_abstol, which this replay leaves at Optim's 0.
#
#   julia --project dev/scripts/av6_89_b7_maptrace.jl <replicate r> <out.tsv>
#
# Reads: per evaluation, its wall time and whether it was a gradient (dual) or a
# value; per iteration, the objective. Many evaluations at the usual cost point at
# the optimiser (Fminbox maps maxiters to outer iterations, leaving 1,000 inner);
# a few very slow ones point at stiff trial points.

using InferCell
using InferCell: mixed_loglik
using Turing
using Turing: @addlogprob!
using ForwardDiff
using LinearAlgebra
using Printf

const R = parse(Int, ARGS[1])
const OUT = ARGS[2]
const MAX_HOURS = parse(Float64, get(ENV, "MAX_HOURS", "6"))
const RELTOL = parse(Float64, get(ENV, "RELTOL", "0"))
const T0 = time()
const IO_ = open(OUT, "w")
println(IO_, "kind\tn\tt_since_start\tseconds\tvalue\tx")
const N = Ref(0)
const NIT = Ref(0)

function traced_loglik(s, x)
    t = time()
    v = mixed_loglik(s, x)
    dt = time() - t
    N[] += 1
    dual = eltype(x) <: ForwardDiff.Dual
    xv = dual ? ForwardDiff.value.(x) : x
    vv = v isa ForwardDiff.Dual ? ForwardDiff.value(v) : v
    println(IO_, join((dual ? "grad" : "value", N[], @sprintf("%.1f", t - T0), @sprintf("%.3f", dt),
                       @sprintf("%.6g", vv), join((@sprintf("%.4f", a) for a in xv), ",")), '\t'))
    flush(IO_)
    return v
end

@model function traced_model(s)
    x ~ MvNormal(mixed_prior_median(s), Diagonal(vcat(s.prior_sd, s.sigma_prior.σ) .^ 2))
    @addlogprob! traced_loglik(s, x)
end

function callback(state, loss)
    NIT[] += 1
    println(IO_, join(("iter", NIT[], @sprintf("%.1f", time() - T0), "", @sprintf("%.6g", loss),
                       join((@sprintf("%.4f", a) for a in state.u), ",")), '\t'))
    flush(IO_)
    return time() - T0 > MAX_HOURS * 3600
end

function main()
    ds = m0_dataset(16_200_000 + R; ncells = 50, horizon = M0_HORIZON_S)   # as the prep
    s = mixed_setup(ds)
    lb, ub = mixed_box(s)
    x0 = mixed_prior_median(s)
    @printf("r%d: box ±%g SD; prior median %s\n", R, s.box_sd, x0)
    res = maximum_a_posteriori(traced_model(s); initial_params = collect(Float64, x0), lb, ub,
                               maxiters = 200, reltol = RELTOL, outer_f_reltol = RELTOL, callback)
    close(IO_)
    println(res.optim_result.original)
    @printf("r%d: %d evaluations, %d callback iterations, %.2f h; lp %.4f\n",
            R, N[], NIT[], (time() - T0) / 3600, res.lp)
end

main()
