# Spec §11 tasks 4.5 and 4.6: what the 60 s rebuild's channel gain is, and what
# the drain granularity costs.
#
# Two measurements on the phase-3/4 toy, both of which phase 10 and phase 13
# inherit as numbers rather than as assumptions:
#
#   4.5  The elasticity of a rebuilt rate constant to its upstream pool,
#        d ln k / d ln pool, reported *as a function of* the quantities we
#        assert rather than as a point value (spec §10 R8).
#
#   4.6  Spec §4 D10's granularity comparison, in miniature: the nominal
#        trajectory under 1 s, 5 s and 60 s drain granularity, and the largest
#        relative difference in any pool.
#
# **On why 4.6 is an ensemble.** The three configurations cannot be compared on
# one seeded path. Coarsening the drain changes how often the hook clears the
# counters, and the driver rebuilds the SSA's propensity aggregation whenever it
# does — which consumes randomness — so the same seed gives three different
# stochastic paths. A single-path difference would therefore be Monte Carlo
# noise wearing the label of a granularity effect. Averaging over seeds and
# quoting the standard error of the mean beside the difference is what makes the
# number attributable, and it is the method phase 13 needs at full scale.
#
# Usage: julia --project dev/scripts/rebuild_channel.jl

using InferCell
using Random
using Printf
using Dates
using Statistics: mean, std

include(joinpath(@__DIR__, "..", "..", "test", "contribution_test_models.jl"))
include(joinpath(@__DIR__, "..", "..", "test", "jump_test_models.jl"))
include(joinpath(@__DIR__, "..", "..", "test", "hybrid_test_models.jl"))

const HORIZON = 600            # s
const SEEDS = 300              # replicates per granularity
const D10_THRESHOLD = 0.01     # spec §4 D10: one percent

# Provenance captured before the measurements, not after: on a long job HEAD can
# move between load and write, and the artefact would then name source the
# measurement never saw.
const COMMIT = strip(read(`git rev-parse --short HEAD`, String))
const DIRTY = !isempty(strip(read(`git status --porcelain -- src test Project.toml Manifest.toml dev/scripts/rebuild_channel.jl`, String)))

# A metabolic half that moves without running dry, so the pool level genuinely
# feeds back into the flux and a drain that arrives late has something to change.
slow_pool(; kwargs...) = ToyPool(; kcat = 0.02, atp0 = 20.0, kwargs...)

# ---------------------------------------------------------------------------
# 4.5 — the channel's gain
# ---------------------------------------------------------------------------

gain_rows = NamedTuple[]
for km in (2.0, 1.0, 0.5, 0.175, 0.05), atp in (20.0, 3.6529, 1.0)
    m = ToyRebuiltExpression(km_tx = km)
    d = build_problem([slow_pool(atp0 = atp), m]; tspan = (0.0, 60.0))
    row = only(rate_constant_elasticity(d))
    analytic = toy_rebuilt_elasticity(m, atp)
    push!(gain_rows, (km = km, atp = atp, k = row.rate_constant,
                      measured = row.elasticity, analytic = analytic,
                      err = abs(row.elasticity - analytic) / analytic))
end
gain_ok = maximum(r -> r.err, gain_rows) < 1e-5

# ---------------------------------------------------------------------------
# 4.6 — the drain granularity
# ---------------------------------------------------------------------------

# One ensemble per granularity: the mean pool trajectory over `SEEDS` replicates,
# and the standard error of that mean. No rebuild here — the question is what the
# *drain* costs, so the stochastic block is left independent of the pools.
function ensemble(drain)
    sums = zeros(Float64, HORIZON, 2)
    sqs = zeros(Float64, HORIZON, 2)
    for s in 1:SEEDS
        Random.seed!(90_000 + s)
        d = build_problem([slow_pool(), ToyExpression()];
                          tspan = (0.0, Float64(HORIZON)), drain_interval = drain)
        out = run_handshake!(d, HORIZON)
        for i in 1:HORIZON, j in 1:2
            sums[i, j] += out.ode[i][j]
            sqs[i, j] += out.ode[i][j]^2
        end
    end
    m = sums ./ SEEDS
    v = max.(sqs ./ SEEDS .- m .^ 2, 0.0)
    return (mean = m, sem = sqrt.(v ./ SEEDS))
end

base = ensemble(1.0)
drain_rows = NamedTuple[]
for drain in (5.0, 60.0)
    e = ensemble(drain)
    rel = abs.(e.mean .- base.mean) ./ max.(abs.(base.mean), eps())
    # The noise floor the difference has to clear: the standard error of the
    # difference of two independent means, in the same relative units.
    noise = sqrt.(e.sem .^ 2 .+ base.sem .^ 2) ./ max.(abs.(base.mean), eps())
    k = argmax(rel)
    push!(drain_rows, (drain = drain, max_rel = rel[k], at_t = k[1],
                       pool = k[2] == 1 ? :M_atp_c : :M_adp_c,
                       noise = noise[k], ratio = rel[k] / noise[k]))
end

verdict(r) = r.max_rel < D10_THRESHOLD ? "below 1%" : "above 1%"

# ---------------------------------------------------------------------------
# Report
# ---------------------------------------------------------------------------

jobid = get(ENV, "SLURM_JOB_ID", "not under Slurm")
io = IOBuffer()
println(io, "# The 60 s rebuild: channel gain and drain granularity (spec tasks 4.5, 4.6)")
println(io)
println(io, "Generated by `dev/scripts/rebuild_channel.jl` at commit $COMMIT",
        DIRTY ? " (WARNING: tree dirty)" : "",
        ", Slurm job $jobid, on $(gethostname()), $(Dates.now()), Julia $(VERSION).")
println(io)
println(io, "Regenerate with `sbatch dev/scripts/rebuild_channel.slurm`.")
println(io)
println(io, "## 4.5 — the elasticity of a rebuilt rate constant to its pool")
println(io)
println(io, "`ToyRebuiltExpression` rebuilds `k_tx_rb = k_max · atp / (km + atp)`, whose")
println(io, "elasticity `d ln k / d ln atp` is `km / (km + atp)` in closed form. The")
println(io, "diagnostic is a central difference in log space, so the analytic column is a")
println(io, "check on the measurement rather than a restatement of it.")
println(io)
println(io, "| km (mM) | pool (mM) | k (1/s) | measured ε | analytic ε | rel. error |")
println(io, "|---|---|---|---|---|---|")
for r in gain_rows
    @printf(io, "| %.3f | %.4f | %.4f | %.6f | %.6f | %.2e |\n",
            r.km, r.atp, r.k, r.measured, r.analytic, r.err)
end
println(io)
println(io, gain_ok ?
        "Every row agrees with the closed form to better than 1e-5 relative." :
        "**WARNING: a row disagrees with the closed form by more than 1e-5.**")
println(io)
println(io, """
**The gain is a function of the pool, not a constant of the model**, which is
what spec §10 R8 requires of every channel gain: the same rebuild law reports an
elasticity spanning more than an order of magnitude across the pool range above.
The toy's law is ours and its `km` is arbitrary; the number that matters is the
*real* transcription channel's, which `dev/notes/reduced-syn3a-scoping.md`
measures at 0.044–0.051 and which phase 10 computes from the published rate law.
What this table establishes is that the diagnostic recovers a known elasticity
exactly, so the phase-10 number will be the model's and not the estimator's.""")
println(io)
println(io, "## 4.6 — spec §4 D10's granularity comparison, in miniature")
println(io)
@printf(io, "%d replicates per granularity, %d s horizon, drain compared against the\n",
        SEEDS, HORIZON)
println(io, "published 1 s drain. `max rel.` is the largest relative difference in the")
println(io, "*mean* trajectory of any pool at any handshake; `noise` is the standard")
println(io, "error of that difference in the same units, so `ratio` below about 2 means")
println(io, "the difference is not resolved above Monte Carlo noise.")
println(io)
println(io, "| drain (s) | max rel. | at t (s) | pool | noise | ratio | vs D10's 1% |")
println(io, "|---|---|---|---|---|---|---|")
for r in drain_rows
    @printf(io, "| %.0f | %.4f | %d | %s | %.4f | %.1f | %s |\n",
            r.drain, r.max_rel, r.at_t, r.pool, r.noise, r.ratio, verdict(r))
end
println(io)
println(io, """
**Read narrowly.** This is the toy, not Core A′. Its ATP pool is ~404,000
particles against a drain of order 40 particles per second, a buffer of hours,
while the pools D10's comparison is really about turn over in 109 s (adenylate)
and 30 s (guanylate) — which is why D10 expects 60 s to fail there and 5 s to
pass. A small difference here is evidence that the *mechanism* aggregates and
pays correctly, and a bound on nothing else. The measurement that decides the
granularity is phase 13's, on the assembled model; what this phase delivers is
the knob, the labelling and the method — including that the comparison has to be
an ensemble, because coarsening the drain changes how often the propensity
aggregation is rebuilt and therefore the seeded path.""")

s = String(take!(io))
print(s)
write(joinpath(@__DIR__, "rebuild_channel_result.md"), s)

gain_ok || error("The elasticity diagnostic disagrees with its closed form; see the table above")
