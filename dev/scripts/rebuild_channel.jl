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
#   4.6  Spec §4 D10's granularity comparison, in miniature: what the nominal
#        trajectory costs under 1 s, 5 s and 60 s drain granularity.
#
# **Two things about 4.6 that a first attempt got wrong, recorded because phase
# 13 inherits the method and not only the number.**
#
# *The comparison must be sampled at drain-aligned instants.* Between drains a
# coarse configuration holds up to `drain - interval` seconds of accrued cost in
# its counters, so its pools sit high by exactly that unpaid debt. Sampling at
# an arbitrary handshake measures that sawtooth — a deterministic bookkeeping
# offset, `(drain - interval) * cost_rate / pool`, closed form and available
# without any simulation — and not whether coarsening the drain changed the
# dynamics. At a handshake that is a multiple of every drain being compared,
# every configuration has just debited the same total accrual, and what remains
# is the dynamical difference D10's one-percent rule is about. Both quantities
# are reported below, separately and labelled.
#
# *The comparison must be an ensemble, and paired.* Coarsening the drain changes
# how often the hook clears the counters, and the driver rebuilds the SSA's
# propensity aggregation whenever it does — which consumes randomness — so the
# same seed gives different stochastic paths. A single-path difference would be
# Monte Carlo noise wearing the label of a granularity effect. The two
# configurations do share seeds, so the difference is taken **per seed** and its
# standard error computed from those paired differences; combining two marginal
# standard errors as if the ensembles were independent would be wrong in an
# unknown direction.
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
const DRAINS = (5.0, 60.0)     # compared against the published 1 s drain
const ALIGN = 60               # sampling stride: a multiple of every drain above
const D10_THRESHOLD = 0.01     # spec §4 D10: one percent
const POOLS = (:M_atp_c, :M_adp_c)

# Provenance captured before the measurements, not after: on a long job HEAD can
# move between load and write, and the artefact would then name source the
# measurement never saw.
const COMMIT = strip(read(`git rev-parse --short HEAD`, String))
const DIRTY = !isempty(strip(read(`git status --porcelain -- src test Project.toml Manifest.toml dev/scripts/rebuild_channel.jl`, String)))

# `toy_slow_pool` lives in `test/hybrid_test_models.jl`, shared with the test
# suite so the two cannot drift apart on what "slow" means: a pool that moves at
# every handshake and ends around two thirds full, so the level genuinely feeds
# back into the flux and a drain that arrives late has something to change.

# ---------------------------------------------------------------------------
# 4.5 — the channel's gain
# ---------------------------------------------------------------------------

gain_rows = NamedTuple[]
for km in (2.0, 1.0, 0.5, 0.175, 0.05), atp in (20.0, 3.6529, 1.0)
    m = ToyRebuiltExpression(km_tx = km)
    d = build_problem([toy_slow_pool(atp0 = atp), m]; tspan = (0.0, 60.0))
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

const MODELS = () -> [toy_slow_pool(), ToyExpression()]

# Per seed, the two pools at every drain-aligned handshake, and — separately —
# the same pools one handshake *before* a drain, which is where the unpaid
# accrual is largest. No rebuild here: the question is what the *drain* costs,
# so the stochastic block is left independent of the pools.
function ensemble(drain)
    n = HORIZON ÷ ALIGN
    aligned = zeros(Float64, SEEDS, n, length(POOLS))
    predrain = zeros(Float64, SEEDS, n, length(POOLS))
    pending = zeros(Float64, SEEDS, n)
    for s in 1:SEEDS
        Random.seed!(90_000 + s)
        d = build_problem(MODELS(); tspan = (0.0, Float64(HORIZON)),
                          drain_interval = drain)
        out = run_handshake!(d, HORIZON)
        for i in 1:n, j in eachindex(POOLS)
            aligned[s, i, j] = out.ode[i * ALIGN][j]
            predrain[s, i, j] = out.ode[i * ALIGN - 1][j]
        end
        for i in 1:n
            pending[s, i] = out.jump[i * ALIGN - 1][3]   # :atp_cost, undebited
        end
    end
    return (aligned = aligned, predrain = predrain, pending = pending)
end

base = ensemble(1.0)
factor = corea_particles_per_mM()

# The paired relative difference at one family of instants: mean over seeds of
# (coarse - base) / mean(base), and the standard error of that mean taken from
# the per-seed differences themselves.
function paired(a, b)
    rel = similar(a)
    for i in axes(a, 2), j in axes(a, 3)
        denom = mean(@view b[:, i, j])
        rel[:, i, j] = (@view(a[:, i, j]) .- @view(b[:, i, j])) ./ denom
    end
    m = [mean(@view rel[:, i, j]) for i in axes(rel, 2), j in axes(rel, 3)]
    se = [std(@view rel[:, i, j]) / sqrt(SEEDS) for i in axes(rel, 2), j in axes(rel, 3)]
    return m, se
end

drain_rows = NamedTuple[]
for drain in DRAINS
    e = ensemble(drain)
    m, se = paired(e.aligned, base.aligned)
    n = size(m, 1)
    last_i = argmax(abs.(m[n, :]))                      # the final instant, no selection
    k = argmax(abs.(m))                                 # the largest of n × 2
    pm, _ = paired(e.predrain, base.predrain)
    sawtooth_i = argmax(abs.(pm))
    # The closed form the mid-period offset should equal: the accrual the coarse
    # configuration is still holding, as a fraction of the pool it has not left.
    held = mean(e.pending) - mean(base.pending)          # particles
    push!(drain_rows,
          (drain = drain,
           final = m[n, last_i], final_se = se[n, last_i], final_pool = POOLS[last_i],
           max = m[k], max_se = se[k], max_pool = POOLS[k[2]], max_t = k[1] * ALIGN,
           n_points = length(m),
           sawtooth = pm[sawtooth_i], sawtooth_pool = POOLS[sawtooth_i[2]],
           closed_form = held / factor / mean(@view base.predrain[:, end, 1])))
end

resolved(r) = abs(r.final) > 2 * r.final_se
verdict(r) = abs(r.final) < D10_THRESHOLD ? "below 1%" : "above 1%"

# The driver's own policy has to reach the report beside the numbers it produced
# (spec §6 T2): it is a departure the sub-models do not declare.
labels_60 = reduction_report(MODELS(),
                             build_problem(MODELS(); tspan = (0.0, 60.0),
                                           drain_interval = 60.0))

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
@printf(io, "%d replicates per granularity, %d s horizon. Differences are **paired per\n",
        SEEDS, HORIZON)
println(io, "seed** against the published 1 s drain and sampled at **drain-aligned**")
@printf(io, "handshakes (every %d s), where every configuration has just debited the same\n", ALIGN)
println(io, "total accrual — so what is left is the difference in the dynamics and not the")
println(io, "unpaid balance. `final` is the difference at the last such instant, chosen")
println(io, "before looking; `max` is the largest over all of them, and its `n` is the")
println(io, "multiplicity that number was selected out of.")
println(io)
println(io, "| drain (s) | final | ± SE | pool | max | ± SE | at t (s) | n | resolved? | vs D10's 1% |")
println(io, "|---|---|---|---|---|---|---|---|---|---|")
for r in drain_rows
    @printf(io, "| %.0f | %+.5f | %.5f | %s | %+.5f | %.5f | %d | %d | %s | %s |\n",
            r.drain, r.final, r.final_se, r.final_pool, r.max, r.max_se,
            r.max_t, r.n_points, resolved(r) ? "yes" : "no", verdict(r))
end
println(io)
println(io, "And the quantity that is **not** a granularity cost, reported separately so it")
println(io, "cannot be mistaken for one: the sawtooth of the outstanding debit, measured one")
println(io, "handshake before a drain, beside its closed form.")
println(io)
println(io, "| drain (s) | mid-period offset | pool | closed form |")
println(io, "|---|---|---|---|")
for r in drain_rows
    @printf(io, "| %.0f | %+.5f | %s | %+.5f |\n",
            r.drain, r.sawtooth, r.sawtooth_pool, r.closed_form)
end
println(io)
println(io, """
That offset is deterministic — `(drain − interval) · cost_rate / pool` — needs
no simulation, and is what an unaligned sampling instant measures. A first
version of this script reported it as the granularity cost. It is bookkeeping:
the cost has been accrued and not yet paid, not a trajectory that went
somewhere else.""")
println(io)
println(io, "The driver's own policy, as `reduction_report(models, driver)` renders it for")
println(io, "the 60 s configuration:")
println(io)
for line in split(labels_60, "\n")
    println(io, "    ", line)
end
println(io)
println(io, """
**Read narrowly.** This is the toy, not Core A′. Its ATP pool is ~404,000
particles against a drain of order 40 particles per second, a buffer of hours,
while the pools D10's comparison is really about turn over in 109 s (adenylate)
and 30 s (guanylate) — which is why D10 expects 60 s to fail there and 5 s to
pass. A small difference here bounds the *mechanism* — that the debit
aggregates, clamps and pays correctly — and nothing else. The measurement that
decides the granularity is phase 13's, on the assembled model; what this phase
delivers is the knob, the labelling and the method: drain-aligned sampling,
paired per-seed differences, the multiplicity of any maximum stated, and the
sawtooth reported apart from the dynamics.""")

s = String(take!(io))
print(s)
write(joinpath(@__DIR__, "rebuild_channel_result.md"), s)

gain_ok || error("The elasticity diagnostic disagrees with its closed form; see the table above")
