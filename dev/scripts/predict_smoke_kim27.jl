# KIM-27's done-when: a smoke run of every stage of `predict`, from a worktree at
# a pushed commit. It is a check that each path runs end to end on the real
# archives and compositions, at smoke size. Nothing here is a result.
#
# Usage (one stage per job; dev/scripts/predict_smoke_kim27.sh submits them all):
#   julia --project dev/scripts/predict_smoke_kim27.jl <stage> <outdir>
# Stages:
#   quickstart    the getting-started "Predict" block, as written
#   ingest        the archived V7 replicate 1 and B7 r1 chains (KIM-28)
#   m0            M0's crossed cell on replicate 1: exact/mixed × hybrid/mean field,
#                 δ none and ensemble (KIM-31, KIM-27)
#   intervention  a ptsG knockdown at dose 0 (ΔQ ≡ 0) and 0.25, both forwards (KIM-32)
#   reduced       d11 completion at widths ½, 1, 2, and under a toy V (KIM-29, KIM-30)
#   wrapper       posterior_predictive's original call sites, test_txl.jl and
#                 test_joint_inference.jl (KIM-34)

using InferCell
using InferCell: AbstractSubModel
using Random
using Serialization
using Statistics: mean
using Printf
using Test

const STAGE = Symbol(ARGS[1])
const OUT = ARGS[2]
const COMMIT = strip(read(`git rev-parse --short HEAD`, String))
slurm_job() = get(ENV, "SLURM_ARRAY_JOB_ID", get(ENV, "SLURM_JOB_ID", "none"))

const V7 = "/fred/oz022/tkimpson/infercell-16a9b-archive-2026-10-17/v7/sbc_v7_16a9b/replicate1.jls"
const B7 = ["/fred/oz022/tkimpson/infercell-89-b7-archive-2026-10-27/run/r1/chain_$c.jls" for c in 1:2]
const H = 120.0

check(ok, what) = ok ? println("  ok: ", what) : error("smoke check failed: $what")

function report(name, p::Prediction)
    println("  $name: $p")
    for (w, s) in pairs(p.summary)
        @printf("    %-12s mean %.6g  var %.4g  mc %.4g  param %.4g%s  intrinsic %.4g  ess %.2f  pairs %d (no Q %d)\n",
                w, s.mean, s.var, s.mc, s.parameter_var, s.flagged ? " (flagged)" : "",
                s.intrinsic_var, s.ess, s.n_pairs, s.n_no_q)
    end
    if p.effect_summary !== nothing
        for (w, s) in pairs(p.effect_summary)
            @printf("    effect %-12s mean %.6g  param %.4g  mc %.4g\n", w, s.mean,
                    s.parameter_var, s.mc)
        end
    end
    println("    statuses ", p.statuses, "; seconds so far ", round(time() - T0; digits = 1))
    serialize(joinpath(OUT, "$(STAGE)_$name.jls"),
              (; prediction = p, commit = COMMIT, job = slurm_job(), stage = STAGE))
end

const T0 = time()

function quickstart()
    # docs/src/getting-started.md, the dataset of "Fit a synthetic dataset", then
    # the "Predict" block as written.
    ds = m0_dataset(16085; ncells = 3, horizon = 180)
    ms = m0_models()
    ch = infer(ms, ds; route = :exact, rng = Xoshiro(3), n_sweeps = 5, N = 4)
    pd = posterior_draws(ch)
    p = predict(ms, pd; horizon = 180, n_cells = 4, n_parameter_draws = 3, rng = Xoshiro(4))
    report("hybrid", p)
    m = predict(ms, pd; horizon = 180, forward = :mean_field, delta = :ensemble, n_ref = 8,
                rng = Xoshiro(4))
    report("mean_field_delta", m)
    println("  m.summary.uniform.mean, m.delta.delta = ", (m.summary.uniform.mean, m.delta.delta))
    k = predict(ms, pd; horizon = 180, n_cells = 4, n_parameter_draws = 3,
                intervention = Knockdown(; dose = 0.25), rng = Xoshiro(4))
    report("knockdown", k)
    check(isfinite(k.effect_summary.uniform.mean), "the quick start's effect is finite")
end

function ingest()
    e = posterior_draws(V7)
    m = posterior_draws(B7)
    println("  ", e, "\n  ", m)
    check(e.names == m.names == M0_TARGETS, "both posteriors name M0_TARGETS")
    check(length(e) == 556 && length(m) == 400, "556 exact draws, 400 mixed (2 × 200)")
    check(all(isapprox(e.truth[k], m.truth[k]; rtol = 1e-12) for k in keys(e.truth)),
          "the V7 and B7 truths agree")
    check(e.sigma_b !== nothing && m.sigma_b !== nothing, "σ_b is carried on both")
end

function m0()
    ms = m0_models()
    for (label, pd) in (("exact", posterior_draws(V7; thin = 10)), ("mixed", posterior_draws(B7)))
        hy = predict(ms, pd; horizon = H, n_cells = 4, n_parameter_draws = 3, rng = Xoshiro(27))
        report("$(label)_hybrid", hy)
        mf = predict(ms, pd; horizon = H, forward = :mean_field, n_parameter_draws = 3,
                     rng = Xoshiro(27))
        report("$(label)_mean_field", mf)
        mc = predict(ms, pd; horizon = H, forward = :mean_field, delta = :ensemble, n_ref = 8,
                     n_parameter_draws = 3, rng = Xoshiro(27))
        report("$(label)_mean_field_delta", mc)
        check(mc.mu == mf.mu, "$label: the corrected run keeps the uncorrected m(θ)")
        check(hy.draws == mf.draws, "$label: the same rng picks the same draws")
    end
end

function intervention()
    ms = m0_models()
    pd = posterior_draws(V7; thin = 10)
    for fwd in (:hybrid, :mean_field)
        z = predict(ms, pd; horizon = H, forward = fwd, n_cells = 3, n_parameter_draws = 2,
                    intervention = Knockdown(; dose = 0.0), rng = Xoshiro(32))
        check(all(==(0.0), z.effect), "$fwd: a zero dose gives ΔQ ≡ 0")
        k = predict(ms, pd; horizon = H, forward = fwd, n_cells = 3, n_parameter_draws = 2,
                    intervention = Knockdown(; dose = 0.25), rng = Xoshiro(32),
                    delta = fwd === :mean_field ? :ensemble : :none, n_ref = 4)
        report("$(fwd)_knockdown", k)
        check(all(isfinite, k.effect), "$fwd: a 25% knockdown's effect is finite")
    end
end

function reduced()
    dms = d11_models()
    pd = posterior_draws(B7)
    for w in (0.5, 1.0, 2.0)
        p = predict(dms, pd; horizon = H, completion = :prior, width = w, n_held = 2,
                    n_parameter_draws = 2, rng = Xoshiro(29))
        report("prior_width_$w", p)
        check(length(p.held) == 15, "width $w: 15 parameters are held (which ones: test_predict_completion.jl)")
    end
    V = ViabilityRule(; label = "smoke toy V, not M4's", allow_clipped = true,
                      min_viable_fraction = 0.5, floors = [:M_gtp_c => 500.0])
    v = predict(dms, pd; horizon = H, completion = :viability_conditioned_prior,
                viability = V, n_viability_cells = 0,   # M3's behaviour; required since M4
                n_held = 2, n_parameter_draws = 2, n_cells = 2,
                intervention = Knockdown(; dose = 0.25), rng = Xoshiro(30))
    report("viability", v)
    println("    conditioning ", v.conditioning)
    println("    decisions ", [d.counts for d in v.decisions])
    mf = predict(dms, pd; horizon = H, forward = :mean_field, completion = :prior, n_held = 2,
                 n_parameter_draws = 2, rng = Xoshiro(29))
    report("prior_mean_field", mf)
end

function wrapper()
    root = joinpath(@__DIR__, "..", "..", "test")
    @testset "posterior_predictive's original call sites" begin
        include(joinpath(root, "test_txl.jl"))
        include(joinpath(root, "test_joint_inference.jl"))
    end
end

function main()
    mkpath(OUT)
    println("KIM-27 smoke, stage $STAGE, at $COMMIT, job $(slurm_job()), ",
            "$(Threads.nthreads()) thread(s)")
    STAGE === :quickstart ? quickstart() :
    STAGE === :ingest ? ingest() :
    STAGE === :m0 ? m0() :
    STAGE === :intervention ? intervention() :
    STAGE === :reduced ? reduced() :
    STAGE === :wrapper ? wrapper() : error("unknown stage $STAGE")
    @printf("stage %s done in %.1f s\n", STAGE, time() - T0)
end

main()
