# KIM-90 item 4b: does sharing cells between the transcript and bulk streams explain
# the m0_v6 gap? (Decision log 2026-11-17: measured by cross-covariance, not a
# disjoint-cell refit.)
#
# The mixed route treats the two streams as independent given θ. In the exact model
# the bulk is the mean of the same 50 cells whose transcripts are observed, so the
# transcripts shift the bulk's expected value. To first order that shift is
#   Δm = log(1 + B (t̄_obs − t̄) / x̄),   B = Cov(x, t) Cov(t)⁻¹,
# with x a cell's panel (pools × saves), t its transcripts (genes × saves), the
# moments from the reference ensemble at θ_ref and t̄_obs the mean over m0_v6's 50
# cells. The mixed posterior ignores Δm, so it is displaced from the exact one by
# about Σ Jᵀ W Δm: Σ is V6's posterior covariance (transcripts' information
# carried), J = ∂m/∂log θ at m0_v6's posterior mean, W = (σ_b² I + S)⁻¹. This is
# jensen_bias_89's "Σ_F10" linear response, applied to Δm instead of the Jensen gap.
#
# Read, fixed before the analysis is run: for each of S_0779, krnadeg and ENO, the
# predicted shift's 90% bootstrap interval (cells resampled, 500 times) either
#   - excludes the observed gap and lies within ±0.15 SD of 0: shared cells ruled out;
#   - contains the observed gap and excludes 0: shared cells named as a cause;
#   - otherwise: indeterminate (the fallback, a disjoint-cell exact run, is Tom's call).
#
# Stages:
#   ensemble <out>  the archived 800 reference cells re-simulated at θ_ref on their
#                   seeds, keeping transcripts and panel; the panel is checked against
#                   the archive's cells. To <out>/crosscov/ensemble.jls.
#   analyse <out>   the regression, Δm, J, the linear response and the bootstrap,
#                   to <out>/crosscov/crosscov.md.
#
# Usage: julia --project dev/scripts/gap_diag_crosscov.jl <stage> <out>
# SMOKE=1: 40 cells, 20 bootstrap draws.

using InferCell
using InferCell: _mixed_solve
using Random
using Serialization
using Statistics
using LinearAlgebra
using Printf
using Dates
using MCMCChains

const SMOKE = get(ENV, "SMOKE", "0") == "1"
const COMMIT = strip(read(`git rev-parse --short HEAD`, String))
const B7_RUN = "/fred/oz022/tkimpson/infercell-89-b7-archive-2026-10-27/run/v6"
const F10_DRAWS = joinpath(@__DIR__, "f10_16a10_draws.tsv")
const NBOOT = SMOKE ? 20 : 500
const NCELL_ENS = SMOKE ? 40 : 800
const GAP = (S_0779 = -0.468, krnadeg = +0.564, ENO = +0.419)   # mean shifts, item 2
const PN = ["S_0607", "S_0779", "krnadeg", "ENO", "FBA", "σ_b"]
slurm_job() = get(ENV, "SLURM_ARRAY_JOB_ID", get(ENV, "SLURM_JOB_ID", "none"))
ccdir(out) = joinpath(out, "crosscov")

function run_ensemble(out)
    mkpath(ccdir(out))
    p = deserialize(joinpath(B7_RUN, "prep.jls"))
    ms = collect(AbstractSubModel, m0_models())
    seeds = p.seeds[1:NCELL_ENS]
    build() = (d = build_problem(ms; tspan = (0.0, Float64(M0_HORIZON_S)), complete = true);
               set_parameters!(d, ms, p.values); (ms, d))
    t = @elapsed e = emit_ensemble(build, seeds)
    keep = findall(>(0), e.times)
    panel = e.counts[:, state_rows(collect(M0_PANEL), e.species), keep]
    tx = e.transcripts[:, :, keep]
    diff = maximum(abs.(panel .- p.cells[1:NCELL_ENS, :, :]))
    serialize(joinpath(ccdir(out), "ensemble.jls"),
              (; seeds, genes = e.genes, panel, tx, max_panel_diff = diff, commit = COMMIT,
               job = slurm_job(), seconds = t))
    @printf("%s ensemble: %d cells in %.0f s; genes %s; max |panel − archived| = %.3g\n",
            now(), length(seeds), t, e.genes, diff)
end

function v6_cov()
    lines = readlines(F10_DRAWS); head = split(lines[1], '\t')
    nm = ["S_JCVISYN3A_0607", "S_JCVISYN3A_0779", "krnadeg", "kcatF_R_ENO", "kcatF_R_FBA", "sigma_b"]
    cols = [findfirst(==(n), head) for n in nm]
    X = reduce(hcat, [r[cols] for r in (parse.(Float64, split(l, '\t')) for l in lines[2:end])
                      if r[1] in 1:4 && r[2] >= 101])'
    return cov(X), vec(std(X; dims = 1))
end

# Δm (pools × saves, vec, pools fastest) from cells `idx` of the ensemble.
function delta_m(panel, tx, idx, tobs)
    n = length(idx)
    X = reshape(panel[idx, :, :], n, :)                  # cells × (pools·saves), pools fastest
    T = reshape(Float64.(tx[idx, :, :]), n, :)
    keepT = findall(>(0), vec(var(T; dims = 1)))         # drop a constant transcript column
    T = T[:, keepT]
    x̄, t̄ = vec(mean(X; dims = 1)), vec(mean(T; dims = 1))
    Cxt = cov(X, T); Ctt = cov(T)
    B = Cxt / Symmetric(Ctt)
    return log.(max.(x̄ .+ B * (tobs[keepT] .- t̄), 1.0)) .- log.(max.(x̄, 1.0)), keepT
end

function run_analyse(out)
    p = deserialize(joinpath(B7_RUN, "prep.jls"))
    e = deserialize(joinpath(ccdir(out), "ensemble.jls"))
    ds = p.ds
    gi = [findfirst(==(g), ds.genes) for g in e.genes]
    keep = findall(>(0), ds.times)                        # the panel's saves, without t = 0
    obs = Float64.(ds.transcripts[:, gi, :])
    size(obs, 3) == size(e.tx, 3) || (obs = obs[:, :, end-size(e.tx, 3)+1:end])
    tobs = vec(mean(reshape(obs, size(obs, 1), :); dims = 1))
    n = size(e.panel, 1)
    Δm, keepT = delta_m(e.panel, e.tx, 1:n, tobs)

    s = with_reference(mixed_setup(ds), p.delta, p.S)
    draws = cat((Array(deserialize(joinpath(B7_RUN, "chain_$c.jls")).chain[:, [Symbol("x[$i]") for i in 1:6], 1])
                 for c in 1:2)...; dims = 3)
    x̄ = vec(mean(draws; dims = (1, 3)))
    K = length(s.names)
    mf(x) = mean_field_log_bulk(s, _mixed_solve(s, exp.(x[1:K])).sol) |> v -> reduce(vcat, v)
    h = 1e-3
    J = zeros(length(Δm), 6)
    for k in 1:K
        xp, xm = copy(x̄), copy(x̄); xp[k] += h; xm[k] -= h
        J[:, k] = (mf(xp) .- mf(xm)) ./ 2h
    end
    σ = exp(x̄[end])
    W = inv(Symmetric(σ^2 * I + s.S))
    Σ, sd = v6_cov()
    resp(dm) = (Σ * (J' * (W * dm))) ./ sd
    pred = resp(Δm)

    rng = Xoshiro(17_090_090)
    boot = reduce(hcat, [resp(first(delta_m(e.panel, e.tx, rand(rng, 1:n, n), tobs))) for _ in 1:NBOOT])
    lo, hi = [quantile(boot[q, :], 0.05) for q in 1:6], [quantile(boot[q, :], 0.95) for q in 1:6]

    io = IOBuffer(); say(s_) = (println(s_); println(io, s_))
    say("# KIM-90 item 4b: transcript–bulk cross-covariance, sized against the gap\n")
    say("Ensemble job $(e.job) at $(e.commit): $n cells at θ_ref on the archived seeds; " *
        @sprintf("max |panel − archived| = %.3g. ", e.max_panel_diff) *
        "Analysis job $(slurm_job()) at $COMMIT. $(length(keepT)) transcript features " *
        "($(e.genes), saves); $NBOOT bootstrap resamples of the cells.\n")
    say(@sprintf("‖Δm‖ = %.4g (RMS %.3g per row); RMS of Δm in units of the bulk noise (W^½ Δm): %.3g.\n",
                 norm(Δm), sqrt(mean(abs2, Δm)), sqrt(dot(Δm, W * Δm) / length(Δm))))
    say("| Parameter (log) | predicted shift (V6 SDs) | 90% bootstrap | observed gap | reading |")
    say("|---|---|---|---|---|")
    obsgap = Dict("S_0779" => GAP.S_0779, "krnadeg" => GAP.krnadeg, "ENO" => GAP.ENO)
    for q in 1:6
        g = get(obsgap, PN[q], nothing)
        rd = g === nothing ? "—" :
             (!(lo[q] <= g <= hi[q]) && -0.15 <= lo[q] && hi[q] <= 0.15) ? "ruled out" :
             (lo[q] <= g <= hi[q] && !(lo[q] <= 0 <= hi[q])) ? "named" : "indeterminate"
        say(@sprintf("| %s | %+.3f | [%+.3f, %+.3f] | %s | %s |", PN[q], pred[q], lo[q], hi[q],
                     g === nothing ? "—" : @sprintf("%+.3f", g), rd))
    end
    write(joinpath(ccdir(out), "crosscov.md"), take!(io))
end

if abspath(PROGRAM_FILE) == (@__FILE__)
    stage, out = ARGS[1], ARGS[2]
    stage == "ensemble" ? run_ensemble(out) :
    stage == "analyse" ? run_analyse(out) : error("unknown stage $stage")
end
