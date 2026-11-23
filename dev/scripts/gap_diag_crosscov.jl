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
# cells. `analyse` carries Δm through V6's posterior covariance as Σ Jᵀ W Δm (Σ V6's
# posterior covariance, J = ∂m/∂log θ at m0_v6's posterior mean, W = (σ_b² I + S)⁻¹),
# jensen_bias_89's "Σ_F10" linear response.
#
# That response is incomplete (PR #189 review, 2026-11-23). Δm depends on θ through
# t̄(θ), G = ∂Δm/∂log θ ≈ −B ∂t̄/∂log θ / x̄, and the first-order displacement is
# Σ (Jᵀ W Δm − Gᵀ W r) with r the bulk residual, with W held. Given the transcripts the
# exact bulk covariance is also smaller (S → S_c), which adds a first-order
# Jᵀ (W_c − W) r. The `full` stage computes both, the bulk's information with and
# without conditioning on the transcripts, and the curvature that says whether a
# one-step response is valid.
# `analyse` is kept as run at 66431db; its reading rule applies to its own output.
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
#   full <out>      the full first-order term, the information shares and the
#                   Hessian check, to <out>/crosscov/full.md. No bootstrap.
#
# Usage: julia --project dev/scripts/gap_diag_crosscov.jl <stage> <out>
# SMOKE=1: 40 cells, 20 bootstrap draws.

using InferCell
using InferCell: _mixed_solve, _state_at
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
const GAP = (S_0779 = -0.468, krnadeg = +0.564, ENO = +0.419)   # item 2; recomputed and asserted below
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
    return cov(X), vec(std(X; dims = 1)), vec(mean(X; dims = 1))
end

function m0v6_draws()
    cat((Array(deserialize(joinpath(B7_RUN, "chain_$c.jls")).chain[:, [Symbol("x[$i]") for i in 1:6], 1])
         for c in 1:2)...; dims = 3)
end

# The gap recomputed from the draws, asserted against GAP to 3 dp.
function check_gap(x̄, v6mean, sd)
    g = (x̄ .- v6mean) ./ sd
    for (q, k) in ((2, :S_0779), (3, :krnadeg), (4, :ENO))
        abs(g[q] - getfield(GAP, k)) < 5e-4 ||
            error(@sprintf("gap on %s recomputes to %+.4f, not %+.3f", k, g[q], getfield(GAP, k)))
    end
    return g
end

# The observed transcripts on the ensemble's saves, matched by time, not position.
function observed_tx(ds, s, e)
    gi = [findfirst(==(g), ds.genes) for g in e.genes]
    ti = [findfirst(==(τ), ds.times) for τ in s.saves]
    any(isnothing, ti) && error("a mixed save is missing from the dataset's times")
    size(e.tx, 3) == length(s.saves) == size(e.panel, 3) ||
        error("ensemble saves $(size(e.tx, 3)) / panel $(size(e.panel, 3)) / mixed saves $(length(s.saves))")
    return Float64.(ds.transcripts[:, gi, ti])
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
    s = with_reference(mixed_setup(ds), p.delta, p.S)
    obs = observed_tx(ds, s, e)
    tobs = vec(mean(reshape(obs, size(obs, 1), :); dims = 1))
    n = size(e.panel, 1)
    Δm, keepT = delta_m(e.panel, e.tx, 1:n, tobs)

    draws = m0v6_draws()
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
    Σ, sd, v6mean = v6_cov()
    check_gap(x̄, v6mean, sd)
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

# The full first-order term, the bulk's information with and without conditioning
# on the transcripts, and the curvature check. Point values only.
function run_full(out)
    p = deserialize(joinpath(B7_RUN, "prep.jls"))
    e = deserialize(joinpath(ccdir(out), "ensemble.jls"))
    ds = p.ds
    s = with_reference(mixed_setup(ds), p.delta, p.S)
    obs = observed_tx(ds, s, e)
    tobs = vec(mean(reshape(obs, size(obs, 1), :); dims = 1))
    n = size(e.panel, 1)
    Δm, keepT = delta_m(e.panel, e.tx, 1:n, tobs)
    X = reshape(e.panel, n, :); T = reshape(Float64.(e.tx), n, :)
    length(keepT) == size(T, 2) || error("$(size(T, 2) - length(keepT)) constant transcript columns; B and Δm would differ")
    xe, te = vec(mean(X; dims = 1)), vec(mean(T; dims = 1))
    Cxt, Ctt = cov(X, T), cov(T)
    B = Cxt / Symmetric(Ctt)

    rt = s.route
    trows = [findfirst(==(InferCell.transcript_state(g)), rt.state_names) for g in e.genes]
    K = length(s.names)
    function mf(x)                                   # mean-field log bulk and transcripts, saves × genes
        sol = _mixed_solve(s, exp.(x[1:K])).sol
        m = reduce(vcat, mean_field_log_bulk(s, sol))
        t = reduce(vcat, [[_state_at(sol, τ)[i] for i in trows] for τ in s.saves])
        (m, t)
    end
    # G is −B ∂t̄/∂log θ / x̄ with ∂t̄ taken from the mean field, so the mean field's
    # transcripts must be on the ensemble's scale (a units check: the mean field and
    # the ensemble mean differ by about 13% at most, and G carries that error).
    _, tref = mf(p.ref.x)
    trel = maximum(abs.(tref .- te) ./ te)
    trel < 0.25 || error(@sprintf("mean-field transcripts at θ_ref differ from the ensemble mean by up to %.3g", trel))

    draws = m0v6_draws()
    x̄ = vec(mean(draws; dims = (1, 3)))
    Σ, sd, v6mean = v6_cov()
    gap = check_gap(x̄, v6mean, sd)
    h = 1e-3
    J = zeros(length(Δm), 6); Tj = zeros(length(te), 6)
    for k in 1:K
        xp, xm = copy(x̄), copy(x̄); xp[k] += h; xm[k] -= h
        (mp, tp), (mm, tm) = mf(xp), mf(xm)
        J[:, k] = (mp .- mm) ./ 2h; Tj[:, k] = (tp .- tm) ./ 2h
    end
    σ = exp(x̄[end])
    W = inv(Symmetric(σ^2 * I + s.S))
    C = size(ds.transcripts, 1)
    maximum(abs.(cov(X) ./ (xe * xe') ./ C .- s.S)) < 1e-10 * maximum(abs.(s.S)) ||
        error("S recomputed from the ensemble panel differs from the archived S")
    Sc = (cov(X) .- Cxt * (Ctt \ Cxt')) ./ (xe * xe') ./ C   # S conditioned on the transcripts
    Wc = inv(Symmetric(σ^2 * I + Symmetric(Sc)))
    G = -(B * Tj) ./ xe
    m̄, _ = mf(x̄)
    r = s.logz .- m̄ .- s.delta
    resp(v) = (Σ * v) ./ sd
    kept = resp(J' * (W * Δm))                       # Δm at θ_ref, as analyse
    kept_x = resp(J' * (W * (Δm .+ G * (x̄ .- p.ref.x))))   # Δm carried to x̄, where J and r are taken
    missing_ = resp(-(G' * (W * r)))
    # With the exact route's bulk covariance given the transcripts (σ²I + Sc) as well:
    # mixed − exact ≈ Σ (∇ log L_mixed − ∇ log L_exact) at x̄, both gradients in full.
    Δmx = Δm .+ G * (x̄ .- p.ref.x)
    full_c = resp(J' * (W * r) .- (J + G)' * (Wc * (r .- Δmx)))
    P = inv(Σ)
    Im = diag(J' * W * J)[1:K] ./ diag(P)[1:K]
    Ic = diag((J + G)' * W * (J + G))[1:K] ./ diag(P)[1:K]
    Icc = diag((J + G)' * Wc * (J + G))[1:K] ./ diag(P)[1:K]
    H = G' * W * J + J' * W * G + G' * W * G             # exact − mixed bulk precision, W held
    Hc = (J + G)' * Wc * (J + G) - J' * W * J             # the same with S conditioned
    Σh = sqrt(Symmetric(Σ))
    ev = eigvals(Symmetric(Σh * H * Σh))
    evc = eigvals(Symmetric(Σh * Hc * Σh))

    io = IOBuffer(); say(s_) = (println(s_); println(io, s_))
    say("# KIM-90 item 4b, the full first-order term and the bulk's information\n")
    say("Ensemble job $(e.job) at $(e.commit) ($n cells). This job $(slurm_job()) at $COMMIT. " *
        "Point values at m0_v6's posterior mean, no bootstrap.\n")
    say(@sprintf("Mean-field transcripts at θ_ref against the ensemble mean: max relative difference %.3g. ", trel) *
        "$(length(keepT)) of $(size(T, 2)) transcript columns kept.\n")
    say(@sprintf("RMS in bulk-noise units: W^½ Δm %.3g, W^½ r %.3g. V6 posterior SDs (log): %s.\n",
                 sqrt(dot(Δm, W * Δm) / length(Δm)), sqrt(dot(r, W * r) / length(r)),
                 join([@sprintf("%s %.4f", PN[q], sd[q]) for q in 1:6], ", ")))
    say("θ_ref − m0_v6's posterior mean, V6 SDs: " *
        join([@sprintf("%s %+.2f", PN[q], (p.ref.x[q] - x̄[q]) / sd[q]) for q in 1:6], ", ") *
        ". (δ and S come from simulated cells, which σ_b does not enter.)\n")
    say("## Displacement, V6 SDs (mixed − exact, first order)\n")
    say("Δm(x̄) = Δm + G (x̄ − x_ref) is Δm carried from θ_ref to m0_v6's posterior mean, " *
        "where J and r are evaluated; the sum is Σ (Jᵀ W Δm(x̄) − Gᵀ W r).\n")
    say("The 'sum' holds W at the mixed route's (σ²I + S)⁻¹. The last column also uses the " *
        "exact route's bulk covariance given the transcripts, σ²I + Sc: " *
        "Σ (JᵀW r − (J+G)ᵀW_c (r − Δm(x̄))).\n")
    say("| Parameter | gap (recomputed) | Σ JᵀWΔm, Δm at θ_ref (analyse) | Σ JᵀWΔm(x̄) | −Σ GᵀWr (omitted by analyse) | sum, W held | full, S conditioned |")
    say("|---|---|---|---|---|---|---|")
    for q in 1:6
        say(@sprintf("| %s | %+.3f | %+.3f | %+.3f | %+.3f | %+.3f | %+.3f |", PN[q], gap[q], kept[q], kept_x[q],
                     missing_[q], kept_x[q] + missing_[q], full_c[q]))
    end
    say(@sprintf("\ntr S_c / tr S = %.3f; tr S / tr σ²I = %.3f.\n", tr(Sc) / tr(s.S), tr(s.S) / (σ^2 * length(r))))
    say("\n## The bulk's information as a share of V6's posterior precision (diagonal)\n")
    say("Mixed: JᵀWJ, the bulk treated as independent of the transcripts. Given the transcripts, " *
        "first order: (J+G)ᵀW(J+G) with W held, and with S conditioned on the transcripts.\n")
    say("| Parameter | ‖J+G‖/‖J‖ | mixed | given transcripts, W held | given transcripts, S conditioned |")
    say("|---|---|---|---|---|")
    for k in 1:K
        say(@sprintf("| %s | %.3f | %.3f | %.3f | %.3f |", PN[k], norm(J[:, k] + G[:, k]) / norm(J[:, k]), Im[k], Ic[k], Icc[k]))
    end
    say("\n## Validity of a one-step response\n")
    say("Eigenvalues of Σ^½ H Σ^½, H the exact route's bulk precision minus the mixed route's " *
        "(W held: GᵀWJ + JᵀWG + GᵀWG; S conditioned: (J+G)ᵀW_c(J+G) − JᵀWJ), relative to V6's " *
        "posterior precision. A negative eigenvalue is a direction the exact posterior is less " *
        "precise in; one near ±1 means the change is as large as the posterior's own curvature, " *
        "so the one-step response does not hold.\n")
    say("W held: " * join([@sprintf("%+.3f", v) for v in ev], ", "))
    say("\nS conditioned: " * join([@sprintf("%+.3f", v) for v in evc], ", "))
    write(joinpath(ccdir(out), "full.md"), take!(io))
end

if abspath(PROGRAM_FILE) == (@__FILE__)
    stage, out = ARGS[1], ARGS[2]
    stage == "ensemble" ? run_ensemble(out) :
    stage == "analyse" ? run_analyse(out) :
    stage == "full" ? run_full(out) : error("unknown stage $stage")
end
