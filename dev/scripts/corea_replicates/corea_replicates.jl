# KIM-178: replicate truths on the reduced model, by MAP and Laplace (Linear M11). The
# design, the seeds, the failure rule and gate 2 are frozen in spec/prereg/kim-178.md at
# tag prereg/kim-178, before any of this ran. This script implements them and changes
# none of them.
#
# Usage: julia --project dev/scripts/corea_replicates/corea_replicates.jl <stage> <out> [<arg>]
# Submit it through dev/scripts/corea_replicates/corea_replicates.sh, from a worktree at a
# pushed commit.
#
# Stages:
#   replicate <out> <seed>
#       One replicate in <out>/seed_<seed>/, in four steps (prereg §2):
#         1. corea_dataset(seed; …)                  -> dataset.jls
#         2. θ_ref, then the 800-cell ensemble, δ, S -> prep.jls
#         3. the fitted-likelihood MAP from θ_ref    -> map.jls
#         4. −H at that MAP                          -> laplace.jls
#       It then judges the fit against the prereg §3 failure rule, writing fit.jls and
#       fit.md. A step whose .jls exists is skipped, so a killed task resumes. A failure
#       is recorded, never retried.
#   kim92 <out> [<src>]
#       Prereg §5: steps 3 and 4 on KIM-92's dataset, read from <src> (default: the KIM-92
#       archive), into <out>/kim92/. Then the MAP check, the Laplace SD ratios against
#       KIM-92's 1,600 pooled NUTS draws, and the promoters' mean z by each, in kim92.md.
#
# SMOKE=1 is a plumbing pass, never a result:
#   - replicate: 3 cells, a 300 s run cut at 150 s, 8 reference cells and 5 MAP iterations;
#   - kim92 (pass <src> = a smoke replicate's directory): 5 MAP iterations.
# A smoke runs on seed 17_200_199 only, which is outside the frozen set (prereg §3).

include(joinpath(@__DIR__, "..", "pilot_common.jl"))
using InferCell: _mixed_solve
using ForwardDiff
using Printf
using LinearAlgebra
using Statistics: mean, std
using DelimitedFiles

# prereg §2 and §3
const SEEDS = 17_200_002:17_200_011
const SPARES = 17_200_012:17_200_016
const SMOKE_SEED = 17_200_199
const KIM92_SEED = 17_200_001
const KIM92_ARCHIVE = "/fred/oz022/tkimpson/infercell-kim92-archive-2026-11-17"
const NCELLS = SMOKE ? 3 : 200
const HORIZON = SMOKE ? 300.0 : COREA_CYCLE_S
const FIT_TO = SMOKE ? 150.0 : COREA_FIT_S
const NREF = SMOKE ? 8 : 800
const MAP_ITERS = SMOKE ? 5 : 200
const FD_STEP = 1e-3                    # × the coordinate's prior SD, as KIM-92's metric
const BOX_MARGIN = 0.01                 # prior SDs from the support box (prereg §3)

pilot_setup(ds; kw...) = mixed_setup(ds; models = d11_models(), targets = COREA_TARGETS, kw...)

function truth_x(ds, s)
    tv = Dict(ds.truth.values)
    return vcat([log(tv[n]) for n in s.names], log(ds.sigma))
end

has(dir, name) = isfile(joinpath(dir, "$name.jls"))
promoters(names) = findall(startswith.(names, "log_S_"))

# KIM-92's `hessian_fd`, copied unchanged: central differences of the gradient, then
# symmetrised. The asymmetry before symmetrising shows the differences' noise.
function hessian_fd(grad, x, h)
    P = length(x)
    H = zeros(P, P)
    for j in 1:P
        e = zeros(P); e[j] = h[j]
        H[:, j] = (grad(x .+ e) .- grad(x .- e)) ./ (2h[j])
        say(@sprintf("laplace: column %d of %d", j, P))
    end
    asym = maximum(abs, H .- H') / maximum(abs, H)
    return (H .+ H') ./ 2, asym
end

# ---------------------------------------------------------------- the four steps

function step_dataset(dir, seed)
    has(dir, "dataset") && return
    mkpath(dir)
    t = @elapsed ds = corea_dataset(seed; ncells = NCELLS, horizon = HORIZON, fit_to = FIT_TO,
                                    targets = COREA_TARGETS)
    save(dir, "dataset", ds; stage = "dataset")
    say(@sprintf("dataset: seed %d, %d cells, %d saves to %.0f s, %d bulk rows, in %.0f s",
                 seed, NCELLS, length(ds.times), ds.times[end], length(ds.bulk), t))
end

# KIM-92's `run_prep`, on this replicate.
function step_prep(dir)
    has(dir, "prep") && return
    ds = load(dir, "dataset")
    s0 = pilot_setup(ds)
    t_map = @elapsed ref = reference_theta(s0; maxiters = MAP_ITERS)
    say(@sprintf("prep: θ_ref by %s in %.0f s, log density %.4f", ref.method, t_map, ref.logdensity))
    values = mixed_theta_values(s0, ref.x)
    m_ref = mean_field_log_bulk(s0, _mixed_solve(s0, exp.(ref.x[1:length(s0.names)])).sol)
    seeds = reference_seeds(ds, NREF)
    t_ens = @elapsed cells = reference_ensemble(values, seeds; models = d11_models(),
                                                horizon = ds.times[end])
    size(cells, 3) == length(ds.times) || error("the ensemble's saves are not the window's")
    d = delta_and_S(cells, m_ref; C = NCELLS)
    save(dir, "prep", (; ref, values, m_ref, seeds, cells, delta = d.delta, S = d.S,
                       seconds = (map = t_map, ensemble = t_ens)); stage = "prep")
    say(@sprintf("prep: %d reference cells to %.0f s in %.0f s; tr S %.4g; max |δ| %.4g",
                 length(seeds), ds.times[end], t_ens, tr(d.S), maximum(abs, d.delta)))
end

# Prereg §2 step 5: the MAP of the fitted likelihood (δ and S in place), from θ_ref.
function step_map(dir, ds, p)
    has(dir, "map") && return
    s = with_reference(pilot_setup(ds), p.delta, p.S)
    lp_ref = mixed_logdensity(s, p.ref.x)
    t = @elapsed fit = reference_theta(s; x0 = p.ref.x, maxiters = MAP_ITERS)
    save(dir, "map", (; fit, lp_ref, seconds = t); stage = "map")
    say(@sprintf("map: by %s in %.0f s; fitted log density %.4f at the MAP, %.4f at θ_ref",
                 fit.method, t, fit.logdensity, lp_ref))
end

# Prereg §2 step 6: −H at the fitted MAP, no floor.
function step_laplace(dir, ds, p)
    has(dir, "laplace") && return
    s = with_reference(pilot_setup(ds), p.delta, p.S)
    x = load(dir, "map").fit.x
    sd = mixed_prior_sd(s)
    grad(z) = ForwardDiff.gradient(y -> mixed_logdensity(s, y), z)
    t = @elapsed (H, asym) = hessian_fd(grad, x, FD_STEP .* sd)
    ev = eigvals(Symmetric(-H))
    pd = ev[1] > 0
    lap_sd = pd ? sqrt.(diag(inv(Symmetric(-H)))) : fill(NaN, length(x))
    save(dir, "laplace", (; H, asym, eigenvalues = ev, positive_definite = pd, lap_sd,
                          fd_step = FD_STEP, seconds = t); stage = "laplace")
    say(@sprintf("laplace: %.0f s; eigenvalues of −H %.4g to %.4g; asymmetry %.2g",
                 t, ev[1], ev[end], asym))
end

# ---------------------------------------------------------------- the judgement

# Prereg §3's failure rule, and the z values of §1 (prereg §1, §6).
function judge(ds, p, m, l)
    s = pilot_setup(ds)
    names = mixed_parameter_names(s)
    sd = mixed_prior_sd(s)
    lb, ub = mixed_box(s)
    x = m.fit.x
    reasons = String[]
    p.ref.method === :map || push!(reasons, "θ_ref fell back to the prior median")
    m.fit.method === :map || push!(reasons, "the fitted-likelihood MAP fell back to the prior median")
    edge = min.(x .- lb, ub .- x) ./ sd
    near = findall(<(BOX_MARGIN), edge)
    isempty(near) || push!(reasons, "within $BOX_MARGIN prior SD of the box: " * join(names[near], ", "))
    l.positive_definite || push!(reasons, @sprintf("−H is not positive definite (smallest eigenvalue %.4g)",
                                                    l.eigenvalues[1]))
    all(isfinite, l.lap_sd) || push!(reasons, "a Laplace SD is not finite")
    xt = truth_x(ds, s)
    z = (x .- xt) ./ l.lap_sd
    z_ref = (p.ref.x .- xt) ./ l.lap_sd          # θ_ref against the same SDs (prereg §6)
    prom = promoters(names)
    # KIM-161's correlation-aware SE of the promoters' mean z, from (−H)⁻¹ (prereg §6).
    se_corr = if l.positive_definite
        C = inv(Symmetric(-l.H))
        Z = Diagonal(1 ./ l.lap_sd[prom]) * C[prom, prom] * Diagonal(1 ./ l.lap_sd[prom])
        sqrt(sum(Z)) / length(prom)
    else
        NaN
    end
    return (; ok = isempty(reasons), reasons, names, truth = xt, x, z, z_ref, prom,
            m = mean(z[prom]), m_ref = mean(z_ref[prom]), nneg = count(<(0), z[prom]),
            se_corr, edge_min = minimum(edge))
end

function report(dir, seed, ds, p, m, l, j)
    io = IOBuffer()
    pr(xs...) = println(io, xs...)
    pr("# KIM-178 replicate, seed ", seed, "\n")
    pr(@sprintf("Commit `%s`, job %s%s. Prereg `prereg/kim-178`.\n", COMMIT, JOB, SMOKE ? ", **SMOKE (not a result)**" : ""))
    pr("**Fit: ", j.ok ? "succeeded" : "FAILED", ".**", j.ok ? "" : " " * join(j.reasons, "; ") * ".", "\n")
    pr(@sprintf("- θ_ref by %s, %.0f s; the ensemble %.0f s", p.ref.method, p.seconds.map, p.seconds.ensemble))
    pr(@sprintf("- the fitted MAP by %s, %.0f s; fitted log density %.4f at the MAP, %.4f at θ_ref",
                m.fit.method, m.seconds, m.fit.logdensity, m.lp_ref))
    pr(@sprintf("- −H: %.0f s, eigenvalues %.4g to %.4g, asymmetry %.2g; nearest the box %.2f prior SD",
                l.seconds, l.eigenvalues[1], l.eigenvalues[end], l.asym, j.edge_min))
    pr(@sprintf("- promoters: mean z %.3f (θ_ref %.3f), %d of %d below the truth, correlation-aware SE %.3f\n",
                j.m, j.m_ref, j.nneg, length(j.prom), j.se_corr))
    pr("| coordinate | truth | MAP | Laplace SD | z | z at θ_ref |")
    pr("|---|---|---|---|---|---|")
    for i in eachindex(j.names)
        pr(@sprintf("| %s | %.4f | %.4f | %.4f | %+.2f | %+.2f |", j.names[i], j.truth[i], j.x[i],
                    l.lap_sd[i], j.z[i], j.z_ref[i]))
    end
    write(joinpath(dir, "fit.md"), String(take!(io)))
    print(read(joinpath(dir, "fit.md"), String))
end

function run_replicate(out, seed)
    seed in SEEDS || seed in SPARES || (SMOKE && seed == SMOKE_SEED) ||
        error("seed $seed is not one of prereg §3's (a smoke uses $SMOKE_SEED only)")
    SMOKE && seed != SMOKE_SEED && error("a smoke runs on seed $SMOKE_SEED only")
    !SMOKE && seed == SMOKE_SEED && error("seed $SMOKE_SEED is for smokes only")
    dir = joinpath(out, "seed_$seed")
    step_dataset(dir, seed)
    step_prep(dir)
    ds, p = load(dir, "dataset"), load(dir, "prep")
    step_map(dir, ds, p)
    step_laplace(dir, ds, p)
    m, l = load(dir, "map"), load(dir, "laplace")
    j = judge(ds, p, m, l)
    save(dir, "fit", (; seed, j...); stage = "fit")
    report(dir, seed, ds, p, m, l, j)
end

# ---------------------------------------------------------------- KIM-92's dataset

function run_kim92(out, src)
    dir = joinpath(out, "kim92")
    mkpath(dir)
    ds = deserialize(joinpath(src, "dataset.jls")).result
    p = deserialize(joinpath(src, "prep.jls")).result
    SMOKE || ds.replicate_seed == KIM92_SEED || error("$src is not KIM-92's dataset")
    step_map(dir, ds, p)
    step_laplace(dir, ds, p)
    m, l = load(dir, "map"), load(dir, "laplace")
    j = judge(ds, p, m, l)
    # KIM-92's 1,600 pooled NUTS draws (prereg §5), in the setup's coordinate order.
    D, hdr = readdlm(joinpath(KIM92_ARCHIVE, "pilot_draws.tsv"), '\t', Float64; header = true)
    vec(hdr) == j.names || error("pilot_draws.tsv's columns are not the setup's coordinates")
    size(D, 1) == 1600 || error("expected 1,600 pooled draws, found $(size(D, 1))")
    pm, ps = vec(mean(D; dims = 1)), vec(std(D; dims = 1))
    z_nuts = (pm .- j.truth) ./ ps
    ratio = l.lap_sd ./ ps
    # The Laplace SD at θ_ref, from KIM-92's warm metric (job 20651329), for the comparison.
    # A smoke's source is a smoke replicate, which has no metric: its own −H stands in.
    mf = joinpath(src, "metric.jls")
    Href = isfile(mf) ? deserialize(mf).result.H :
           SMOKE ? deserialize(joinpath(src, "laplace.jls")).result.H :
           error("$mf is missing")
    lap_ref = sqrt.(diag(inv(Symmetric(-Href))))
    save(dir, "kim92", (; j..., pm, ps, z_nuts, ratio, lap_ref, ratio_ref = lap_ref ./ ps,
                        lp_map = m.fit.logdensity, lp_ref = m.lp_ref); stage = "kim92")
    report(dir, KIM92_SEED, ds, p, m, l, j)
    io = IOBuffer()
    pr(xs...) = println(io, xs...)
    pr("# KIM-178 prereg §5: MAP and Laplace on KIM-92's dataset\n")
    pr(@sprintf("Commit `%s`, job %s%s. Source `%s`.\n", COMMIT, JOB, SMOKE ? ", **SMOKE (not a result)**" : "", src))
    pr(@sprintf("- **MAP check:** fitted log density %.4f at the MAP against %.4f at θ_ref (Δ %+.4f): %s.",
                m.fit.logdensity, m.lp_ref, m.fit.logdensity - m.lp_ref,
                m.fit.logdensity >= m.lp_ref ? "holds" : "FAILS"))
    pr(@sprintf("- **Promoters' mean z:** %.3f by MAP and Laplace, %.3f by NUTS (pooled draws), %.3f at θ_ref by its Laplace SD.",
                j.m, mean(z_nuts[j.prom]), j.m_ref))
    pr(@sprintf("- **Laplace SD / NUTS SD:** %.3f to %.3f at the fitted MAP; %.3f to %.3f at θ_ref.",
                extrema(ratio)..., extrema(lap_ref ./ ps)...))
    pr(@sprintf("- The fit %s under prereg §3.%s\n", j.ok ? "succeeds" : "FAILS",
                j.ok ? "" : " " * join(j.reasons, "; ")))
    pr("| coordinate | NUTS mean | NUTS SD | MAP | Laplace SD | ratio | ratio at θ_ref | z NUTS | z MAP |")
    pr("|---|---|---|---|---|---|---|---|---|")
    for i in eachindex(j.names)
        pr(@sprintf("| %s | %.4f | %.4f | %.4f | %.4f | %.3f | %.3f | %+.2f | %+.2f |", j.names[i],
                    pm[i], ps[i], j.x[i], l.lap_sd[i], ratio[i], lap_ref[i] / ps[i], z_nuts[i], j.z[i]))
    end
    write(joinpath(dir, "kim92.md"), String(take!(io)))
    print(read(joinpath(dir, "kim92.md"), String))
end

# ---------------------------------------------------------------- main

function main(args)
    length(args) >= 2 || error("usage: corea_replicates.jl <replicate|kim92> <out> [<arg>]")
    stage, out = args[1], abspath(args[2])
    if stage == "replicate"
        run_replicate(out, parse(Int, args[3]))
    elseif stage == "kim92"
        run_kim92(out, length(args) >= 3 ? abspath(args[3]) : KIM92_ARCHIVE)
    else
        error("unknown stage $stage")
    end
end

main(ARGS)
