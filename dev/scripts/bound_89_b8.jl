# Issue #89 B8.5: AV7's information bound per truth, matched to the mixed fit's
# covariance (spec/phases/approximate-method.md §12 2026-10-29, "B8 as planned").
#
# A2's bound (`sensitivity_88_merge.jl`, on 15.8's sampler model at the nominal θ)
# weights each metabolite row by max(SE over C cells, σ_b · mean), diagonal, with
# the saves independent. The mixed likelihood weights the bulk record by
# σ_b² I + S, with S the delta method's full covariance of log x̄_C. This script
# keeps A2's Jacobian, central differences of ensemble means at ln θ ± 0.2 under
# common random numbers, and replaces the metabolite rows' weighting:
#
#   J_m = (log max(x̄₊, 1) − log max(x̄₋, 1)) / 2δ    per panel pool and save,
#   Σ   = σ_b² I + Cov(nominal cells) / (x̄ x̄ᵀ C),      x̄ the nominal mean,
#   I_m = J_mᵀ Σ⁻¹ J_m.
#
# σ_b is the injection's truth, and the bound is for D11's six at that known σ_b.
# Transcript rows keep A2's rule: ensemble means at a 1% floor, which understates
# what per-cell counts carry. The split-half bracket is A2's: the same calculation
# with the noise matrix (J_A − J_B)/2 from disjoint halves of the seeds in place of
# J. A target is informed if its ratio is below its bracket.
#
# The cells are the reference model AV7's data come from, `d11_models()` built as
# `generate_dataset` builds it, at injection r's truth (corea_draws(17_200_000 + r)),
# 1470 paired seeds as A2 (20001–21470).
#
# Stages:
#   cells <out> <r> <task> <ntask>   one array task of the sensitivity runs, to
#                                    <out>/r<r>/task_<task>.jls;
#   bound <out> <r>                  the bound from r's task files, to
#                                    <out>/r<r>/bound.jls and bound.md;
#   a2check <A2 archive>             the calculation under A2's rule on A2's own
#                                    task files, held to A2's merged.jls (B8.5's
#                                    verification clause).
#
# Usage: julia --project dev/scripts/bound_89_b8.jl <stage> ...
# SMOKE=1: a 120 s horizon and NSEEDS (default 4) seeds.

using InferCell
using Random
using Serialization
using LinearAlgebra
using Printf
using Dates
using Statistics: mean, std, cov
using Distributions: LogNormal, params

const SMOKE = get(ENV, "SMOKE", "0") == "1"
const COMMIT = strip(read(`git rev-parse --short HEAD`, String))
const DELTA = 0.2
const NSEEDS = parse(Int, get(ENV, "NSEEDS", SMOKE ? "4" : "1470"))
const SEEDS = 20_000 .+ (1:NSEEDS)
const HORIZON = SMOKE ? 120.0 : COREA_CYCLE_S
const NCELLS = 200
const DATA_SEED0 = 17_200_000
const DESIGNS = ("metabolites only", "transcripts only", "joint")
job() = get(ENV, "SLURM_ARRAY_JOB_ID", get(ENV, "SLURM_JOB_ID", "none"))
istx(o) = startswith(string(o), "mRNA_")

# ---------------------------------------------------------------- cells

# The writes for column c at sign s about the truth: c scaled by e^{sδ}, and a
# forward constant's reverse re-derived through Haldane, as m0_dataset's
# `forwards` derives it.
function config_values(ms, truth, c, s)
    s == 0 && return truth
    tv = Dict(truth)
    v = tv[c] * exp(s * DELTA)
    moved = c in (r.forward for r in haldane_relations(ms)) ?
            Dict(derived_ode_values(ms, [c => v])) : Dict(c => v)
    return [n => get(moved, n, x) for (n, x) in truth]
end

function run_cells(out, r, task, ntask)
    dir = joinpath(out, "r$r")
    mkpath(dir)
    ms = d11_models()
    truth = corea_draws(DATA_SEED0 + r; models = ms).truth
    configs = vcat([(:nominal, 0)], [(c, s) for c in D11_TARGETS for s in (+1, -1)])
    runs = Dict{Tuple{Int, Symbol, Int}, Matrix{Float64}}()
    species, times = Symbol[], Float64[]
    for (j, seed) in enumerate(SEEDS)
        (j - 1) % ntask == task || continue
        for (c, s) in configs
            Random.seed!(seed)
            d = build_problem(ms; tspan = (0.0, HORIZON), complete = true)
            set_parameters!(d, ms, config_values(ms, truth.values, c, s))
            Random.seed!(seed)
            e = emit_observables!(d, ms; every = 60)
            proteins = Set(protein_count_states(ms))
            keep = [i for (i, o) in enumerate(e.species) if !(o in proteins)]
            species, times = e.species[keep], e.times
            runs[(seed, c, s)] = e.counts[keep, :]
        end
        println("$(now()) r$r seed $seed done"); flush(stdout)
    end
    serialize(joinpath(dir, "task_$task.jls"),
              (; r, truth, commit = COMMIT, job = job(), smoke = SMOKE, delta = DELTA,
               horizon = HORIZON, species, times, columns = D11_TARGETS, runs))
end

# ---------------------------------------------------------------- the bound

function load_tasks(dir)
    files = sort(filter(f -> occursin(r"^task_\d+\.jls$", f), readdir(dir)))
    isempty(files) && error("no task files in $dir")
    parts = [deserialize(joinpath(dir, f)) for f in files]
    p1 = first(parts)
    all(p -> p.species == p1.species && p.times == p1.times && p.columns == p1.columns &&
             p.commit == p1.commit && p.job == p1.job, parts) || error("tasks disagree")
    return p1, merge((p.runs for p in parts)...)
end

function prior_sds(ms, cols)
    map(cols) do c
        pr = only(p for m in ms for p in parameters(m) if p.name === c).prior
        pr isa LogNormal || error("$c's prior is $(typeof(pr)), not lognormal")
        params(pr)[2]
    end
end

# Posterior over prior SD from a Fisher information F in ln θ.
ratio(F, prior_sd) = sqrt.(diag(inv(F + Diagonal(1 ./ prior_sd .^ 2)))) ./ prior_sd

# The bound's pieces from one run's layout: replicates × rows matrices per config.
function layout(p1, runs)
    seeds = sort(unique(k[1] for k in keys(runs)))
    nsp, nt = length(p1.species), length(p1.times)
    sp(r) = p1.species[(r - 1) % nsp + 1]
    ti(r) = (r - 1) ÷ nsp + 1
    live = [r for r in 1:(nsp * nt) if p1.times[ti(r)] > 0]
    met = Set(M0_PANEL)
    # Pools fastest within a save, in M0_PANEL's order, as the bulk record's vec.
    metrows = [r for t in sort(unique(ti.(live))) for o in M0_PANEL
               for r in live if ti(r) == t && sp(r) === o]
    length(metrows) == length(M0_PANEL) * length(unique(ti.(live))) ||
        error("the runs do not hold every panel pool at every save")
    txrows = [r for r in live if istx(sp(r))]
    mat(ss, c, s) = permutedims(reduce(hcat, [vec(runs[(k, c, s)]) for k in ss]))
    h = length(seeds) ÷ 2
    return (; seeds, A = seeds[1:h], B = seeds[(h + 1):(2h)], sp, metrows, txrows, mat)
end

# The metabolite rows' log-mean Jacobian on seeds ss, floored at one particle as the
# likelihood floors x̄.
function log_jacobian(L, cols, ss)
    lm(c, s) = log.(max.(vec(mean(L.mat(ss, c, s)[:, L.metrows]; dims = 1)), 1.0))
    return reduce(hcat, [(lm(c, +1) .- lm(c, -1)) ./ (2DELTA) for c in cols])
end

# The matched bound: metabolite rows under σ² I + S, transcript rows under A2's rule.
function matched_bound(p1, runs, cols, prior_sd; σ, ncells = NCELLS)
    L = layout(p1, runs)
    X = L.mat(L.seeds, :nominal, 0)[:, L.metrows]
    x̄ = vec(mean(X; dims = 1))
    S = cov(X) ./ (x̄ * x̄') ./ ncells
    W = cholesky(Symmetric(S + σ^2 * I))
    whiten(J) = W.L \ J
    Jm, JmA, JmB = (log_jacobian(L, cols, ss) for ss in (L.seeds, L.A, L.B))
    nominal = L.mat(L.seeds, :nominal, 0)
    tx(ss) = ensemble_jacobian([L.mat(ss, c, +1) for c in cols], [L.mat(ss, c, -1) for c in cols],
                               nominal; delta = p1.delta, ncells, rel = 0.01,
                               floors = ones(size(nominal, 2)))[L.txrows, :]
    Jt, JtA, JtB = tx(L.seeds), tx(L.A), tx(L.B)
    Jw, Nw = whiten(Jm), whiten((JmA .- JmB) ./ 2)
    Nt = (JtA .- JtB) ./ 2
    out = Dict{String, Any}()
    for (name, J, N) in (("metabolites only", Jw, Nw), ("transcripts only", Jt, Nt),
                         ("joint", vcat(Jw, Jt), vcat(Nw, Nt)))
        out[name] = (ratio = ratio(J' * J, prior_sd), bracket = ratio(N' * N, prior_sd),
                     rows = size(J, 1))
    end
    return out
end

# A2's rule (sensitivity_88_merge.jl), on the same layout: metabolite rows at
# resolution max(SE, cv · max(|mean|, 500)), transcript rows at 1% with a floor of 1.
function a2_bound(p1, runs, cols, prior_sd; cv, ncells = NCELLS)
    L = layout(p1, runs)
    nominal = L.mat(L.seeds, :nominal, 0)
    floors = [o === :volume_litres ? 0.0 : istx(o) ? 1.0 : 500.0
              for o in (L.sp(r) for r in 1:size(nominal, 2))]
    jac(ss, rel) = ensemble_jacobian([L.mat(ss, c, +1) for c in cols], [L.mat(ss, c, -1) for c in cols],
                                     nominal; delta = p1.delta, ncells, rel, floors)
    J1, NA1, NB1 = jac(L.seeds, 0.01), jac(L.A, 0.01), jac(L.B, 0.01)
    Jc, NAc, NBc = cv == 0.01 ? (J1, NA1, NB1) : (jac(L.seeds, cv), jac(L.A, cv), jac(L.B, cv))
    # A2 orders metabolite rows as they lie in the record; the bound does not depend on it.
    mr = sort(L.metrows)
    J = vcat(Jc[mr, :], J1[L.txrows, :])
    N = vcat((NAc[mr, :] .- NBc[mr, :]) ./ 2, (NA1[L.txrows, :] .- NB1[L.txrows, :]) ./ 2)
    nm = length(mr)
    out = Dict{String, Any}()
    for (name, rows) in (("metabolites only", 1:nm), ("transcripts only", (nm + 1):size(J, 1)),
                         ("joint", 1:size(J, 1)))
        out[name] = (ratio = ratio(J[rows, :]' * J[rows, :], prior_sd),
                     bracket = ratio(N[rows, :]' * N[rows, :], prior_sd), rows = length(rows))
    end
    return out
end

function table(io, b, cols)
    println(io, "| Design | Rows | ", join(String.(cols), " | "), " |")
    println(io, "|---|---|", repeat("---|", length(cols)))
    for name in DESIGNS
        x = b[name]
        println(io, "| $name | $(x.rows) | ",
                join([@sprintf("%.3f (%.3f)", a, c) for (a, c) in zip(x.ratio, x.bracket)], " | "), " |")
    end
end

function run_bound(out, r)
    dir = joinpath(out, "r$r")
    p1, runs = load_tasks(dir)
    ms = d11_models()
    d = corea_draws(DATA_SEED0 + r; models = ms)
    p1.truth.values == d.truth.values || error("r$r's task files are at another truth")
    cols = p1.columns
    psd = prior_sds(ms, cols)
    b = matched_bound(p1, runs, cols, psd; σ = d.σ)
    nseeds = length(unique(k[1] for k in keys(runs)))
    serialize(joinpath(dir, "bound.jls"),
              (; r, cols, sigma = d.σ, bound = b, nseeds, cells_commit = p1.commit,
               cells_job = p1.job, commit = COMMIT, job = job(), smoke = p1.smoke))
    io = IOBuffer()
    println(io, "# AV7's bound at r$r (σ_b = $(round(d.σ; sigdigits = 4))), $nseeds paired seeds\n")
    println(io, "Cells at `$(p1.commit)` (job $(p1.job)); bound at `$COMMIT` (job $(job())). ",
            "Posterior over prior SD, linearised, with the split-half bracket.\n")
    table(io, b, cols)
    write(joinpath(dir, "bound.md"), String(take!(io)))
    print(read(joinpath(dir, "bound.md"), String))
end

# B8.5's clause: A2's rule on A2's task files reproduces A2's merged.jls.
function run_a2check(a2dir)
    p1, runs = load_tasks(a2dir)
    p1.model === :corea || error("$a2dir is not A2's Core A′ run")
    m = deserialize(joinpath(a2dir, "merged.jls"))
    cols = p1.columns
    psd = prior_sds(d11_models(), cols)
    worst = 0.0
    for cv in (0.01, 0.05, 0.10, 0.20)
        b = a2_bound(p1, runs, cols, psd; cv)
        for name in DESIGNS
            ref = m.bound[(cv, name)]
            e = max(maximum(abs.(b[name].ratio .- ref.ratio)), maximum(abs.(b[name].bracket .- ref.bracket)))
            worst = max(worst, e)
            @printf("σ_b %.0f%% %-17s ratio %s  bracket %s  max |Δ| %.2e\n", 100cv, name,
                    join((@sprintf("%.2f", v) for v in b[name].ratio), " "),
                    join((@sprintf("%.2f", v) for v in b[name].bracket), " "), e)
        end
    end
    @printf("a2check at %s (job %s): worst |Δ| against A2's merged.jls %.3e: %s\n", COMMIT, job(),
            worst, worst <= 1e-12 ? "PASS" : "FAIL")
    worst <= 1e-12 || exit(1)
end

function main(args)
    stage = args[1]
    if stage == "cells"
        run_cells(args[2], parse(Int, args[3]), parse(Int, args[4]), parse(Int, args[5]))
    elseif stage == "bound"
        run_bound(args[2], parse(Int, args[3]))
    elseif stage == "a2check"
        run_a2check(args[2])
    else
        error("unknown stage $stage")
    end
end

main(ARGS)
