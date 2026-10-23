# Issue #89 B5: the mixed likelihood end to end on one M0 dataset
# (spec/phases/approximate-method.md §12 2026-10-23, "B5 as planned").
#
# On m0_v6, V6's replicate 16085:
#   1. theta_ref: the MAP at δ = 0 and S = 0 (pass 1), and the cost of one gradient;
#   2. ensemble:  the reference composition at θ_ref, 800 cells from the dataset's
#                 own seed range, in chunks;
#   3. nuts:      one NUTS chain on m(θ) + δ with S, from θ_ref;
#   4. merge:     δ, S's spectrum, R̂, ESS, divergences, tree depth, the cost per
#                 gradient against KC's bound, and, as information only, the
#                 marginals beside F10's.
# No pass line on the posterior: AV6 (B7) is the test.
#
# Usage: julia --project dev/scripts/mixed_nuts_89.jl theta_ref <outdir>
#        julia --project dev/scripts/mixed_nuts_89.jl ensemble <outdir> <chunk>
#        julia --project dev/scripts/mixed_nuts_89.jl nuts <outdir> <chain>
#        julia --project dev/scripts/mixed_nuts_89.jl merge <outdir>
# Regenerate: dev/scripts/mixed_nuts_89.sh <outdir outside the worktrees>.
# SMOKE=1: 3 cells, a 120 s horizon, 8 reference cells, 10 + 10 NUTS steps.

using InferCell
using Random
using Serialization
using Statistics
using LinearAlgebra
using Printf
using Dates
using ForwardDiff
using Turing: sample, NUTS
using MCMCChains

const SMOKE = get(ENV, "SMOKE", "0") == "1"
const COMMIT = strip(read(`git rev-parse --short HEAD`, String))
const REPLICATE = 16_085
const NCELLS = SMOKE ? 3 : 50
const HORIZON = SMOKE ? 120.0 : M0_HORIZON_S
const NREF = SMOKE ? 8 : 800
const NCHUNK = SMOKE ? 2 : 16
const NADAPT = parse(Int, get(ENV, "NADAPT", SMOKE ? "10" : "300"))
const NDRAW = parse(Int, get(ENV, "NDRAW", SMOKE ? "10" : "300"))
const TARGET_ACCEPT = 0.8
# KC: a tenth of F10's production posterior, 507.5 CPU-h (phase16_ledger.md, 18751635).
const KC_BOUND_CPUH = 50.75
const F10_DRAWS = joinpath(@__DIR__, "f10_16a10_draws.tsv")
slurm_job() = get(ENV, "SLURM_ARRAY_JOB_ID", get(ENV, "SLURM_JOB_ID", "none"))

ds_file(out) = joinpath(out, "dataset.jls")
ref_file(out) = joinpath(out, "theta_ref.jls")
chunk_file(out, k) = joinpath(out, "ensemble_$k.jls")
chain_file(out, c) = joinpath(out, "chain_$c.jls")

dataset(out) = deserialize(ds_file(out))

function run_theta_ref(out)
    mkpath(out)
    t_ds = @elapsed ds = m0_dataset(REPLICATE; ncells = NCELLS, horizon = HORIZON)
    serialize(ds_file(out), ds)
    s = mixed_setup(ds)
    x0 = mixed_prior_median(s)
    ℓ(z) = mixed_logdensity(s, z)
    ForwardDiff.gradient(ℓ, x0)                                   # compiles
    grad = [@timed(ForwardDiff.gradient(ℓ, x0)) for _ in 1:3]
    dens = [@timed(ℓ(x0)) for _ in 1:3]
    t_map = @elapsed ref = reference_theta(s)
    values = mixed_theta_values(s, ref.x)
    sref = InferCell._mixed_solve(s, exp.(ref.x[1:length(s.names)]))
    m_ref = mean_field_log_bulk(s, sref.sol)
    serialize(ref_file(out),
              (; ref..., values, m_ref, names = mixed_parameter_names(s), commit = COMMIT,
               smoke = SMOKE, job = slurm_job(), dataset_seconds = t_ds, map_seconds = t_map,
               grad_seconds = [g.time for g in grad], grad_bytes = [g.bytes for g in grad],
               density_seconds = [d.time for d in dens], density_bytes = [d.bytes for d in dens],
               truth = ds.truth.values, sigma = ds.sigma))
    @printf("%s theta_ref: %s in %.0f s, log density %.3f; gradient %.2f s, %.2f GiB\n",
            now(), ref.method, t_map, ref.logdensity, median(g.time for g in grad),
            median(g.bytes for g in grad) / 2^30)
end

function run_ensemble(out, k)
    ds = dataset(out)
    ref = deserialize(ref_file(out))
    seeds = reference_seeds(ds, NREF)
    per = cld(NREF, NCHUNK)
    mine = seeds[((k - 1) * per + 1):min(k * per, NREF)]
    t = @elapsed cells = reference_ensemble(ref.values, mine; horizon = HORIZON)
    serialize(chunk_file(out, k), (; k, seeds = mine, cells, commit = COMMIT, job = slurm_job(),
                                   seconds = t))
    @printf("%s ensemble chunk %d: %d cells in %.0f s\n", now(), k, length(mine), t)
end

function reference(out)
    ds = dataset(out)
    ref = deserialize(ref_file(out))
    chunks = [deserialize(chunk_file(out, k)) for k in 1:NCHUNK]
    cells = cat((c.cells for c in chunks)...; dims = 1)
    seeds = reduce(vcat, (c.seeds for c in chunks))
    seeds == reference_seeds(ds, NREF) || error("the ensemble's chunks are not the reference seeds")
    d = delta_and_S(cells, ref.m_ref; C = NCELLS)
    return (; ds, ref, cells, d..., ensemble_seconds = sum(c.seconds for c in chunks))
end

function run_nuts(out, c)
    r = reference(out)
    s = with_reference(mixed_setup(r.ds), r.delta, r.S)
    model = mixed_model(s)
    t = @elapsed chain = sample(Xoshiro(1000 + c), model, NUTS(NADAPT, TARGET_ACCEPT), NDRAW;
                                initial_params = r.ref.x, progress = false)
    serialize(chain_file(out, c), (; c, chain, seconds = t, nadapt = NADAPT, ndraw = NDRAW,
                                   commit = COMMIT, job = slurm_job()))
    @printf("%s chain %d: %d + %d steps in %.0f s\n", now(), c, NADAPT, NDRAW, t)
end

# F10's reference: chains 1 to 4, sweeps 101 to 1,650, on the log scale
# (f10_16a10_result.md).
function f10_draws()
    lines = readlines(F10_DRAWS)
    head = split(lines[1], '\t')
    rows = [parse.(Float64, split(l, '\t')) for l in lines[2:end]]
    filter!(r -> r[1] in 1:4 && r[2] > 100, rows)
    return Dict(Symbol(h) => [r[j] for r in rows] for (j, h) in enumerate(head))
end

function run_merge(out)
    r = reference(out)
    files = sort(filter(f -> startswith(f, "chain_"), readdir(out)))
    runs = [deserialize(joinpath(out, f)) for f in files]
    chains = chainscat((x.chain for x in runs)...)
    names = r.ref.names
    X = Array(chains[:, [Symbol("x[$i]") for i in eachindex(names)], :])   # draws × params × chains
    io = IOBuffer()
    p(args...) = (println(io, args...); println(args...))
    p("# B5 smoke: the mixed likelihood on m0_v6 (#89)\n")
    p("Commit `$(COMMIT)` (merge job $(slurm_job())); stages at `$(r.ref.commit)`. ",
      "SMOKE = $(SMOKE). $(NCELLS) cells, $(HORIZON) s, $(NREF) reference cells, ",
      "$(length(runs)) chains of $(NADAPT) adaptation + $(NDRAW) draws, target ",
      "acceptance $(TARGET_ACCEPT). No pass line on the posterior (AV6, B7).\n")
    p("## The reference θ\n")
    p("Method: `$(r.ref.method)`, log density $(round(r.ref.logdensity; digits = 3)) ",
      "($(round(r.ref.map_seconds; digits = 0)) s).\n")
    p("| coordinate | θ_ref (log) | truth (log) |")
    p("|---|---|---|")
    tv = Dict(r.ref.truth)
    for (i, n) in enumerate(names)
        tr = n === :log_sigma_b ? log(r.ref.sigma) : log(tv[Symbol(String(n)[5:end])])
        p(@sprintf("| %s | %.4f | %.4f |", n, r.ref.x[i], tr))
    end
    p("\n## δ and S\n")
    P = length(M0_PANEL)
    Δ = reshape(r.delta, P, :)
    p("δ per pool over the saves: mean, min, max (nats).\n")
    p("| pool | mean | min | max |")
    p("|---|---|---|---|")
    for (j, o) in enumerate(M0_PANEL)
        p(@sprintf("| %s | %+.4f | %+.4f | %+.4f |", o, mean(Δ[j, :]), minimum(Δ[j, :]),
                   maximum(Δ[j, :])))
    end
    ev = sort(eigvals(Symmetric(r.S)); rev = true)
    p(@sprintf("\nS (C = %d): trace %.4g; largest eigenvalues %s; smallest %.3g; ",
               NCELLS, tr_(r.S), join((@sprintf("%.3g", e) for e in ev[1:5]), ", "), ev[end]),
      @sprintf("σ_b² at θ_ref %.3g.", exp(2 * r.ref.x[end])))
    p("\n## The chains\n")
    ess_ = ess(chains)
    rh = rhat(chains)
    p("| coordinate | mean | sd | ESS (bulk) | R̂ |")
    p("|---|---|---|---|---|")
    for (i, n) in enumerate(names)
        v = vec(X[:, i, :])
        sym = Symbol("x[$i]")
        p(@sprintf("| %s | %.4f | %.4f | %.0f | %.3f |", n, mean(v), std(v),
                   ess_[sym, :ess], rh[sym, :rhat]))
    end
    nerr = sum(chains[:numerical_error])
    p(@sprintf("\nDivergent or non-finite steps: %d of %d draws. Tree depth: median %.0f, max %.0f. ",
               nerr, length(chains[:numerical_error]), median(chains[:tree_depth]),
               maximum(chains[:tree_depth])),
      @sprintf("Leapfrog steps per draw: median %.0f.", median(chains[:n_steps])))
    g = median(r.ref.grad_seconds)
    steps = sum(chains[:n_steps])
    cpu_chains = sum(x.seconds for x in runs) / 3600
    p("\n## Cost (KC)\n")
    p(@sprintf("One gradient: %.2f s, %.2f GiB allocated (median of 3, at the prior median). ",
               g, median(r.ref.grad_bytes) / 2^30),
      @sprintf("One density: %.2f s. ", median(r.ref.density_seconds)),
      @sprintf("The chains took %.2f CPU-h for %d leapfrog steps in their draws.", cpu_chains, steps))
    p(@sprintf("The reference ensemble took %.2f CPU-h, the MAP %.2f.", r.ensemble_seconds / 3600,
               r.ref.map_seconds / 3600))
    total = cpu_chains + r.ensemble_seconds / 3600 + r.ref.map_seconds / 3600
    p(@sprintf("Total %.2f CPU-h against KC's bound of %.2f (a tenth of F10's 507.5). ",
               total, KC_BOUND_CPUH),
      "A smoke's chains are short, so this is a lower bound on a production fit.")
    p("\n## Beside F10's exact posterior (information only; its reference chains 1–4, sweeps 101–1,650)\n")
    f = f10_draws()
    p("| coordinate | mixed: mean ± sd | F10: mean ± sd | difference / F10 sd |")
    p("|---|---|---|---|")
    for (i, n) in enumerate(names)
        key = n === :log_sigma_b ? :sigma_b : Symbol(String(n)[5:end])
        haskey(f, key) || continue
        v = vec(X[:, i, :])
        a, b = mean(f[key]), std(f[key])
        p(@sprintf("| %s | %.4f ± %.4f | %.4f ± %.4f | %+.2f |", n, mean(v), std(v), a, b,
                   (mean(v) - a) / b))
    end
    write(joinpath(out, "summary.md"), String(take!(io)))
end

tr_(A) = sum(diag(A))

function main()
    stage = ARGS[1]
    out = ARGS[2]
    stage == "theta_ref" ? run_theta_ref(out) :
    stage == "ensemble" ? run_ensemble(out, parse(Int, ARGS[3])) :
    stage == "nuts" ? run_nuts(out, parse(Int, ARGS[3])) :
    stage == "merge" ? run_merge(out) : error("unknown stage $stage")
end

main()
