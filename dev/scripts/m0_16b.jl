# spec/phases/16-recovery.md tasks 16b.1 and 16b.4, at the M0 rung (§12 "16a.11",
# approved 2026-10-22): three production-length chains of the three-block Gibbs
# sampler, at F10's settings (50 cells, 600 s, N = 20, 656 sweeps, 5 threads).
#
# Modes (RUN_MODE):
#   perturbed  16b.1: M0 replicate 16085 (V6's) with ENO's forward constant moved
#              by +1 prior SD, its reverse derived; every other draw is V6's
#   tx_only    16b.4: V6's dataset, transcripts only (σ_b held at its start)
#   met_only   16b.4: V6's dataset, bulk metabolites only
# Each starts from θ and σ_b drawn from their priors with Xoshiro(16_400_000 + i),
# i the mode's index, and from paths drawn by `gibbs_init`.
#
# Stages:
#   run <sweeps>               resumable: it checkpoints after every sweep
#   rates                      per-window update rates so far (the early check)
#   summary <v6dir> <dir>...   the posterior tables and F6, after burn-in
# Regenerate: m0_16b.sh, from a worktree at a pushed commit.

using InferCell
using Random
using Serialization
using Statistics
using Printf
using Distributions: params
using LinearAlgebra: BLAS
using MCMCChains: MCMCDiagnosticTools
const MDT = MCMCDiagnosticTools

BLAS.set_num_threads(1)

const MODES = (:perturbed, :tx_only, :met_only)
const MODE = Symbol(get(ENV, "RUN_MODE", "perturbed"))
MODE in MODES || error("RUN_MODE must be one of $MODES, not $MODE")
const DIR = get(ENV, "PILOT_DIR", joinpath(@__DIR__, "m0_16b_$(MODE)"))
mkpath(DIR)
const REPLICATE = 16085
const NCELLS = parse(Int, get(ENV, "PILOT_NCELLS", "50"))     # 3 for a smoke run
const HORIZON = parse(Float64, get(ENV, "PILOT_HORIZON", string(M0_HORIZON_S)))
const N = 20
const BURN = parse(Int, get(ENV, "BURN", "100"))       # 0 for a smoke run
const THRESHOLD = 0.10

"ENO's prior SD on the log scale, ln of its geometric SD (central_glycolysis.tsv)."
eno_prior_sd(models) = params(InferCell._parameter_by_name(models, :kcatF_R_ENO).prior)[2]

function dataset(mode)
    mode === :perturbed || return m0_dataset(REPLICATE; ncells = NCELLS, horizon = HORIZON)
    ms = m0_models()
    v6 = Dict(m0_dataset(REPLICATE; ncells = 1, horizon = 60.0).truth.values)
    v = v6[:kcatF_R_ENO] * exp(eno_prior_sd(ms))
    return m0_dataset(REPLICATE; ncells = NCELLS, horizon = HORIZON,
                      forwards = [:kcatF_R_ENO => v])
end

streams(mode) = mode === :tx_only ? [:transcripts] :
                mode === :met_only ? [:metabolites] : [:transcripts, :metabolites]

setup(mode = MODE) = (ds = dataset(mode);
                      (ds, gibbs_setup(ds; N, horizon = HORIZON, streams = streams(mode))))

chainfile(dir = DIR) = joinpath(dir, "chain.jls")

function cycle_seconds(g)
    run1() = (d = restore(g.pristine[1]); Random.seed!(1);
              for _ in 1:round(Int, HORIZON); handshake_step!(d); end)
    run1()
    return @elapsed run1()
end

function save(s, meta)
    tmp = chainfile() * ".tmp"
    serialize(tmp, (; state = gibbs_checkpoint(s), meta))
    mv(tmp, chainfile(); force = true)
end

function run_chain(sweeps)
    ds, g = setup()
    cyc = cycle_seconds(g)
    @printf("16b %s: M0 replicate %d, %d cells, streams %s, N = %d, %d threads; a plain cycle takes %.2f s\n",
            MODE, REPLICATE, NCELLS, g.streams, N, Threads.nthreads(), cyc)
    if isfile(chainfile())
        saved = deserialize(chainfile())
        saved.meta.mode === MODE || error("$(chainfile()) is a $(saved.meta.mode) chain")
        s = gibbs_resume(saved.state, g)
        meta = saved.meta
        println("  resumed at sweep $(s.sweep)")
    else
        rng = Xoshiro(16_400_000 + findfirst(==(MODE), MODES))
        t0 = draw_truth(rng, g.models; names = M0_TARGETS, purpose = :calibration)
        v = Dict(t0.values)
        σ0 = rand(rng, SIGMA_MET_PRIOR)
        cme0 = [v[n] for n in g.cme]
        u0 = [log(v[f]) for f in g.forwards]
        tinit = @elapsed s = gibbs_init(rng, g; cme = cme0, u = u0, σ = σ0)
        meta = (; mode = MODE, streams = g.streams, forwards = ds.forwards,
                start = (cme = cme0, u = u0, σ = σ0), init_seconds = tinit,
                names = vcat(g.cme, g.forwards, :sigma_b), truth = ds.truth.values,
                sigma_true = ds.sigma, ncells = NCELLS, N, replicate = REPLICATE)
        @printf("  start: cme = %s, u = %s, σ_b = %.3g; init took %.0f s\n",
                round.(cme0; sigdigits = 3), round.(u0; digits = 3), σ0, tinit)
    end
    meta = merge(meta, (; cycle = cyc, threads = Threads.nthreads()))
    save(s, meta)
    while s.sweep < sweeps
        gibbs_sweep!(s, g)
        save(s, meta)
        t = s.trace[end]
        @printf("  sweep %d: %.0f s (B3 %.0f, B1 %.0f, B2 %.2f, rebuild %.1f); ln θ = %s; σ_b = %.4f; windows changed %.2f\n",
                s.sweep, sum(t.seconds), t.seconds.b3, t.seconds.b1, t.seconds.b2,
                t.seconds.rebuild, round.(vcat(log.(t.cme), t.u); digits = 3), t.σ, mean(t.changed))
        flush(stdout)
    end
end

# Per cell-window update rates over every sweep so far. A window whose observed
# counts are 0 → 0 for every gene is reported apart, as the pilot's merge does;
# the dataset's own transcripts decide it in every mode.
function rates()
    r = deserialize(chainfile())
    ds, g = setup()
    T = size(g.bulk, 2)
    S = length(r.state.trace)
    rate = sum(t.changed for t in r.state.trace) ./ S
    tx = [Matrix{Int}(ds.transcripts[c, :, :]) for c in 1:r.meta.ncells]
    y0 = [g.pristine[c].jump.u[s] for s in g.tmap.states, c in 1:r.meta.ncells]
    live = falses(T, r.meta.ncells)
    for c in 1:r.meta.ncells, m in 1:T
        prev = m == 1 ? y0[:, c] : tx[c][:, m - 1]
        live[m, c] = any(prev .> 0) || any(tx[c][:, m] .> 0)
    end
    lv = rate[live]
    @printf("16b %s after %d sweeps: over %d live cell-windows the lowest update rate is %.2f and the median %.2f; %d are below %.2f. The %d empty ones average %.2f.\n",
            r.meta.mode, S, length(lv), minimum(lv), median(lv), count(<(THRESHOLD), lv),
            THRESHOLD, count(.!live), count(.!live) > 0 ? mean(rate[.!live]) : NaN)
end

# Kept draws, sweeps × parameters, on the log scale.
draws(trace, upto = length(trace)) =
    reduce(vcat, [permutedims(vcat(log.(t.cme), t.u, log(t.σ))) for t in trace[(BURN + 1):upto]])

# Bulk ESS per parameter of draws × chains × parameters.
bulk_ess(x::AbstractArray{<:Real, 3}) = MDT.ess(x; kind = :bulk)
chains3(xs) = cat([reshape(x, size(x, 1), 1, size(x, 2)) for x in xs]...; dims = 2)

function summary(v6dir, dirs)
    v6 = [deserialize(joinpath(v6dir, "chain$k.jls")) for k in 1:4]
    runs = [deserialize(chainfile(d)) for d in dirs]
    ms = m0_models()
    names = string.(first(runs).meta.names)
    P = length(names)
    prior_sd = vcat([params(InferCell._parameter_by_name(ms, Symbol(n)).prior)[2]
                     for n in names[1:(P - 1)]], params(SIGMA_MET_PRIOR)[2])
    string.(first(v6).meta.names) == names || error("V6's parameters differ from $names")
    # V6, the joint fit: each chain's sweeps 101 to 1,650, as in its merge.
    S = minimum(length(r.state.trace) for r in runs)
    joint = chains3([draws(r.state.trace) for r in v6])
    tv(meta) = vcat([log(Dict(meta.truth)[Symbol(n)]) for n in names[1:(P - 1)]],
                    log(meta.sigma_true))
    function rows(lab, x, truth)
        ess = bulk_ess(x)
        for p in 1:P
            if lab == "tx_only" && names[p] == "sigma_b"
                println("| $lab | sigma_b | held at its start, not inferred | — | — | — | — |")
                continue
            end
            v = vec(x[:, :, p])
            @printf("| %s | %s | %.3f ± %.3f | %.3f | %.0f | %.3f | %.2f |\n", lab, names[p],
                    mean(v), std(v), std(v) / prior_sd[p], ess[p], truth[p],
                    mean(v .<= truth[p]))
        end
    end
    println("Posteriors after $BURN sweeps of burn-in: sweeps $(BURN + 1) to $S of each new chain, ",
            "and $(BURN + 1) to $(length(first(v6).state.trace)) of V6's four (the joint fit). ",
            "SD / prior SD is the shrinkage. σ_b is held, not inferred, under transcripts only.\n")
    println("| Fit | Parameter (log) | Mean ± SD | SD / prior SD | Bulk ESS | Truth | Truth's quantile |")
    println("|---|---|---|---|---|---|---|")
    rows("joint (V6, 4 chains)", joint, tv(first(v6).meta))
    for r in runs
        rows(string(r.meta.mode), chains3([draws(r.state.trace, S)]), tv(r.meta))
    end
    p = findfirst(==("kcatF_R_ENO"), names)
    for r in runs
        r.meta.mode === :perturbed || continue
        x = chains3([draws(r.state.trace, S)])
        Δtruth = tv(r.meta)[p] - tv(first(v6).meta)[p]
        Δ = mean(x[:, :, p]) - mean(joint[:, :, p])
        se = hypot(std(x[:, :, p]) / sqrt(bulk_ess(x)[p]),
                   std(joint[:, :, p]) / sqrt(bulk_ess(joint)[p]))
        @printf("\n16b.1: ENO's truth moved by %+.4f and its posterior mean by %+.4f (MC SE %.4f), a ratio of %.3f ± %.3f. Shrinkage %.3f; the moved truth's quantile %.3f; V6's truth's quantile %.3f.\n",
                Δtruth, Δ, se, Δ / Δtruth, se / Δtruth, std(x[:, :, p]) / prior_sd[p],
                mean(x[:, :, p] .<= tv(r.meta)[p]), mean(x[:, :, p] .<= tv(first(v6).meta)[p]))
    end
end

if ARGS[1] == "run"
    run_chain(parse(Int, ARGS[2]))
elseif ARGS[1] == "rates"
    rates()
elseif ARGS[1] == "summary"
    summary(ARGS[2], ARGS[3:end])
else
    error("stage must be run, rates or summary")
end
