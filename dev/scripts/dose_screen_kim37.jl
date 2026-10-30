# KIM-37: the ptsG knockdown dose screen (M4). It picks the dose by a declared
# criterion on the validity of the regime: the pools the rebuilt rate constants
# read must stay away from K_d, since the one-way reduction rests on
# elasticities of 0.04–0.09. It never computes the effect ΔQ, its size or its
# coverage. θ comes from fresh prior draws, the population the SBC truths came
# from, with every seed from M4's reserved block 18_437_xxx; no evidence dataset
# or posterior is read.
#
# Usage (dev/scripts/dose_screen_kim37.sh submits every stage):
#   julia --project dev/scripts/dose_screen_kim37.jl <stage> <outdir>
# Stages:
#   m0      M0 at 600 s: 50 prior draws of M0's five targets, 4 cells each
#   d11     the reduced model at 6,300 s: 20 prior draws of the five targets and
#           the 15 held promoters, 2 cells each (for a reduced-model intervention,
#           which the O2 claim names; the freeze decides whether it runs)
#   merge   the criterion and the selection, below
# Each cell runs at doses 0, 0.1, 0.2, 0.3, 0.4 and 0.5 (`Knockdown`, ptsG), on
# the same seed at every dose. The driver is stepped here, as `predict`'s hybrid
# cell steps it, on the metered composition; after every 60th handshake (the
# rebuild cadence) it records, for each watched pool, its level in mM and the
# largest elasticity d ln k / d ln pool over the rebuilt rate constants that read
# it (`rate_constant_elasticity`'s central difference, pool by pool, so a pool at
# zero is its own category rather than an error).
#
# Watched pools and their K_d (spec §12 2026-10-29; transcription.jl, translation.jl):
#   :M_atp_c, :M_gtp_c   transcription, K_d = 0.1 mM
#   :M_trna_chg_c        translation, per-amino-acid share against K_d = 1e-3 mM,
#                        so 0.02 mM on the lumped pool
# CTP and UTP are chemostats, which a knockdown cannot move.
#
# The criterion, declared before the run (M4 plan, approved 2026-10-30). A dose
# is admissible when, for each watched pool,
#   (1) the 95th percentile, over draws × cells × instants, of the pool's largest
#       elasticity is at most max(0.1, 1.1 × its dose-0 value), and
#   (2) the fraction of instants with the pool below its K_d (a pool at zero
#       counts) is at most dose 0's fraction plus 0.01.
# The frozen dose is the largest admissible dose in {0.1, 0.2, 0.3}, judged on
# M0, where the panel is scored. Doses 0.4 and 0.5 only show where the regime
# ends. If no dose up to 0.3 is admissible, the panel is dropped (the cut line)
# and §12 records it. The reduced model is judged the same way, for information,
# with and without the draws whose baseline cells starve by 14c's rule (a counter
# clipped at more than half its drains), since V would reject those before any
# intervention.
#
# Seeds: draw i's parameters come from Xoshiro(18_437_100 + i) on M0 (i = 1…50)
# and Xoshiro(18_437_500 + i) on the reduced model (i = 1…20). Cell c of draw i
# has seed 1000 × that base + 100i + c: 18_437_100_101 to 18_437_105_004 on M0,
# and 18_437_500_101 to 18_437_502_002 on the reduced model.

using InferCell
using InferCell: AbstractSubModel
using Random
using Statistics: mean, quantile
using Printf
using Dates

const STAGE = Symbol(ARGS[1])
const OUT = ARGS[2]
const COMMIT = strip(read(`git rev-parse --short HEAD`, String))
const JOB = get(ENV, "SLURM_JOB_ID", "none")
const DOSES = (0.0, 0.1, 0.2, 0.3, 0.4, 0.5)
const WATCHED = (:M_atp_c, :M_gtp_c, :M_trna_chg_c)
const KD = Dict(:M_atp_c => 0.1, :M_gtp_c => 0.1, :M_trna_chg_c => 0.02)
const CANDIDATES = (0.1, 0.2, 0.3)
const SMOKE = get(ENV, "SMOKE", "0") == "1"      # a plumbing pass at tiny sizes, not a result

say(x...) = (println(now(), " ", x...); flush(stdout))

# Each watched pool's level (mM) and its largest elasticity over the rate
# constants rebuilt from it, at the driver's current pools; NaN at a pool ≤ 0.
function watched(d)
    level = Dict{Symbol, Float64}()
    eps = Dict{Symbol, Float64}(s => -Inf for s in WATCHED)
    for r in d.rebuilds
        pools = InferCell._rebuild_pools(d, r)
        p_local = collect(Float64, view(d.jump.p, r.p_idxs))
        k_at(v) = InferCell.rate_constants(p_local, d.ode.t, r.model, InferCell._fvec(v))
        base = collect(Float64, k_at(pools))
        for (j, s) in enumerate(r.species)
            s in WATCHED || continue
            x = pools[j]
            level[s] = x
            if x > 0 && all(>(0), base)
                up, dn = copy(pools), copy(pools)
                up[j], dn[j] = x * (1 + 1e-4), x * (1 - 1e-4)
                ku, kd = k_at(up), k_at(dn)
                e = maximum((log(ku[i]) - log(kd[i])) / (log(up[j]) - log(dn[j]))
                            for i in eachindex(base))
                eps[s] = isnan(eps[s]) ? NaN : max(eps[s], e)
            else
                eps[s] = NaN
            end
        end
    end
    return level, eps
end

# One cell at one dose: the watched pools at every 60th handshake, whether it
# finished, and whether it starved by 14c's rule.
function screen_cell(ms, vals, seed, H)
    Random.seed!(seed)
    d = build_problem(ms; tspan = (0.0, H), complete = true)
    set_parameters!(d, ms, vals)
    Random.seed!(seed)
    rec = ClipRecord(d)
    rows = NamedTuple[]
    failure = :none
    try
        for k in 1:round(Int, H / d.interval)
            handshake_step!(d)
            record_clips!(rec, d)
            if k % 60 == 0
                level, eps = watched(d)
                for s in WATCHED
                    push!(rows, (; t = d.ode.t, pool = s, mM = level[s], eps = eps[s]))
                end
            end
        end
    catch e
        e isa MetabolicBlockStalled || rethrow()
        failure = :not_advanced
    end
    starved = any(r -> r.drains > 0.5 * rec.drains, clip_summary(rec))
    return (; rows, failure, starved)
end

function screen(label, models, names, n_draws, n_cells, H, seed0)
    ms = metered(models)
    draws = [draw_parameters(Xoshiro(seed0 + i), ms, names) for i in 1:n_draws]
    jobs = [(i, c, dose) for i in 1:n_draws for c in 1:n_cells for dose in DOSES]
    res = Vector{Any}(undef, length(jobs))
    done = Threads.Atomic{Int}(0)
    InferCell._restoring_global_rng() do
        InferCell._each(length(jobs), Threads.nthreads() > 1) do j
            i, c, dose = jobs[j]
            vals = intervene(Knockdown(; dose), ms, draws[i])
            res[j] = screen_cell(ms, vals, seed0 * 1000 + 100i + c, H)
            n = Threads.atomic_add!(done, 1) + 1
            n % 20 == 0 && say("  $label: $n of $(length(jobs)) cells")
        end
    end
    open(joinpath(OUT, "$(label)_cells.tsv"), "w") do io
        println(io, "# commit $COMMIT, job $JOB, $(now()), H = $H, seeds from $seed0")
        println(io, "draw\tcell\tdose\tfailure\tstarved")
        for (j, (i, c, dose)) in enumerate(jobs)
            println(io, "$i\t$c\t$dose\t$(res[j].failure)\t$(res[j].starved)")
        end
    end
    open(joinpath(OUT, "$(label)_pools.tsv"), "w") do io
        println(io, "# commit $COMMIT, job $JOB, $(now()), H = $H, seeds from $seed0")
        println(io, "draw\tcell\tdose\tt\tpool\tmM\teps")
        for (j, (i, c, dose)) in enumerate(jobs), r in res[j].rows
            println(io, "$i\t$c\t$dose\t$(r.t)\t$(r.pool)\t$(r.mM)\t$(r.eps)")
        end
    end
    say("$label: $(count(r -> r.failure !== :none, res)) failed, ",
        "$(count(r -> r.starved, res)) starved, of $(length(res))")
end

m0() = screen("m0", m0_models(), M0_TARGETS, SMOKE ? 2 : 50, SMOKE ? 2 : 4,
              SMOKE ? 120.0 : 600.0, 18_437_100)

function d11()
    dm = metered(d11_models())
    held = held_parameters(dm, M0_TARGETS).held
    screen("d11", d11_models(), vcat(M0_TARGETS, held), SMOKE ? 2 : 20, 2,
           SMOKE ? 120.0 : 6300.0, 18_437_500)
end

function readtsv(path)
    lines = [l for l in eachline(path) if !startswith(l, "#")]
    head = Symbol.(split(lines[1], '\t'))
    return [NamedTuple{Tuple(head)}(Tuple(split(l, '\t'))) for l in lines[2:end]]
end

function judge(io, label; exclude_starved = false)
    cells = readtsv(joinpath(OUT, "$(label)_cells.tsv"))
    starved = Set((c.draw, c.cell) for c in cells if c.dose == "0.0" && c.starved == "true")
    bad_draws = Set(d for (d, _) in starved)
    keep(r) = !(exclude_starved && r.draw in bad_draws)
    rows = [r for r in readtsv(joinpath(OUT, "$(label)_pools.tsv")) if keep(r)]
    stat = Dict{Tuple{Float64, Symbol}, NamedTuple}()
    for dose in DOSES, s in WATCHED
        rs = [r for r in rows if parse(Float64, r.dose) == dose && Symbol(r.pool) === s]
        e = [parse(Float64, r.eps) for r in rs]
        x = [parse(Float64, r.mM) for r in rs]
        fin = filter(isfinite, e)
        stat[(dose, s)] = (; e95 = isempty(fin) ? NaN : quantile(fin, 0.95),
                           below = mean(x .< KD[s]), zero = mean(x .<= 0), n = length(rs))
    end
    nfail = Dict(dose => count(c -> parse(Float64, c.dose) == dose && c.failure != "none" &&
                                    keep(c), cells) for dose in DOSES)
    tag = exclude_starved ? " (without the $(length(bad_draws)) draws that starve at baseline)" : ""
    println(io, "\n## $label$tag\n")
    println(io, "| dose | pool | e95 | limit | below K_d | limit | zero | admissible |")
    println(io, "|---|---|---|---|---|---|---|---|")
    ok = Dict{Float64, Bool}()
    for dose in DOSES
        ok[dose] = true
        for s in WATCHED
            b, z = stat[(0.0, s)], stat[(dose, s)]
            le, lb = max(0.1, 1.1b.e95), b.below + 0.01
            a = z.e95 <= le && z.below <= lb
            ok[dose] &= a
            println(io, @sprintf("| %.1f | %s | %.4f | %.4f | %.4f | %.4f | %.4f | %s |", dose, s,
                                 z.e95, le, z.below, lb, z.zero, a))
        end
    end
    println(io, "\nFailed cells by dose: ", join(["$d: $(nfail[d])" for d in DOSES], ", "), ".")
    adm = [d for d in CANDIDATES if ok[d]]
    pick = isempty(adm) ? nothing : maximum(adm)
    println(io, "Admissible doses: ", [d for d in DOSES if ok[d] && d > 0], ". ",
            "Selected (largest admissible in $(CANDIDATES)): ",
            pick === nothing ? "none, so the panel is dropped" : pick, ".")
    return pick
end

function merge()
    open(joinpath(OUT, "dose_screen.md"), "w") do io
        println(io, "# KIM-37 dose screen at $COMMIT, job $JOB, $(now())")
        pick = judge(io, "m0")
        judge(io, "d11")
        judge(io, "d11"; exclude_starved = true)
        println(io, "\n**The rule's pick, on M0:** ", pick === nothing ? "none" : pick)
    end
    print(read(joinpath(OUT, "dose_screen.md"), String))
end

function main()
    mkpath(OUT)
    say("KIM-37 dose screen, stage $STAGE, at $COMMIT, job $JOB, $(Threads.nthreads()) thread(s)")
    STAGE === :m0 ? m0() : STAGE === :d11 ? d11() : STAGE === :merge ? merge() :
    error("unknown stage $STAGE")
    say("done")
end

main()
