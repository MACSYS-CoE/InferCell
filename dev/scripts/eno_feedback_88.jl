# Issue #88 task A3: what does reducing the reverse feedback (pools → rebuilt
# rate constants) cost, and does saturation or the 60 s cadence explain the ENO
# null of `eno_path_information_16a7.jl`?
#
# Four arms on the reduced syn3A model (published drains, `d11_models`):
# - `live60`: the rebuild reads the live pools every 60 s (the published model);
# - `nominal`: it reads the build's initial ATP, GTP and charged tRNA, so the
#   reverse edge is cut at fixed reference pools (`override_rebuild_pools!`);
# - `popmean`: it reads the population mean of those pools at the truth, from
#   `live60`'s 200-cell ensemble, a fixed trajectory. In the ensembles, run at the
#   truth, that is mean field's reverse edge: only the cell-to-cell part of the
#   feedback is dropped (note §7). In the path replays at ln θ ± δ the pools stay
#   frozen at the truth's mean, where a mean-field loop would move them with θ, so
#   there it is a frozen reference like `nominal`;
# - `live6`: live pools, rebuilt every 6 s (spec R11's comparison).
# CTP and UTP are chemostats, read as before in every arm.
#
# Two scores, at each of two truths (15.7's `Xoshiro(1507)` and `Xoshiro(8802)`):
# 1. Path information, as the 16a7 script: one cell's full jump path recorded
#    and replayed at ln θ ± δ (δ ∈ {0.05, 0.2}) for ENO, FBA and D11's three
#    promoters; complete-data Fisher information of 200 cells, expected and
#    observed, as posterior over prior SD. Cells are seeds 30000 + k.
#    `live60` against `live6` tests the cadence; `nominal` and `popmean` against
#    `live60` give frozen pools' cost for the promoters (ENO and FBA reach the
#    path only through the pools, so under frozen pools their ratio should be 1,
#    a check; a mean-field path would keep a θ → mean pool → rate constant
#    channel that neither arm has).
# 2. The bulk metabolite likelihood: the 200-cell mean of each M0_PANEL pool at
#    each 60 s save, per arm, against `live60`'s. Rebuilds change how the jump
#    block consumes random numbers, so arms are not paired path by path; the
#    score is distributional: Σ_t (log x̄_arm − log x̄_live60)² / (2σ_b²), the KL
#    between the two arms' bulk likelihoods, beside the same sum between two
#    halves of `live60`'s cells (its sampling noise). Cells are seeds 40000 + k.
#
# Usage: julia --project dev/scripts/eno_feedback_88.jl ensemble <arm> <truth> <task> <ntasks> <outdir>
#        julia --project dev/scripts/eno_feedback_88.jl path <arm> <truth> <cell> <outdir>
#        julia --project dev/scripts/eno_feedback_88.jl merge <outdir>
# Regenerate: dev/scripts/eno_feedback_88.sh <outdir>

using InferCell
using Random
using Serialization
using Statistics
using Printf
using LinearAlgebra
using Dates
using Distributions: params

const ARMS = (:live60, :nominal, :popmean, :live6)
const TRUTHS = Dict(1 => 1507, 2 => 8802)
const FW = vcat([:kcatF_R_ENO, :kcatF_R_FBA], D11_TARGETS[1:3])
const SHORT = ["ENO", "FBA", "S_0607 (GAPD)", "S_0445 (PGI)", "S_0779 (ptsG)"]
const DELTAS = [0.05, 0.2]
const NCELLS = 200
const NENSEMBLE = 200
const FED = (:M_atp_c, :M_gtp_c, :M_trna_chg_c)   # the live pools a rebuild reads
const SMOKE = get(ENV, "SMOKE", "0") == "1"
const HORIZON = SMOKE ? 120.0 : COREA_CYCLE_S
const COMMIT = strip(read(`git rev-parse --short HEAD`, String))

function models(arm)
    ms = d11_models()
    arm === :live6 || return ms
    drain = (clip = :clamped_deficit_carried, smoothing = nothing)
    i = findfirst(m -> m isa CoreATranscription, ms)
    j = findfirst(m -> m isa CoreATranslation, ms)
    ms[i] = CoreATranscription(; drain..., interval = 6.0)
    ms[j] = CoreATranslation(; drain..., interval = 6.0)
    return ms
end

truth_of(ms, t) = draw_truth(Xoshiro(TRUTHS[t]), ms; purpose = :recovery).values

ensemble_dir(out, arm, t) = joinpath(out, "ensemble_$(arm)_truth$(t)")
path_dir(out, arm, t) = joinpath(out, "path_$(arm)_truth$(t)")

# An arm's ensemble: every cell's matrix, the species and the save times.
function load_ensemble(out, arm, t)
    ps = load_parts(ensemble_dir(out, arm, t))
    return vcat((p.cells for p in ps)...), first(ps).species, first(ps).times
end

# The population-mean fed pools (mM) at each 60 s save, from live60's ensemble.
const MEAN_POOLS = Dict{Tuple{String, Int}, Any}()
mean_pools(out, t) = get!(() -> _mean_pools(out, t), MEAN_POOLS, (out, t))
function _mean_pools(out, t)
    cells, species, times = load_ensemble(out, :live60, t)
    v = findfirst(==(:volume_litres), species)
    pools = Dict{Symbol, Vector{Float64}}()
    for s in FED
        i = findfirst(==(s), species)
        pools[s] = vec(mean(reduce(vcat, [permutedims(c[i, :] ./ particles_per_mM.(c[v, :]))
                                          for c in cells]); dims = 1))
    end
    return times, pools
end

# A fresh driver at the truth, with the arm's rebuild in place.
function fresh(arm, ms, truth, out, t)
    d = build_problem(ms; tspan = (0.0, HORIZON), complete = true)
    set_parameters!(d, ms, truth)
    if arm === :nominal
        ref = Dict(s => d.ode.u[k] for r in d.rebuilds for (s, k) in zip(r.species, r.pool_idxs)
                   if k != 0)
        Set(keys(ref)) == Set(FED) || error("the rebuilds read $(keys(ref)), not $FED")
        override_rebuild_pools!(d, (s, _) -> get(ref, s, nothing))
    elseif arm === :popmean
        times, pools = mean_pools(out, t)
        override_rebuild_pools!(d, (s, τ) -> haskey(pools, s) ?
                                pools[s][argmin(abs.(times .- τ))] : nothing)
    end
    return d
end

function run_ensemble(arm, t, task, ntask, out)
    arm in ARMS || error("arm must be one of $ARMS")
    dir = ensemble_dir(out, arm, t); mkpath(dir)
    ms = models(arm)
    truth = truth_of(ms, t)
    cells = Matrix{Float64}[]; seeds = Int[]
    species, times, labels = Symbol[], Float64[], Symbol[]
    for k in 1:(SMOKE ? 4 : NENSEMBLE)
        (k - 1) % ntask == task || continue
        seed = 40_000 + k
        println("$(now()) ensemble $arm truth $t seed $seed"); flush(stdout)
        Random.seed!(seed)
        d = fresh(arm, ms, truth, out, t)
        Random.seed!(seed)
        r = emit_observables!(d, ms; every = 60)
        labels = [l.category for l in driver_declarations(d)]
        keep = [i for (i, s) in enumerate(r.species) if s in Set(vcat(M0_PANEL, collect(FED)))]
        species, times = vcat(r.species[keep], [:volume_litres]), r.times
        push!(cells, vcat(r.counts[keep, :], permutedims(r.volume))); push!(seeds, seed)
    end
    serialize(joinpath(dir, "task_$(task).jls"),
              (; arm, truth = t, values = truth, commit = COMMIT, smoke = SMOKE, job = slurm_job(),
               species, times, seeds, cells, labels))
end

slurm_job() = get(ENV, "SLURM_ARRAY_JOB_ID", get(ENV, "SLURM_JOB_ID", "none"))

function run_path(arm, t, cell, out)
    arm in ARMS || error("arm must be one of $ARMS")
    dir = path_dir(out, arm, t); mkpath(dir)
    ms = models(arm)
    truth = truth_of(ms, t)
    seed = 30_000 + cell
    Random.seed!(seed)
    base = snapshot(fresh(arm, ms, truth, out, t))
    d = restore(base)
    Random.seed!(seed)
    record_path!(d)
    for _ in 1:round(Int, HORIZON)
        handshake_step!(d)
    end
    path = recorded_path(d)
    u0 = [log(Dict(truth)[f]) for f in FW]
    function lp(i, x)
        e = restore(base)
        set_parameters!(e, ms, derived_ode_values(ms, [FW[i] => exp(x)]))
        return path_logdensity(replay!(e, path, round(Int, HORIZON); density = true))
    end
    e0 = restore(base)
    lp0 = path_logdensity(replay!(e0, path, round(Int, HORIZON); density = true))
    res = Dict{Tuple{Int, Float64, Int}, Float64}()
    for i in eachindex(FW), δ in DELTAS, s in (-1, 1)
        res[(i, δ, s)] = lp(i, u0[i] + s * δ)
        @printf("%s cell %d: %s %+g → Δ log p = %+.4f\n", now(), seed, FW[i], s * δ,
                res[(i, δ, s)] - lp0)
        flush(stdout)
    end
    serialize(joinpath(dir, "cell_$(cell).jls"),
              (; arm, truth = t, values = truth, seed, lp0, out = res, events = length(path.times),
               commit = COMMIT, smoke = SMOKE, job = slurm_job()))
end

function load_parts(dir)
    files = sort(filter(f -> occursin(r"^(task|cell)_\d+\.jls$", f), readdir(dir)))
    isempty(files) && error("no files in $dir")
    return [deserialize(joinpath(dir, f)) for f in files]
end

# The mean over `cells` of each panel pool at each save after t = 0, floored at 1.
function bulk(cells, species, times)
    live = findall(>(0), times)
    [log(max(mean(c[i, j] for c in cells), 1.0)) for i in [findfirst(==(o), species)
                                                         for o in M0_PANEL], j in live]
end

# Every input file must come from one commit and one kind of run, so a moved
# HEAD or a smoke file left in `out` cannot be merged silently.
function check_provenance(out)
    dirs = filter(isdir, [f(out, a, t) for f in (ensemble_dir, path_dir) for a in ARMS
                          for t in keys(TRUTHS)])
    parts = [x for d in dirs for x in load_parts(d)]
    commits, smokes = unique(x.commit for x in parts), unique(x.smoke for x in parts)
    length(commits) == 1 || error("the input files come from commits $commits")
    length(smokes) == 1 || error("the input files mix smoke and full runs")
    jobs = Dict(basename(d) => only(unique(x.job for x in load_parts(d))) for d in dirs)
    return only(commits), only(smokes), jobs
end

function merge(out)
    io = IOBuffer()
    p(args...) = (println(io, args...); println(args...))
    ms = d11_models()
    s0 = [params(InferCell._parameter_by_name(ms, f).prior)[2] for f in FW]
    # Posterior over prior SD from the information I about the parameters `ix`.
    ratio(I, ix = eachindex(FW)) = sqrt.(diag(inv(I[ix, ix] + Diagonal(1 ./ s0[ix] .^ 2)))) ./ s0[ix]
    # From one parameter's own information alone.
    alone(c, i) = c >= 0 ? 1 / sqrt(1 + NCELLS * c * s0[i]^2) : NaN
    rng = Xoshiro(88_003)
    commit, smoke, jobs = check_provenance(out)
    p("# Issue #88 A3: the cost of reducing the reverse feedback", smoke ? " — SMOKE RUN, NOT A RESULT" : "")
    p()
    p("Runs at $commit, merged at $COMMIT. Jobs: ",
      join(["$k $(jobs[k])" for k in sort(collect(keys(jobs)))], ", "), ".")
    p()
    p("Prior SD in ln θ: ", join([@sprintf("%s %.3g", n, s) for (n, s) in zip(SHORT, s0)], ", "), ".")
    p()
    summary = Dict{Any, Any}()
    for t in sort(collect(keys(TRUTHS)))
        p("## Truth $t (`Xoshiro($(TRUTHS[t]))`)")
        p()
        vals = [x.values for a in ARMS for d in (path_dir(out, a, t), ensemble_dir(out, a, t))
                if isdir(d) for x in load_parts(d)]
        all(==(first(vals)), vals) || error("truth $t differs across arms or stages")
        p(join([@sprintf("%s %.4g", k, v) for (k, v) in first(vals)], ", "),
          " — the same in every arm and stage ($(length(vals)) files).")
        p()
        p("### Path information ($NCELLS cells)")
        p()
        p("| Arm | δ | cells | Expected-information ratio: ", join(SHORT, " | "),
          " | ENO, FBA alone | Observed-information ratio: ", join(SHORT, " | "), " |")
        p("|---|---|---|", repeat("---|", 2length(FW) + 1))
        for arm in ARMS
            dir = path_dir(out, arm, t)
            isdir(dir) || (p("| $arm | – | missing | ", repeat(" | ", 2length(FW) + 1)); continue)
            rs = load_parts(dir)
            for δ in DELTAS
                S = reduce(hcat, [[(r.out[(i, δ, 1)] - r.out[(i, δ, -1)]) / 2δ for i in eachindex(FW)]
                                  for r in rs])'
                C = reduce(hcat, [[-(r.out[(i, δ, 1)] - 2r.lp0 + r.out[(i, δ, -1)]) / δ^2
                                   for i in eachindex(FW)] for r in rs])'
                Iexp(rows) = NCELLS .* (S[rows, :]' * S[rows, :]) ./ length(rows)
                rE, rE2 = ratio(Iexp(1:size(S, 1))), ratio(Iexp(1:size(S, 1)), 1:2)
                boot = reduce(hcat, [ratio(Iexp(rand(rng, 1:size(S, 1), size(S, 1)))) for _ in 1:2000])
                lo, hi = quantile.(eachrow(boot), 0.025), quantile.(eachrow(boot), 0.975)
                c̄, se = vec(mean(C; dims = 1)), vec(std(C; dims = 1)) ./ sqrt(size(C, 1))
                rO = [alone(c, i) for (i, c) in enumerate(c̄)]
                rOhi = [alone(c + 2s, i) for (i, (c, s)) in enumerate(zip(c̄, se))]
                summary[(t, arm, δ)] = (; rE, rE2, lo, hi, rO, rOhi, n = length(rs))
                p("| $arm | $δ | $(length(rs)) | ",
                  join([@sprintf("%.3f (%.3f–%.3f)", rE[i], lo[i], hi[i]) for i in eachindex(FW)], " | "), " | ",
                  @sprintf("%.3f, %.3f", rE2...), " | ",
                  join([@sprintf("%.3f (%.3f)", rO[i], rOhi[i]) for i in eachindex(FW)], " | "), " |")
            end
        end
        p()
        p("Expected-information cells are joint over the five parameters, with a 95% bootstrap over cells. ",
          "\"ENO, FBA alone\" is the same from their 2 × 2 block only, the 16a7 script's quantity. Observed-information ",
          "cells are each parameter's own curvature alone, at the mean and, in brackets, at mean + 2 SE (the most ",
          "information the cells allow). NaN is a negative mean curvature; its bracket then assumes the curvature ",
          "is 2 SE more positive than measured, and is not evidence of information.")
        p()

        p("### Bulk metabolite likelihood ($NENSEMBLE-cell means, 60 s saves)")
        p()
        rc, sp, ts = load_ensemble(out, :live60, t)
        Lref = bulk(rc, sp, ts)
        h = length(rc) ÷ 2
        noise = bulk(rc[1:h], sp, ts) .- bulk(rc[(h + 1):(2h)], sp, ts)
        # Two halves differ by twice the variance of one full-size mean's
        # difference from another, so halve the squared half-split difference.
        noise2 = vec(sum(noise .^ 2; dims = 2)) ./ 2
        p("Per pool: Σ_t (log x̄_arm − log x̄_live60)² / (2σ_b²) at σ_b = 10% (it scales as 1/σ_b²), with ",
          "`live60`'s split-half noise on the same scale in the first column. The arms share cell seeds with ",
          "`live60`, so until their paths part the difference is below independent noise, and the noise column ",
          "is conservative. Above the noise is a real shift; above ~1 nat the reduction moves the bulk ",
          "likelihood materially.")
        p()
        arms = [a for a in ARMS if a !== :live60 && isdir(ensemble_dir(out, a, t))]
        p("| Pool | noise | ", join(string.(arms), " | "), " |")
        p("|---|---|", repeat("---|", length(arms)))
        shifts = Dict(a => vec(sum((bulk(load_ensemble(out, a, t)...) .- Lref) .^ 2; dims = 2))
                      for a in arms)
        sb = 0.10
        for (i, o) in enumerate(M0_PANEL)
            p("| `$o` | ", @sprintf("%.3g", noise2[i] / (2sb^2)), " | ",
              join([@sprintf("%.3g", shifts[a][i] / (2sb^2)) for a in arms], " | "), " |")
        end
        p("| **total** | ", @sprintf("**%.3g**", sum(noise2) / (2sb^2)), " | ",
          join([@sprintf("**%.3g**", sum(shifts[a]) / (2sb^2)) for a in arms], " | "), " |")
        p()
        summary[(t, :bulk)] = (; noise2, shifts, sb)
    end
    write(joinpath(out, "merged.md"), take!(io))
    serialize(joinpath(out, "merged.jls"), summary)
    println("Wrote ", joinpath(out, "merged.md"))
end

function main()
    stage = ARGS[1]
    if stage == "ensemble"
        run_ensemble(Symbol(ARGS[2]), parse(Int, ARGS[3]), parse(Int, ARGS[4]), parse(Int, ARGS[5]), ARGS[6])
    elseif stage == "path"
        run_path(Symbol(ARGS[2]), parse(Int, ARGS[3]), parse(Int, ARGS[4]), ARGS[5])
    elseif stage == "merge"
        merge(ARGS[2])
    else
        error("stage must be ensemble, path or merge")
    end
end

main()
