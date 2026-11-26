# KIM-116 item 7: the modular arm beside the joint, from the two archives alone.
# Reads KIM-92's archive (the joint: pilot_draws.tsv, pilot_coordinates.tsv, jobs.tsv,
# sacct.tsv), KIM-116's (run/cut_draws_<v>.tsv, run/jobs.tsv, predict/jobs.tsv,
# sacct.tsv) and KIM-95's (jobs.tsv, sacct.tsv, for the D3 cost), and writes
# compare.md to <out>. It loads no InferCell code and calls no sacct: every number
# comes from a file in an archive.
#
# The shift and contraction tables are checked against the merge's own
# run/modular_coordinates_<v>.tsv (job 21627295) to 1e-10, so the comparison and the
# merge cannot drift apart.
#
# Usage: julia --project dev/scripts/corea_modular_kim116_compare.jl <kim116-archive> <out>

using DelimitedFiles: readdlm
using Printf
using Statistics: mean, std, quantile, cor

const KIM92 = "/fred/oz022/tkimpson/infercell-kim92-archive-2026-11-17"
const KIM95 = "/fred/oz022/tkimpson/infercell-kim95-archive-2026-11-19"
const VARIANTS = ("cut", "plugin")
const CLOSE_SHIFT, CLOSE_CONTRACTION = 0.5, 0.05       # prereg/kim-116 §5
const CI_LO, CI_HI = 0.05, 0.95
const TOL = 1e-10
# Job groups that make one posterior, per arm. What the arms share (KIM-92's dataset,
# its δ and S reference ensemble, the probe and the profile) is listed apart.
const JOINT_GROUPS = ["metric", "chain", "check", "extend", "ppc", "merge"]
const STAGE1 = ["metric", "chain", "check", "extend"]
const STAGE2 = ["stage2"]
const AFTER = ["ppc"]
const PARALLEL = ("chain", "extend", "stage2", "run_posterior")   # jobs in a group run side by side

# ---------------------------------------------------------------- reading

function table(path)
    A = readdlm(path, '\t', Any; header = false)
    return String.(A[1, :]), A[2:end, :]
end

function draws(path)
    h, X = table(path)
    return Symbol.(h), Float64.(X)
end

# Slurm durations: [D-]HH:MM:SS or MM:SS[.mmm].
function seconds(s)
    s = String(s)
    d = 0
    if occursin('-', s)
        a, s = split(s, '-'; limit = 2)
        d = parse(Int, a)
    end
    p = parse.(Float64, split(s, ':'))
    length(p) == 2 && (p = [0.0; p])
    return d * 86400 + p[1] * 3600 + p[2] * 60 + p[3]
end

# One row per allocation (array tasks included, steps dropped): id, its job (the array's
# id for a task), name, state, start, end, elapsed s, CPUs.
function allocations(path)
    lines = readlines(path)
    h = split(lines[1], '|')
    col(n) = findfirst(==(n), h)
    out = NamedTuple[]
    for l in lines[2:end]
        f = split(l, '|')
        id = f[col("JobID")]
        occursin('.', id) && continue
        push!(out, (; id = String(id), job = String(first(split(id, '_'))), name = String(f[col("JobName")]),
                    state = String(f[col("State")]), start = String(f[col("Start")]), stop = String(f[col("End")]),
                    elapsed = seconds(f[col("Elapsed")]), cpus = parse(Int, f[col("AllocCPUS")])))
    end
    return out
end

function groups(paths...)
    g = Dict{String, String}()
    for p in paths, l in readlines(p)
        name, id = split(l, '\t')
        g[String(id)] = String(name)
    end
    return g
end

# CPU-h and critical-path hours by group. Within a PARALLEL group the jobs overlap, so
# the path takes the longest; otherwise the jobs ran one after another and it takes
# their sum.
function cost(rows, g; label = "")
    cpuh, path = Dict{String, Float64}(), Dict{String, Float64}()
    for r in rows
        name = get(g, r.job, nothing)
        if name === nothing
            println(@sprintf("%s: allocation %s (%s, %s) is not in jobs.tsv; left out, %.2f CPU-h", label, r.id,
                             r.name, r.state, r.elapsed * r.cpus / 3600))
            continue
        end
        r.state == "COMPLETED" || error("job $(r.id) ($name) is $(r.state)")
        cpuh[name] = get(cpuh, name, 0.0) + r.elapsed * r.cpus / 3600
    end
    for name in keys(cpuh)
        hs = [r.elapsed / 3600 for r in rows if get(g, r.job, "") == name]
        path[name] = any(p -> startswith(name, p), PARALLEL) ? maximum(hs) : sum(hs)
    end
    return cpuh, path
end

total(d, ks) = sum(get(d, k, 0.0) for k in ks; init = 0.0)

# The Q rows and the family rows of a KIM-95 score.md, by their first cell.
function score_rows(path)
    rows = Dict{String, Vector{String}}()
    for l in readlines(path)
        startswith(l, "| ") || continue
        f = String.(strip.(split(l, '|')[2:end-1]))
        f[1] in ("Q", "Q_kd", "dQ", "pool", "transcript") && (rows[f[1]] = f)
    end
    length(rows) == 5 || error("$path: expected the Q, Q_kd, dQ, pool and transcript rows")
    return rows
end

# Stage 2's per-chain summary lines in run/logs/kim116_stage2_<v>_*.log.
function stage2_log(A, v)
    acc, err, depth = Float64[], 0, 0
    for f in readdir(joinpath(A, "run", "logs"); join = true)
        occursin("kim116_stage2_$(v)_", basename(f)) || continue
        for l in readlines(f)
            m = match(r"stage 2 chain \d+ of \d+: .*mean acceptance ([0-9.]+), (\d+) numerical errors, max depth (\d+)", l)
            m === nothing && continue
            push!(acc, parse(Float64, m[1]))
            err += parse(Int, m[2])
            depth = max(depth, parse(Int, m[3]))
        end
    end
    return (; n = length(acc), acc = mean(acc), minacc = minimum(acc), err, depth)
end

# ---------------------------------------------------------------- statistics

function summarise(X, psd)
    μ, sd = vec(mean(X; dims = 1)), vec(std(X; dims = 1))
    return (; μ, sd, contraction = 1 .- sd .^ 2 ./ psd .^ 2,
            lo = [quantile(X[:, q], CI_LO) for q in axes(X, 2)],
            hi = [quantile(X[:, q], CI_HI) for q in axes(X, 2)], R = cor(X))
end

function check_against_merge(path, names, J, M, shift)
    h, T = table(path)
    c(n) = findfirst(==(n), h)
    Symbol.(T[:, c("coordinate")]) == names || error("$path: coordinates differ")
    for (n, v) in (("joint_mean", J.μ), ("joint_sd", J.sd), ("cut_mean", M.μ), ("cut_sd", M.sd),
                   ("shift", shift), ("joint_contraction", J.contraction), ("cut_contraction", M.contraction))
        e = maximum(abs.(Float64.(T[:, c(n)]) .- v))
        e <= TOL * max(1, maximum(abs.(v))) || error("$path: $n differs from the merge by $e")
    end
end

# ---------------------------------------------------------------- main

function main(args)
    length(args) == 2 || error("usage: corea_modular_kim116_compare.jl <kim116-archive> <out>")
    A, out = args
    mkpath(out)
    names, J0 = draws(joinpath(KIM92, "pilot_draws.tsv"))
    h, C92 = table(joinpath(KIM92, "pilot_coordinates.tsv"))
    Symbol.(C92[:, 1]) == names || error("KIM-92's coordinates are not its draws' columns")
    psd = Float64.(C92[:, findfirst(==("prior_sd"), h)])
    truth = Float64.(C92[:, findfirst(==("truth"), h)])
    ess = Float64.(C92[:, findfirst(==("ess"), h)])        # the joint's bulk ESS, gate_400
    P = length(names)
    eno, fba = findfirst(==(:log_kcatF_R_ENO), names), findfirst(==(:log_kcatF_R_FBA), names)
    phi = setdiff(1:P, [eno, fba, findfirst(==(:log_sigma_b), names)])
    J = summarise(J0, psd)

    io = IOBuffer()
    p(xs...) = (println(io, xs...); println(xs...))
    p("# KIM-116: the modular arm beside the joint\n")
    p("From the archives alone: ", KIM92, " (the joint, ", size(J0, 1), " draws), ", A,
      " (the modular arm) and ", KIM95, " (D3). Written by `dev/scripts/corea_modular_kim116_compare.jl`. ",
      "\"Close\" (prereg/kim-116 §5): every |shift| ≤ ", CLOSE_SHIFT, " joint SD and every |Δ contraction| ≤ ",
      CLOSE_CONTRACTION, ". The shift and contraction columns match the merge's `modular_coordinates_<v>.tsv` ",
      "(job 21627295) to ", TOL, ".\n")

    # ------------------------------------------------ cost
    r92 = allocations(joinpath(KIM92, "sacct.tsv"))
    g92 = groups(joinpath(KIM92, "jobs.tsv"))
    cpu92, path92 = cost(r92, g92; label = "KIM-92")
    # KIM-92's first pass (precompile, dataset, the profile before KIM-99), in its own
    # sacct.tsv: these make T1's 491.9 with the merge's 486.60.
    first92 = sum(r.elapsed * r.cpus / 3600 for r in allocations(joinpath(KIM92, "pilot-2026-11-10", "sacct.tsv")))
    r116 = allocations(joinpath(A, "sacct.tsv"))
    g116 = groups(joinpath(A, "run", "jobs.tsv"), joinpath(A, "predict", "jobs.tsv"))
    cpu116, path116 = cost(r116, g116; label = "KIM-116")
    r95 = allocations(joinpath(KIM95, "sacct.tsv"))
    cpu95, path95 = cost(r95, groups(joinpath(KIM95, "jobs.tsv")); label = "KIM-95")

    p("## Cost\n")
    p("CPU-h is Elapsed × AllocCPUS from each archive's `sacct.tsv`, summed over allocations (array tasks included). ",
      "The critical path adds the longest job of each group that ran side by side (chains, extensions, stage-2 tasks) ",
      "and the sum of the rest; it leaves out time in the queue.\n")
    p("| arm | stage 1 CPU-h | stage 2 CPU-h | PPC + merge CPU-h | total CPU-h | ÷ joint | critical path, h |")
    p("|---|---|---|---|---|---|---|")
    jt = total(cpu92, JOINT_GROUPS)
    jpath = total(path92, JOINT_GROUPS)
    p(@sprintf("| joint (KIM-92) | %.1f | — | %.1f | %.1f | 1.00 | %.1f |", total(cpu92, STAGE1),
               total(cpu92, ["ppc", "merge"]), jt, jpath))
    rows = Dict{String, Any}()
    for v in VARIANTS
        k(gs) = [g * "_" * v for g in gs]
        s1, s2, pm = total(cpu116, k(STAGE1)), total(cpu116, k(STAGE2)), total(cpu116, k(AFTER)) + get(cpu116, "merge", 0.0)
        t = s1 + s2 + pm
        cp = total(path116, k([STAGE1; STAGE2; AFTER])) + get(path116, "merge", 0.0)
        rows[v] = (; s1, s2, pm, t, cp)
        p(@sprintf("| %s | %.1f | %.1f | %.1f | %.1f | %.2f | %.1f |", v, s1, s2, pm, t, t / jt, cp))
    end
    shared = setdiff(unique(values(g92)), JOINT_GROUPS)
    sh = total(cpu92, shared)
    p(@sprintf("\nShared by every arm, paid once in KIM-92: %s, %.1f CPU-h. The merge (%.2f CPU-h) serves both variants and is counted in each.",
               join(sort(shared), ", "), sh, get(cpu116, "merge", 0.0)))
    p(@sprintf("The joint's %.1f here, plus the shared %.1f, is KIM-92's merge total (%.2f); KIM-92's first pass (`pilot-2026-11-10/sacct.tsv`, %.2f) brings it to T1's %.1f.",
               jt, sh, jt + sh, first92, jt + sh + first92))
    p(@sprintf("Against T1's %.1f (prereg §5's comparator), the totals are cut %.2f× and plug-in %.2f×; the table's ÷ joint uses %.1f, the joint's own posterior jobs, so that the shared jobs count in neither arm.",
               jt + sh + first92, rows["cut"].t / (jt + sh + first92), rows["plugin"].t / (jt + sh + first92), jt))
    p(@sprintf("Stage 1 against the joint's chains, check and extension: cut %.2f×, plug-in %.2f×. Stage 2's share of each variant: cut %.0f%%, plug-in %.0f%%; its size is the frozen K = 200 × NUTS(100, 0.8), not tuned for cost.",
               rows["cut"].s1 / total(cpu92, STAGE1), rows["plugin"].s1 / total(cpu92, STAGE1),
               100 * rows["cut"].s2 / rows["cut"].t, 100 * rows["plugin"].s2 / rows["plugin"].t))
    p(@sprintf("\nD3 scores (the posterior arm only, 100 draws × 973 cells × 2 arms): plug-in %.1f CPU-h, critical path %.1f h; joint (KIM-95) %.1f CPU-h, %.1f h. ",
               get(cpu116, "run_posterior", 0.0), get(path116, "run_posterior", 0.0),
               get(cpu95, "run_posterior", 0.0), get(path95, "run_posterior", 0.0)),
      "(`collect.md` counts the chunks' timed seconds × threads instead, and gives a little less.) ",
      "Wall-clock here is the critical path only; time in the queue is not counted.\n")

    # ------------------------------------------------ posterior
    for v in VARIANTS
        nv, X = draws(joinpath(A, "run", "cut_draws_$v.tsv"))
        nv == names || error("$v's cut draws are not in the joint's coordinate order")
        M = summarise(X, psd)
        shift = (M.μ .- J.μ) ./ J.sd
        dcon = M.contraction .- J.contraction
        check_against_merge(joinpath(A, "run", "modular_coordinates_$v.tsv"), names, J, M, shift)
        far = [names[q] for q in 1:P if abs(shift[q]) > CLOSE_SHIFT]
        wide = [names[q] for q in 1:P if abs(dcon[q]) > CLOSE_CONTRACTION]
        inside = M.lo .<= truth .<= M.hi
        jin = J.lo .<= truth .<= J.hi
        # φ is constant within each stage-1 draw's block of inner draws, so the cut's φ
        # rests on its distinct stage-1 draws, not on all its rows. Monte Carlo SE of a
        # shift: the cut's mean on that many draws, the joint's on its ESS, in joint SDs.
        nphi = size(unique(X[:, phi]; dims = 1), 1)
        se = sqrt.(M.sd .^ 2 ./ nphi .+ J.sd .^ 2 ./ ess) ./ J.sd
        st2 = stage2_log(A, v)
        p("## Variant: ", v, "\n")
        p(@sprintf("%d draws on %d distinct stage-1 draws. Close to the joint: **%s**.", size(X, 1), nphi,
                   isempty(far) && isempty(wide) ? "yes" : "no"),
          isempty(far) ? "" : " |shift| > $(CLOSE_SHIFT): " * join(string.(far), ", ") * ".",
          isempty(wide) ? "" : " |Δ contraction| > $(CLOSE_CONTRACTION): " * join(string.(wide), ", ") * ".")
        p(@sprintf("Stage 2 (prereg §4, from `run/logs`): %d inner chains, mean acceptance %.2f (lowest chain %.2f), %d numerical errors, max tree depth %d.",
                   st2.n, st2.acc, st2.minacc, st2.err, st2.depth))
        p("Not pre-registered, descriptive only: the Monte Carlo SE, the mean shift and its count, the sd ratio, and the correlation scale below.")
        qm = argmax(abs.(shift))
        p(@sprintf("Shift: max |shift| %.2f joint SD (%s, Monte Carlo SE %.2f); mean shift over the 17 promoters %+.2f; %d of 17 promoters shift up. ",
                   abs(shift[qm]), names[qm], se[qm], mean(shift[1:17]), count(>(0), shift[1:17])),
          @sprintf("Contraction: max |Δ| %.3f. Width (this variant's sd ÷ joint sd): median %.2f, range %.2f–%.2f. ",
                   maximum(abs.(dcon)), quantile(M.sd ./ J.sd, 0.5), minimum(M.sd ./ J.sd), maximum(M.sd ./ J.sd)),
          @sprintf("Truth in the 90%% interval: this variant %d of %d, joint %d of %d.", count(inside), P, count(jin), P))
        p(@sprintf("Promoters' (mean − truth)/sd: this variant's mean %+.2f, %d of 17 negative; joint mean %+.2f, %d of 17 negative.\n",
                   mean(((M.μ .- truth) ./ M.sd)[1:17]), count(<(0), ((M.μ .- truth) ./ M.sd)[1:17]),
                   mean(((J.μ .- truth) ./ J.sd)[1:17]), count(<(0), ((J.μ .- truth) ./ J.sd)[1:17])))
        p("| coordinate | shift (joint SD) | MC SE | Δ contraction | sd ratio | truth in this variant's 90% | truth in joint 90% |")
        p("|---|---|---|---|---|---|---|")
        for q in 1:P
            p(@sprintf("| %s | %+.2f | %.2f | %+.4f | %.2f | %s | %s |", names[q], shift[q], se[q], dcon[q],
                       M.sd[q] / J.sd[q], inside[q] ? "yes" : "no", jin[q] ? "yes" : "no"))
        end
        rj = [J.R[phi, eno]; J.R[phi, fba]]
        rc = [M.R[phi, eno]; M.R[phi, fba]]
        # The cut keeps φ→ψ dependence (stage 2 is fitted at each stage-1 draw) and drops
        # only the bulk's feedback into φ. 1/√n is a scale for sampling noise, not a test.
        nj = quantile(ess, 0.5)
        p(@sprintf("\nENO/FBA–φ correlations (36 pairs): joint max |r| %.3f, rms %.3f; this variant max |r| %.3f, rms %.3f. For scale, 1/√n is %.3f on its %d distinct stage-1 draws and %.3f on the joint's median ESS %.0f.\n",
                   maximum(abs.(rj)), sqrt(mean(rj .^ 2)), maximum(abs.(rc)), sqrt(mean(rc .^ 2)),
                   1 / sqrt(nphi), nphi, 1 / sqrt(nj), nj))
    end
    # ------------------------------------------------ D3
    d116, d95 = score_rows(joinpath(A, "predict", "score.md")), score_rows(joinpath(KIM95, "score.md"))
    p("## D3 scores: the plug-in beside the joint\n")
    p("Read from KIM-95's scorer's `score.md` in each archive (", joinpath(A, "predict"), ", job 21408654; ",
      KIM95, "). The truth and the prior arm are the same files (SHA256-checked in the plug-in's collect). ",
      "Only the plug-in was scored (prereg/kim-116 §6); the cut, the variant close to the joint, has no D3 scores.\n")
    p("| target | truth ± SE | plug-in 5–50–95% | joint 5–50–95% | in 50/90: plug-in | joint | CRPS plug-in | joint | skill plug-in | joint |")
    p("|---|---|---|---|---|---|---|---|---|---|")
    for t in ("Q", "Q_kd", "dQ")
        a, b = d116[t], d95[t]
        a[2] == b[2] || error("$t: the truth differs between the two score.md")
        p("| ", join((t, a[2], a[3], b[3], a[5], b[5], a[7], b[7], a[9], b[9]), " | "), " |")
    end
    p("\n| family | hit 50: plug-in | joint | hit 90: plug-in | joint | median species skill: plug-in | joint |")
    p("|---|---|---|---|---|---|---|")
    for t in ("pool", "transcript")
        a, b = d116[t], d95[t]
        p("| ", join((t, a[4], b[4], a[5], b[5], a[8], b[8]), " | "), " |")
    end
    p()
    write(joinpath(out, "compare.md"), String(take!(io)))
end

main(ARGS)
