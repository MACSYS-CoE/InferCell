# Merge tasks 14c.5 and 14c.7's runs into their result (spec §11 task 14c.5;
# §12 2026-09-28, planning C and the split gate). Reads dev/scripts/smoothing_14c/,
# which `corea_smoothing_14c.jl` writes, and writes
# dev/scripts/corea_smoothing_14c_result.md.
#
# Every difference is paired per seed and taken relative to max(|published|,
# floor), with a floor of 500 particles for an ODE state (check 1b's) and one
# copy for a transcript. The gate has two halves:
# - coupled: every paired difference at a save point before the pair's first
#   jump-state mismatch must be within 1%;
# - decoupled: at the end of the cycle, the mean over seeds of the paired
#   difference must be within 1% for every observable, with its SE; an SE above
#   1% is unresolved, not passing.
# The sampler pair is read beside the control pair (the published model with
# `krnadeg` scaled by 1 + 1e-5), which shows how often a pair decouples with no
# model change.
#
# Usage: julia --project dev/scripts/corea_smoothing_14c_merge.jl

using Printf
using Statistics: mean, std, median

const DIR = joinpath(@__DIR__, "smoothing_14c")
const OUT = joinpath(@__DIR__, "corea_smoothing_14c_result.md")
const THRESHOLD = 0.01
const ODE_FLOOR = 500.0
const TX_FLOOR = 1.0

files = sort(filter(f -> occursin(r"^agree_task_\d+\.tsv$", f), readdir(DIR)))
isempty(files) && error("no agreement TSVs in $DIR")
headers = String[]
census = NamedTuple[]
couple = Dict{Tuple{Int, Symbol}, Int}()
vals = Dict{Tuple{Int, Symbol, Int, Symbol}, Float64}()
for f in files, line in eachline(joinpath(DIR, f))
    if startswith(line, "# census")
        c = split(line, '\t')
        push!(census, (seed = parse(Int, c[2]), variant = Symbol(c[3]),
                       drains = parse(Int, c[4]), clipped = parse(Int, c[5])))
        continue
    elseif startswith(line, "# couple")
        c = split(line, '\t')
        couple[(parse(Int, c[2]), Symbol(c[3]))] = parse(Int, c[4])
        continue
    end
    startswith(line, "#") && (push!(headers, line[3:end]); continue)
    startswith(line, "seed\t") && continue
    c = split(line, '\t')
    vals[(parse(Int, c[1]), Symbol(c[2]), parse(Int, c[3]), Symbol(c[4]))] = parse(Float64, c[5])
end
smoke = any(h -> occursin("SMOKE", h), headers)

seeds = sort(unique(k[1] for k in keys(vals)))
times = sort(unique(k[3] for k in keys(vals)))
obs = sort(unique(k[4] for k in keys(vals)))
tend = last(times)
floor_of(o) = o === :volume_litres ? 0.0 :
              startswith(String(o), "mRNA_") ? TX_FLOOR : ODE_FLOOR
reldiff(s, v, t, o) = (vals[(s, v, t, o)] - vals[(s, :clamped, t, o)]) /
                      max(abs(vals[(s, :clamped, t, o)]), floor_of(o))
# The last save point strictly before a pair's first jump-state mismatch. A
# mismatch at handshake k means the states after handshake k differ, so the
# save point at k itself is already decoupled.
coupled_times(s, v) = (k = couple[(s, v)]; k < 0 ? times : [t for t in times if t < k])

out = IOBuffer()
p(args...) = (println(out, args...); println(args...))
p("# 14c.5: the sampler model against the published one", smoke ? " — SMOKE RUN, NOT A RESULT" : "")
p()
p("Merged by `dev/scripts/corea_smoothing_14c_merge.jl` from $(length(files)) task files. ",
  "Regenerate with `sbatch dev/scripts/corea_smoothing_14c.slurm`, then this script. The tasks' headers:")
p()
for h in unique(replace.(headers, r", task \d+ of \d+" => ""))
    p("- ", h)
end
p()
p("$(length(seeds)) seeds, each a full cycle under the published model (clamped drain, ",
  "fractional carry), the sampler model (drain smoothed at one particle, pools carried ",
  "continuously), and the control (the published model with `krnadeg` scaled by 1 + 1e-5). ",
  "$(length(obs)) observables (every ODE state in particles with its carried remainder, every ",
  "transcript, the volume) at $(length(times)) save points. Each difference is paired per seed and ",
  "relative to max(|published|, floor): 500 particles for an ODE state, one copy for a transcript, ",
  "none for the volume.")
p()

verdicts = Dict{Symbol, Any}()
for (v, title) in ((:smoothed, "The sampler model"), (:control, "The control"))
    p("## ", title, " against the published model")
    p()
    ks = [couple[(s, v)] for s in seeds]
    dec = [k for k in ks if k >= 0]
    p(@sprintf("**Coupling.** %d of %d pairs decouple within the cycle", length(dec), length(seeds)),
      isempty(dec) ? "." : @sprintf(", at a median handshake of %.0f (earliest %d, latest %d).",
                                    median(dec), minimum(dec), maximum(dec)))
    p()
    # Coupled half.
    worst = (r = 0.0, s = 0, t = 0, o = :none)
    n = 0
    for s in seeds, t in coupled_times(s, v), o in obs
        r = reldiff(s, v, t, o)
        n += 1
        abs(r) > abs(worst.r) && (worst = (r = r, s = s, t = t, o = o))
    end
    cpass = abs(worst.r) <= THRESHOLD
    p(@sprintf("**Coupled half.** %d paired comparisons at save points before decoupling. ", n),
      @sprintf("The largest is %.3g%% (`%s`, seed %d, t = %d s). ", 100worst.r, worst.o, worst.s, worst.t),
      "Gate: every one within 1% — **", cpass ? "passes" : "fails", "**.")
    p()
    # Decoupled half: the ensemble at the end of the cycle.
    st = NamedTuple[]
    for o in obs
        r = [reldiff(s, v, tend, o) for s in seeds]
        push!(st, (o = o, mean = mean(r), se = std(r) / sqrt(length(r)), maxabs = maximum(abs, r)))
    end
    unresolved = [x for x in st if x.se > THRESHOLD]
    failing = [x for x in st if x.se <= THRESHOLD && abs(x.mean) > THRESHOLD]
    p(@sprintf("**Decoupled half, at t = %d s.** Mean over %d seeds of the paired difference, per observable. ",
               tend, length(seeds)),
      "$(length(obs) - length(unresolved) - length(failing)) observables pass, ",
      "$(length(failing)) fail, and $(length(unresolved)) are unresolved (SE above 1%).")
    p()
    p("The mean over seeds, its SE, and z = mean/SE, so a mean far from zero is ",
      "visible even when the SE leaves it unresolved.")
    p()
    p("| observable | mean | SE | z | largest single seed | verdict |")
    p("|---|---|---|---|---|---|")
    for x in sort(st; by = x -> -max(abs(x.mean), x.se))[1:min(12, end)]
        verdict = x.se > THRESHOLD ? "unresolved" : abs(x.mean) > THRESHOLD ? "**fails**" : "passes"
        p(@sprintf("| `%s` | %.3g%% | %.2g%% | %.2f | %.3g%% | %s |", x.o, 100x.mean, 100x.se,
                   x.se > 0 ? x.mean / x.se : 0.0, 100x.maxabs, verdict))
    end
    p()
    # A diagnostic beside the gate, not part of it. The gate's per-seed ratio is
    # asymmetric once a pair decouples: a pool the published run holds at its
    # floor can come out several hundred percent high in the other run, but never
    # below −100%, so decoupled noise averages positive. The difference of the
    # ensemble means, over the published mean floored, has no such skew.
    sym = NamedTuple[]
    for o in obs
        a = [vals[(s, :clamped, tend, o)] for s in seeds]
        d = [vals[(s, v, tend, o)] for s in seeds] .- a
        den = max(abs(mean(a)), floor_of(o))
        den == 0 && (den = 1.0)
        m, se = mean(d) / den, std(d) / sqrt(length(d)) / den
        push!(sym, (o = o, mean = m, se = se, z = se > 0 ? abs(m) / se : 0.0))
    end
    top = sort(sym; by = x -> -x.z)
    p("*Diagnostic, not the gate:* the difference of the ensemble means at t = $tend s, over ",
      "max(|published mean|, floor), which is not skewed by the gate's per-seed ratio. ",
      @sprintf("The largest |z| over %d observables is %.2f (`%s`, %.3g%% ± %.2g%%); ",
               length(obs), top[1].z, top[1].o, 100top[1].mean, 100top[1].se),
      "$(count(x -> x.z > 3, sym)) exceed 3, and $(count(x -> x.se <= THRESHOLD && abs(x.mean) > THRESHOLD, sym)) ",
      "of the $(count(x -> x.se <= THRESHOLD, sym)) with SE within 1% move by more than 1%.")
    p()
    verdicts[v] = (coupled = cpass, failing = length(failing), unresolved = length(unresolved),
                   decoupled = length(dec))
end

s, c = verdicts[:smoothed], verdicts[:control]
p("## Verdict")
p()
p("| | pairs decoupled | coupled half | end-of-cycle observables failing | unresolved |")
p("|---|---|---|---|---|")
for (lab, x) in (("sampler model", s), ("control", c))
    p("| $lab | $(x.decoupled) of $(length(seeds)) | ", x.coupled ? "passes" : "**fails**",
      " | $(x.failing) | $(x.unresolved) |")
end
p()
nobs = length(obs)
p("14c.5's agreement gate: the coupled half ", s.coupled ? "**passes**" : "**fails**",
  ". The end-of-cycle half resolves $(nobs - s.unresolved) of $nobs observables, of which ",
  "$(s.failing) fail", s.unresolved == 0 ? "." :
  ", and it can neither pass nor fail the other $(s.unresolved): their SE exceeds 1%, " *
  "because once a pair decouples the per-seed ratio is Monte Carlo noise, skewed upward. " *
  "For those observables the evidence is the symmetric diagnostic above, which is not the gate.")
p()
p("The control was specified as matched to the sampler model's nudge in size. It is not ",
  "matched in effect: it decouples $(c.decoupled) of $(length(seeds)) pairs against the sampler ",
  "model's $(s.decoupled), so it shows that the model is deterministic enough to stay paired ",
  "under a small nudge, and cannot serve as a null for the sampler model's decoupled half.")
p()
p("## Clipping under each run")
p()
p("The two models' counts use different predicates in effect. A smoothed debit leaves a ",
  "residue `w·exp(−gap/w)`, which the census counts as a clip while the pool is within about ",
  "14 particles above the accrual, although it paid in full to 1e-4 of a particle. So the ",
  "sampler row is not a like-for-like count, and it is recorded, not compared.")
p()
p("| | median drains clipped | seeds carrying a deficit |")
p("|---|---|---|")
for v in (:clamped, :smoothed, :control)
    cs = [x.clipped for x in census if x.variant === v]
    p("| `$v` | $(median(cs)) | $(count(>(0), cs)) of $(length(cs)) |")
end
p()
deriv = joinpath(DIR, "deriv.md")
if isfile(deriv)
    for line in eachline(deriv)
        p(line)
    end
    p()
    p("*Note added at merge.* The clamped jump reads 0.915, not one, because the ",
      "sensitivity it is divided by is taken at 1.5 × the clip, where it is larger. The ",
      "paragraph above, as the driver at the run's commit wrote it, says \"the whole ",
      "sensitivity\"; the driver now says so.")
end
write(OUT, take!(out))
println("Wrote $OUT")
