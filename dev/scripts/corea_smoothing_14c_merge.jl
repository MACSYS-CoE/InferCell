# Merge task 14c.5's runs into its result (spec §11 task 14c.5, §12 2026-09-28
# planning C). Reads dev/scripts/smoothing_14c/, which `corea_smoothing_14c.jl`
# writes, and writes dev/scripts/corea_smoothing_14c_result.md.
#
# The gate, fixed at planning: for every candidate observable at every 60 s save
# point, the sampler-minus-published difference is taken per seed relative to
# max(|published|, floor), with a floor of 500 particles for an ODE state (check
# 1b's) and one copy for a transcript. The mean over seeds must be within 1%
# (D10's threshold), and is reported with its standard error.
#
# Usage: julia --project dev/scripts/corea_smoothing_14c_merge.jl

using Printf
using Statistics: mean, std

const DIR = joinpath(@__DIR__, "smoothing_14c")
const OUT = joinpath(@__DIR__, "corea_smoothing_14c_result.md")
const THRESHOLD = 0.01
const ODE_FLOOR = 500.0
const TX_FLOOR = 1.0

files = sort(filter(f -> occursin(r"^agree_task_\d+\.tsv$", f), readdir(DIR)))
isempty(files) && error("no agreement TSVs in $DIR")
headers = String[]
census = NamedTuple[]
vals = Dict{Tuple{Int, Symbol, Int, Symbol}, Float64}()
for f in files, line in eachline(joinpath(DIR, f))
    if startswith(line, "# census")
        c = split(line, '\t')
        push!(census, (seed = parse(Int, c[2]), drain = Symbol(c[3]),
                       drains = parse(Int, c[4]), clipped = parse(Int, c[5]),
                       wall = parse(Float64, c[6])))
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
floor_of(o) = o === :volume_litres ? 0.0 :
              startswith(String(o), "mRNA_") ? TX_FLOOR : ODE_FLOOR

stats = NamedTuple[]
for o in obs, t in times
    r = Float64[]
    for s in seeds
        a = get(vals, (s, :clamped, t, o), NaN)
        b = get(vals, (s, :smoothed, t, o), NaN)
        push!(r, (b - a) / max(abs(a), floor_of(o)))
    end
    any(isnan, r) && error("missing a run for $o at t = $t")
    m = mean(r)
    se = length(r) > 1 ? std(r) / sqrt(length(r)) : NaN
    push!(stats, (obs = o, t = t, mean = m, se = se, maxabs = maximum(abs, r)))
end

out = IOBuffer()
p(args...) = (println(out, args...); println(args...))
p("# 14c.5: the smoothed model against the clipped one", smoke ? " — SMOKE RUN, NOT A RESULT" : "")
p()
p("Merged by `dev/scripts/corea_smoothing_14c_merge.jl` from $(length(files)) task files. ",
  "Regenerate with `sbatch dev/scripts/corea_smoothing_14c.slurm`, then this script. The tasks' headers:")
p()
for h in unique(replace.(headers, r", task \d+ of \d+" => ""))
    p("- ", h)
end
p()
p("## Agreement on the candidate observables")
p()
p("$(length(seeds)) seeds ($(join(seeds, ", "))), each a full cycle under the published ",
  "model (clamped drain, fractional carry) and under the sampler model (drain smoothed at ",
  "one particle, pools carried continuously), paired per seed. ",
  "$(length(obs)) observables (every ODE state in particles with its carried remainder, ",
  "every transcript, the volume) ",
  "at $(length(times)) save points, so $(length(stats)) comparisons. Each difference is ",
  "relative to max(|published|, floor): 500 particles for an ODE state, one copy for a ",
  "transcript, none for the volume. The gate is the mean over seeds within 1%.")
p()
worst = sort(stats; by = s -> -abs(s.mean))
final = [s for s in stats if s.t == last(times)]
wfinal = sort(final; by = s -> -abs(s.mean))
fails = count(s -> abs(s.mean) > THRESHOLD, stats)
p(@sprintf("**Largest |mean| over all %d comparisons: %.3g%%** (`%s` at t = %d s, SE %.2g%%). ",
           length(stats), 100worst[1].mean, worst[1].obs, worst[1].t, 100worst[1].se),
  @sprintf("At the pre-chosen instant, the end of the cycle (t = %d s), the largest is %.3g%% (`%s`, SE %.2g%%). ",
           last(times), 100wfinal[1].mean, wfinal[1].obs, 100wfinal[1].se),
  "Comparisons outside 1%: **$fails of $(length(stats))**.")
p()
p("The ten largest, by |mean|:")
p()
p("| observable | t (s) | mean relative difference | SE | largest single seed |")
p("|---|---|---|---|---|")
for s in worst[1:min(10, end)]
    p(@sprintf("| `%s` | %d | %.3g%% | %.2g%% | %.3g%% |", s.obs, s.t, 100s.mean, 100s.se, 100s.maxabs))
end
p()
p(@sprintf("Verdict against the 1%% threshold: **%s**.", fails == 0 ? "passes" : "fails"))
p()
p("## Clipping under each model")
p()
p("| seed | published: drains clipped | sampler: drains clipped |")
p("|---|---|---|")
for s in seeds
    a = only(c for c in census if c.seed == s && c.drain === :clamped)
    b = only(c for c in census if c.seed == s && c.drain === :smoothed)
    p("| $s | $(a.clipped) of $(a.drains) | $(b.clipped) of $(b.drains) |")
end
p()
deriv = joinpath(DIR, "deriv.md")
if isfile(deriv)
    for line in eachline(deriv)
        p(line)
    end
end
write(OUT, take!(out))
println("Wrote $OUT")
