# Merge the census array's TSVs into check 7's result and K5's score (spec §11
# task 14.8, §8 K5). Reads dev/scripts/census/task_*.tsv, which
# `corea_census.jl` writes, and writes dev/scripts/corea_census_result.md.
#
# K5 fires when any deficit is carried at published parameters, or when more
# than 5% of prior draws clip. It is recorded whatever it shows: check 7 is
# scored, not gated (§12, 2026-09-25 C).
#
# Usage: julia --project dev/scripts/corea_census_merge.jl

using Printf
using Statistics: median

const DIR = joinpath(@__DIR__, "census")
const OUT = joinpath(@__DIR__, "corea_census_result.md")
const K5_DRAW_FRACTION = 0.05
const DRAW_SEED0 = 14_800

files = sort(filter(f -> occursin(r"^task_\d+\.tsv$", f), readdir(DIR)))
isempty(files) && error("no census TSVs in $DIR")
headers = String[]
rows = NamedTuple[]
for f in files
    for line in eachline(joinpath(DIR, f))
        startswith(line, "#") && (push!(headers, line[3:end]); continue)
        startswith(line, "kind\t") && continue
        c = split(line, '\t'; keepempty = true)
        clips = isempty(c[9]) ? NamedTuple[] :
            [(counter = x[1], drains = parse(Int, x[2]), first = parse(Float64, x[3]),
              max_deficit = parse(Float64, x[4])) for x in split.(split(c[9], ';'), ':')]
        push!(rows, (kind = Symbol(c[1]), id = parse(Int, c[2]), seed = parse(Int, c[3]),
                     status = c[4], drains = c[5] == "" ? 0 : parse(Int, c[5]),
                     clipped = c[6] == "" ? 0 : parse(Int, c[6]),
                     fraction = c[7] == "" ? NaN : parse(Float64, c[7]),
                     wall = parse(Float64, c[8]), clips = clips, message = c[10]))
    end
end

pub = sort([r for r in rows if r.kind === :published]; by = r -> r.seed)
draws = sort([r for r in rows if r.kind === :draw]; by = r -> r.id)
nominal = [r for r in draws if r.id == 0]
prior = [r for r in draws if r.id > 0]
ok = [r for r in prior if r.status == "ok"]
failed = [r for r in prior if r.status != "ok"]
clipping = [r for r in ok if r.clipped > 0]

out = IOBuffer()
p(args...) = (println(out, args...); println(args...))
smoke = any(h -> occursin("SMOKE", h), headers)

p("# Check 7: the clipping census, scored for K5", smoke ? " — SMOKE RUN, NOT A RESULT" : "")
p()
p("Merged by `dev/scripts/corea_census_merge.jl` from $(length(files)) task files. ",
  "Regenerate with `sbatch dev/scripts/corea_census.slurm`, then this script. The tasks' headers:")
p()
for h in unique(replace.(headers, r", task \d+ of \d+" => ""))
    p("- ", h)
end
p()
p("A drain clips when a consumer counter's pool pays less than was accrued, so a deficit is carried. ",
  "The fraction is per drain, at the published 1 s drain.")
p()

p("## At published parameters")
p()
p("| seed | status | drains | drains clipped | fraction | counters that clipped (drains, first s) |")
p("|---|---|---|---|---|---|")
for r in pub
    p(@sprintf("| %d | %s | %d | %d | %.4g | %s |", r.seed, r.status, r.drains, r.clipped,
               r.fraction, isempty(r.clips) ? "none" :
               join((@sprintf("`%s` (%d, %.0f)", c.counter, c.drains, c.first) for c in r.clips), ", ")))
end
p()
npub = count(r -> r.clipped > 0, pub)
p("Seeds carrying a deficit: **$npub of $(length(pub))**. K5's threshold here is zero.")
p()

p("## Draw 0: freeing the constants changes nothing")
p()
same = r -> [r2 for r2 in pub if r2.seed == r.seed]
for r in nominal
    q = same(r)
    if isempty(q)
        p("Draw 0 (seed $(r.seed)) has no published run at its seed to compare with.")
    else
        eq = q[1].drains == r.drains && q[1].clipped == r.clipped && q[1].clips == r.clips
        p("Draw 0 at seed $(r.seed): $(r.clipped) drains clipped, against $(q[1].clipped) for the ",
          "published build at the same seed. Per-counter records ", eq ? "**identical**" : "**DIFFER**", ".")
    end
end
p()

p("## Across prior draws")
p()
n = length(prior)
p("$n draws; $(length(ok)) completed, $(length(failed)) failed.")
if !isempty(ok)
    f = length(clipping) / length(ok)
    # Wilson 95% interval on the fraction of draws that clip.
    z = 1.96; m = length(ok)
    c = (f + z^2 / 2m) / (1 + z^2 / m)
    h = z * sqrt(f * (1 - f) / m + z^2 / 4m^2) / (1 + z^2 / m)
    p(@sprintf("Draws carrying any deficit: **%d of %d, %.1f%%** (Wilson 95%% interval %.1f–%.1f%%). K5's threshold is 5%%.",
               length(clipping), m, 100f, 100(c - h), 100(c + h)))
    fr = [r.fraction for r in ok]
    p(@sprintf("Per-drain clipping fraction over completed draws: median %.3g, max %.3g.", median(fr), maximum(fr)))
    p()
    counters = sort(unique(c.counter for r in ok for c in r.clips))
    if !isempty(counters)
        p("| counter | draws in which it clipped | median drains clipped when it did | earliest first clip (s) |")
        p("|---|---|---|---|")
        for k in counters
            ds = [only(c for c in r.clips if c.counter == k) for r in ok if any(c -> c.counter == k, r.clips)]
            p(@sprintf("| `%s` | %d | %.0f | %.0f |", k, length(ds), median([c.drains for c in ds]),
                       minimum(c.first for c in ds)))
        end
        p()
    end
end
if !isempty(failed)
    p("Failed draws, with the error each raised:")
    p()
    for r in failed
        p("- draw $(r.id) (seed $(r.seed)): ", r.message)
    end
    p()
end

p("## K5 on the scoreboard (T3)")
p()
pub_fires = npub > 0
draw_frac = isempty(ok) ? NaN : length(clipping) / length(ok)
draw_fires = draw_frac > K5_DRAW_FRACTION
p("| criterion | threshold | measured | verdict |")
p("|---|---|---|---|")
p(@sprintf("| K5, published parameters | 0 seeds carrying a deficit | %d of %d | %s |",
           npub, length(pub), pub_fires ? "**fires**" : "passes"))
p(@sprintf("| K5, prior draws | ≤ 5%% of draws clip | %.1f%% of %d | %s |",
           100draw_frac, length(ok), draw_fires ? "**fires**" : "passes"))
p()
p("K5 ", pub_fires || draw_fires ? "**fires**" : "does not fire", ". What follows from that — smoothing ",
  "the drain or resizing a pool — is a separate decision (spec §12, 2026-09-25 C).")

write(OUT, take!(out))
println("Wrote $OUT")
