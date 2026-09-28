# Merge phase 14c's three census reruns into one result (spec §11 tasks 14c.2,
# 14c.3 and 14c.6). Reads dev/scripts/census_14c/<mode>/task_*.tsv, which
# `corea_census_14c.jl` writes, and 14b's dev/scripts/census/task_*.tsv, which
# it sets beside them. Writes dev/scripts/corea_census_14c_result.md.
#
# Only 14c.2 is scored: K5's prior half, rescored over D11's targets, fires when
# more than 5% of draws clip. 14c.3 and 14c.6 are recorded, not gated.
#
# A draw is *starved* when some counter clips at more than half its drains
# (§11 task 14c.3): its pool sits at zero for most of the cycle, which is past
# the kink rather than near it.
#
# 14b's task files are gitignored. In a fresh clone, run 14b's census first
# (`sbatch dev/scripts/corea_census.slurm`), which writes dev/scripts/census/.
#
# Usage: julia --project dev/scripts/corea_census_14c_merge.jl

using Printf
using Statistics: median

const DIR = joinpath(@__DIR__, "census_14c")
const DIR_14B = joinpath(@__DIR__, "census")
const OUT = joinpath(@__DIR__, "corea_census_14c_result.md")
const K5_DRAW_FRACTION = 0.05

function read_rows(dir)
    isdir(dir) || error("$dir does not exist. 14b's census is gitignored: run " *
                        "`sbatch dev/scripts/corea_census.slurm` first")
    files = sort(filter(f -> occursin(r"^task_\d+\.tsv$", f), readdir(dir)))
    isempty(files) && error("no census TSVs in $dir")
    headers = String[]
    rows = NamedTuple[]
    for f in files, line in eachline(joinpath(dir, f))
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
    return headers, rows
end

starved(r) = any(c -> c.drains > r.drains / 2, r.clips)

function wilson(k, m; z = 1.96)
    f = k / m
    c = (f + z^2 / 2m) / (1 + z^2 / m)
    h = z * sqrt(f * (1 - f) / m + z^2 / 4m^2) / (1 + z^2 / m)
    return c - h, c + h
end

# The summary of a set of prior draws, as 14b's result reports it plus starved.
function draw_summary(rows)
    prior = sort([r for r in rows if r.kind === :draw && r.id > 0]; by = r -> r.id)
    ok = [r for r in prior if r.status == "ok"]
    m = length(ok)
    k = count(r -> r.clipped > 0, ok)
    s = count(starved, ok)
    lo, hi = m == 0 ? (NaN, NaN) : wilson(k, m)
    slo, shi = m == 0 ? (NaN, NaN) : wilson(s, m)
    return (n = length(prior), ok = ok, failed = [r for r in prior if r.status != "ok"],
            clipping = k, frac = k / m, lo = lo, hi = hi, starved = s, sfrac = s / m,
            slo = slo, shi = shi, median_fraction = median([r.fraction for r in ok]),
            median_starved = median([r.fraction for r in ok if starved(r)]),
            median_other = median([r.fraction for r in ok if !starved(r)]),
            max_fraction = maximum(r.fraction for r in ok))
end

out = IOBuffer()
p(args...) = (println(out, args...); println(args...))

modes = [m for m in ("d11", "broad", "gtp") if isdir(joinpath(DIR, m))]
data = Dict(m => read_rows(joinpath(DIR, m)) for m in modes)
headers_14b, rows_14b = read_rows(DIR_14B)
smoke = any(h -> occursin("SMOKE", h), vcat((first(data[m]) for m in modes)...))

p("# Phase 14c: check 7's census, rerun three ways", smoke ? " — SMOKE RUN, NOT A RESULT" : "")
p()
p("Merged by `dev/scripts/corea_census_14c_merge.jl`. Regenerate with ",
  "`sbatch dev/scripts/corea_census_14c.slurm <mode>` for each of `d11`, `broad` and `gtp`, then ",
  "this script. 14b's task files, which this merge also reads, are gitignored: in a fresh ",
  "clone run `sbatch dev/scripts/corea_census.slurm` first. Every run is at the published ",
  "clamped drain. The tasks' headers:")
p()
for m in modes, h in unique(replace.(first(data[m]), r", task \d+ of \d+" => ""))
    p("- ", h)
end
p()
p("14b's census, set beside these, is `dev/scripts/corea_census_result.md` at ",
  "`7e3c1c9`, recomputed here from its task files. A drain clips when a consumer ",
  "counter's pool pays less than was accrued. A draw is starved when some counter ",
  "clips at more than half its drains.")
p()

s14b = draw_summary(rows_14b)

function draw_section(mode, title)
    rows = last(data[mode])
    s = draw_summary(rows)
    ref = [r for r in rows if r.kind === :published]
    d0 = [r for r in rows if r.kind === :draw && r.id == 0]
    for r in ref, z in d0
        eq = r.drains == z.drains && r.clipped == z.clipped && r.clips == z.clips
        p("Draw 0 (the freed build, no value written) at seed ",
          "$(z.seed): $(z.clipped) drains clipped, against $(r.clipped) for the published ",
          "build at the same seed. Per-counter records ", eq ? "**identical**" : "**DIFFER**", ".")
        p()
    end
    p("$(s.n) draws; $(length(s.ok)) completed, $(length(s.failed)) failed.")
    p()
    p("| | draws clipping | Wilson 95% | starved draws | Wilson 95% | median per-drain fraction | starved draws' median | the others' median | max |")
    p("|---|---|---|---|---|---|---|---|---|")
    for (lab, x) in ((title, s), ("14b broad, independent draws", s14b))
        p(@sprintf("| %s | %d of %d, %.1f%% | %.1f–%.1f%% | %d of %d, %.1f%% | %.1f–%.1f%% | %.3g | %.3g | %.3g | %.3g |",
                   lab, x.clipping, length(x.ok), 100x.frac, 100x.lo, 100x.hi,
                   x.starved, length(x.ok), 100x.sfrac, 100x.slo, 100x.shi,
                   x.median_fraction, x.median_starved, x.median_other, x.max_fraction))
    end
    p()
    p("The per-drain fraction is bimodal: a starved draw clips on most drains and the ",
      "others on almost none. The median over all draws therefore falls in whichever mode ",
      "holds more than half of them, and moves with the starved count rather than with how ",
      "long a starved draw clips. Read the split columns.")
    p()
    counters = sort(unique(c.counter for r in s.ok for c in r.clips))
    if !isempty(counters)
        p("| counter | draws in which it clipped | median drains clipped when it did | earliest first clip (s) |")
        p("|---|---|---|---|")
        for k in counters
            ds = [only(c for c in r.clips if c.counter == k) for r in s.ok
                  if any(c -> c.counter == k, r.clips)]
            p(@sprintf("| `%s` | %d | %.0f | %.0f |", k, length(ds),
                       median([c.drains for c in ds]), minimum(c.first for c in ds)))
        end
        p()
    end
    for r in s.failed
        p("- draw $(r.id) (seed $(r.seed)) failed: ", r.message)
    end
    isempty(s.failed) || p()
    return s
end

if "d11" in modes
    p("## 14c.2: K5's prior half, over D11's six targets")
    p()
    p("Each draw varies the three promoter strengths, `krnadeg`, and ENO's and FBA's ",
      "forward constants from their priors, with the two reverse constants derived ",
      "(task 14c.1). Everything else is at its published value.")
    p()
    s = draw_section("d11", "14c.2, D11's six, Haldane-consistent")
    fires = s.frac > K5_DRAW_FRACTION
    p("| criterion | threshold | measured | verdict |")
    p("|---|---|---|---|")
    p(@sprintf("| K5, prior draws over D11's targets (rescored, §12 2026-09-28 B) | ≤ 5%% of draws clip | %.1f%% of %d (Wilson %.1f–%.1f%%) | %s |",
               100s.frac, length(s.ok), 100s.lo, 100s.hi, fires ? "**fires**" : "passes"))
    p(@sprintf("| K5, prior draws over every informed constant (14b, kept as the stress result) | ≤ 5%% of draws clip | %.1f%% of %d | %s |",
               100s14b.frac, length(s14b.ok), s14b.frac > K5_DRAW_FRACTION ? "**fires**" : "passes"))
    p()
end

if "broad" in modes
    p("## 14c.3: 14b's broad census, drawn Haldane-consistently")
    p()
    p("The same informed constants as 14b's census, at the same seed numbers, with every ",
      "reverse constant derived rather than drawn. The draws are not paired with 14b's: a ",
      "derived reverse constant consumes no random number, so the streams diverge after the ",
      "first. Recorded, not gated: the difference from 14b is what the broken equilibrium ",
      "constants did.")
    p()
    draw_section("broad", "14c.3, broad, Haldane-consistent")
end

if "gtp" in modes
    rows = last(data["gtp"])
    pub14b = Dict(r.seed => r for r in rows_14b if r.kind === :published)
    p("## 14c.6: the published-parameter census with the GTP branch loosened")
    p()
    p("PGK3's and PYK3's catalytic constants, forward and reverse, are scaled by s, ",
      "which scales each reaction's capacity by s at every handshake and keeps its ",
      "equilibrium constant. Recorded, not gated.")
    p()
    p("| scale | seeds carrying a deficit | drains clipped, median (max) | `GTP_translat` drains, median (max) | `tRNA_translat` drains, median (max) | failed |")
    p("|---|---|---|---|---|---|")
    cd(k, r) = sum((c.drains for c in r.clips if c.counter == k); init = 0)
    for sc in sort(unique(r.id for r in rows))
        rs = [r for r in rows if r.id == sc && r.status == "ok"]
        nf = count(r -> r.id == sc && r.status != "ok", rows)
        isempty(rs) && (p("| $sc | — | — | — | — | $nf |"); continue)
        g = [cd("GTP_translat", r) for r in rs]
        t = [cd("tRNA_translat", r) for r in rs]
        c = [r.clipped for r in rs]
        p(@sprintf("| %d | %d of %d | %.0f (%d) | %.0f (%d) | %.0f (%d) | %d |", sc,
                   count(>(0), c), length(rs), median(c), maximum(c), median(g), maximum(g),
                   median(t), maximum(t), nf))
    end
    p()
    s1 = [r for r in rows if r.id == 1 && r.status == "ok"]
    same = all(r -> haskey(pub14b, r.seed) && pub14b[r.seed].clipped == r.clipped &&
                    pub14b[r.seed].clips == r.clips, s1)
    p("At scale 1 the ", length(s1), " seeds' records ", same ? "are **identical**" : "**DIFFER**",
      " to 14b's published-parameter census.")
    p()
end

write(OUT, take!(out))
println("Wrote $OUT")
