# Merge task 14c.5's 5,000-seed run into its result (spec §11 task 14c.5; §12
# 2026-09-29, the later entry). Reads dev/scripts/smoothing_14c5/, which
# `corea_smoothing_14c5.jl` writes, and writes
# dev/scripts/corea_smoothing_14c5_result.md.
#
# - Coupled half: every paired difference at a save point before a pair's first
#   jump-state mismatch within 1%, relative to max(|published|, floor). The
#   driver computes it in the run. This script takes the largest.
# - Decoupled half: `ensemble_agreement` per observable at the end of the cycle.
#   An observable passes when |d| + 2·SE ≤ max(R, 1%), where R is the 200-cell
#   dataset's resolution on it.
# - Null: the published runs split into two halves by seed order and compared
#   unpaired. The new statistic must put no observable beyond |z| = 3.29 (5%,
#   Bonferroni over 50). The old per-seed ratio is reported beside it.
#
# Usage: julia --project dev/scripts/corea_smoothing_14c5_merge.jl

using InferCell
using Printf
using Statistics: mean, std, median

const DIR = joinpath(@__DIR__, "smoothing_14c5")
const OUT = joinpath(@__DIR__, "corea_smoothing_14c5_result.md")
const THRESHOLD = 0.01
const NCELLS = 200
const ZNULL = 3.29
const NSEEDS = 5000
# The driver's floors, which its coupled half used; the two must not drift.
const ODE_FLOOR = 500.0
const TX_FLOOR = 1.0
floor_of(o) = o === :volume_litres ? 0.0 :
              startswith(String(o), "mRNA_") ? TX_FLOOR : ODE_FLOOR
const MERGED_AT = strip(read(`git rev-parse --short HEAD`, String))

files = sort(filter(f -> occursin(r"^task_\d+\.tsv$", f), readdir(DIR)))
isempty(files) && error("no task TSVs in $DIR")
headers = String[]
census = NamedTuple[]
couple = NamedTuple[]
vals = Dict{Tuple{Int, Symbol, Symbol}, Float64}()
for f in files, line in eachline(joinpath(DIR, f))
    c = split(line, '\t')
    if startswith(line, "# census")
        push!(census, (seed = parse(Int, c[2]), variant = Symbol(c[3]),
                       drains = parse(Int, c[4]), clipped = parse(Int, c[5]),
                       wall = parse(Float64, c[6])))
    elseif startswith(line, "# couple")
        push!(couple, (seed = parse(Int, c[2]), t_dec = parse(Int, c[3]),
                       ncomp = parse(Int, c[4]), r = parse(Float64, c[5]),
                       t = parse(Int, c[6]), o = Symbol(c[7])))
    elseif startswith(line, "#")
        push!(headers, line[3:end])
    elseif !startswith(line, "seed\t")
        vals[(parse(Int, c[1]), Symbol(c[2]), Symbol(c[3]))] = parse(Float64, c[4])
    end
end
smoke = any(h -> occursin("SMOKE", h), headers)
seeds = sort([x.seed for x in couple])
length(unique(seeds)) == length(seeds) || error("a seed appears twice")
smoke || length(seeds) == NSEEDS ||
    error("read $(length(seeds)) seeds, not $NSEEDS: a task is missing or incomplete")
obs = sort(unique(k[3] for k in keys(vals)))
col(v, o, ss = seeds) = [vals[(s, v, o)] for s in ss]

out = IOBuffer()
p(args...) = (println(out, args...); println(args...))
p("# 14c.5: the sampler model against the published one, over $(length(seeds)) seeds",
  smoke ? " — SMOKE RUN, NOT A RESULT" : "")
p()
p("Merged by `dev/scripts/corea_smoothing_14c5_merge.jl` from $(length(files)) task files. ",
  "Regenerate with `sbatch dev/scripts/corea_smoothing_14c5.slurm`, then this script. ",
  "The gate is §12 2026-09-29's. Merged at commit $MERGED_AT. The tasks' headers:")
p()
for h in unique(replace.(headers, r", task \d+ of \d+, [0-9T:.-]+" => ""))
    p("- ", h)
end
p()
p("Each seed is a full cycle under the published model (clamped drain, fractional carry) ",
  "and the sampler model (drain smoothed at one particle, pools carried continuously). ",
  "$(length(obs)) observables: every ODE state in particles with its carried remainder, ",
  "every transcript, and the volume. Floors: 500 particles for an ODE state, one copy for a ",
  "transcript, none for the volume.")
p()

# Coupled half.
dec = [x.t_dec for x in couple if x.t_dec >= 0]
w = couple[argmax([abs(x.r) for x in couple])]
cpass = abs(w.r) <= THRESHOLD
p("## Coupled half")
p()
p(@sprintf("%d of %d pairs decouple within the cycle", length(dec), length(seeds)),
  isempty(dec) ? ". " : @sprintf(", at a median handshake of %.0f (earliest %d). ",
                                 median(dec), minimum(dec)),
  @sprintf("Before decoupling there are %d paired comparisons. The largest is %.3g%% ",
           sum(x.ncomp for x in couple), 100w.r),
  "(`$(w.o)`, seed $(w.seed), t = $(w.t) s). Gate: every one within 1% — **",
  cpass ? "passes" : "fails", "**.")
p()

# Decoupled half.
st = [(o = o, a = ensemble_agreement(col(:published, o), col(:sampler, o);
                                      floor = floor_of(o), ncells = NCELLS,
                                      threshold = THRESHOLD)) for o in obs]
npass = count(x -> x.a.verdict === :pass, st)
nfail = count(x -> x.a.verdict === :fail, st)
nunres = count(x -> x.a.verdict === :unresolved, st)
p("## Decoupled half, at the end of the cycle")
p()
p("Per observable: d is the difference of the ensemble means over max(|published mean|, floor), ",
  "SE is from the paired differences, and R is the $NCELLS-cell dataset's standard error on ",
  "the mean. The observable passes when |d| + 2·SE ≤ max(R, 1%). ",
  "**$npass pass, $nfail fail, $nunres unresolved.**")
p()
p("| observable | d | SE | R | (\\|d\\| + 2·SE) / tolerance | verdict |")
p("|---|---|---|---|---|---|")
for x in sort(st; by = x -> -(abs(x.a.diff) + 2x.a.se) / x.a.tolerance)
    a = x.a
    p(@sprintf("| `%s` | %+.3g%% | %.2g%% | %.2g%% | %.2f | %s |", x.o, 100a.diff, 100a.se,
               100a.resolution, (abs(a.diff) + 2a.se) / a.tolerance,
               a.verdict === :pass ? "passes" : a.verdict === :fail ? "**fails**" : "unresolved"))
end
p()

# Null.
h = length(seeds) ÷ 2
h >= 2 || @warn "the null needs two seeds per half; smoke output only"
sa, sb = seeds[1:h], seeds[h+1:2h]
nul = NamedTuple[]
for o in obs
    h >= 2 || break
    x, y = col(:published, o, sb), col(:published, o, sa)
    std(x) == 0 && std(y) == 0 && continue
    a = ensemble_agreement(y, x; floor = floor_of(o), ncells = NCELLS, threshold = THRESHOLD)
    r = (x .- y) ./ max.(abs.(y), floor_of(o))
    oldse = std(r) / sqrt(h)
    push!(nul, (o = o, z = a.se > 0 ? a.diff / a.se : 0.0, diff = a.diff, se = a.se,
                old = mean(r), oldse = oldse, oldz = mean(r) / max(oldse, eps())))
end
zmax = isempty(nul) ? NaN : maximum(abs(x.z) for x in nul)
npass_null = zmax <= ZNULL
p("## The null: the published model against itself on independent seeds")
p()
p("The published runs at the first $h seeds against those at the other $h, unpaired. ",
  "$(length(nul)) observables vary across seeds. For the new statistic the largest |z| is ",
  @sprintf("%.2f", zmax), ", against $(ZNULL) (5%, Bonferroni over 50) — **",
  npass_null ? "passes" : "fails", "**. The old per-seed ratio, on the same pairs:")
p()
p("| observable | old: mean per-seed ratio | its z | new: d | its z |")
p("|---|---|---|---|---|")
for x in sort(nul; by = x -> -abs(x.oldz))[1:min(12, end)]
    p(@sprintf("| `%s` | %+.3g%% ± %.2g%% | %.2f | %+.3g%% ± %.2g%% | %.2f |", x.o,
               100x.old, 100x.oldse, x.oldz, 100x.diff, 100x.se, x.z))
end
p()

# Census, recorded not compared (see 14c's result for why the predicates differ).
p("## Clipping and wall-clock")
p()
p("| | median drains clipped | seeds carrying a deficit | median wall-clock per cycle (s) |")
p("|---|---|---|---|")
for v in (:published, :sampler)
    cs = [x for x in census if x.variant === v]
    p("| `$v` | $(median([x.clipped for x in cs])) | $(count(x -> x.clipped > 0, cs)) of ",
      "$(length(cs)) | $(@sprintf("%.1f", median([x.wall for x in cs]))) |")
end
p()

p("## Verdict")
p()
closed = cpass && nfail == 0 && nunres == 0 && npass_null
p("Coupled half: ", cpass ? "passes" : "**fails**", ". Decoupled half: $npass of ",
  "$(length(obs)) pass, $nfail fail, $nunres unresolved. Null: ",
  npass_null ? "passes" : "**fails**", ". 14c.5 is ",
  closed ? "**closed**." : "**not closed**.")
cpass || (p(); p("The coupled half's failure is recorded in spec §12 (2026-09-29, the entry ",
                 "on the coupled half) and carried as the measured cost of continuous pools; ",
                 "`corea_smoothing_14c5_decomp_result.md` attributes it. 14c.5 closes on that record."))
write(OUT, take!(out))
println("Wrote $OUT")
