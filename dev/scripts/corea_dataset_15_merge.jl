# Merge task 15.7's dataset run (corea_dataset_15.jl) into
# dev/data/corea_dataset_15/. The panel comes from 15.6's merged sensitivity
# (dev/scripts/sensitivity_15/merged.jls). Metabolites in it are observed
# lognormally at the truth's noise scale, and the 17 transcripts are kept as
# exact counts (§4 D8, D11).
#
# Writes `observations.tsv` (cell seed, time, observable, latent, observed), which
# is regenerable and gitignored, and `meta.md`, which is tracked. `meta.md`
# carries the truth, the seeds, the noise scale and seed, the run's commit and
# job, the file's SHA-256, and the reduction_report of the model that produced it.
#
# The lognormal floor is one particle, the resolution of the count round trip:
# a pool below one particle cannot be observed as a concentration. The t = 0
# save point is the fixed initial condition, and is left out (§12 2026-09-29).
#
# Regeneration order, all in one tree: the sensitivity run
#   sbatch --array=0-1469 --export=ALL,NSEEDS=1470 dev/scripts/corea_sensitivity_15.slurm
# then corea_sensitivity_15_merge.jl, which writes the panel to merged.jls; then
#   sbatch dev/scripts/corea_dataset_15.slurm
# and then this script.
#
# Usage: julia --project dev/scripts/corea_dataset_15_merge.jl

using InferCell
using Printf
using Serialization
using SHA

const DIR = joinpath(@__DIR__, "dataset_15")
const SENS = joinpath(@__DIR__, "sensitivity_15", "merged.jls")
const OUT = joinpath(@__DIR__, "..", "data", "corea_dataset_15")
const SIGMA_MET = 0.1      # the truth's metabolite noise scale, ours
const NOISE_SEED = 1508
const NCELLS = 200
const OBS_FLOOR = 1.0      # particles

files = sort(filter(f -> occursin(r"^task_\d+\.jls$", f), readdir(DIR)))
parts = [deserialize(joinpath(DIR, f)) for f in files]
p1 = first(parts)
all(p -> p.ds.truth.values == p1.ds.truth.values && p.ds.species == p1.ds.species &&
         p.ds.times == p1.ds.times, parts) || error("tasks disagree on truth or layout")
all(p -> p.commit == p1.commit && p.job == p1.job && p.smoke == p1.smoke, parts) ||
    error("the task files come from more than one run")
order = sortperm(reduce(vcat, [p.ds.seeds for p in parts]))
seeds = reduce(vcat, [p.ds.seeds for p in parts])[order]
latent = cat((p.ds.latent for p in parts)...; dims = 1)[order, :, :]
p1.smoke || length(seeds) == NCELLS || error("read $(length(seeds)) cells, not $NCELLS")
allunique(seeds) || error("a cell seed repeats")

sens = deserialize(SENS)
panel = sens.panel
p1.smoke || !sens.smoke || error("the panel comes from a smoke sensitivity run")
keep_t = findall(>(0), p1.ds.times)      # t = 0 is the fixed initial condition
latent = latent[:, :, keep_t]
times = p1.ds.times[keep_t]
tx = [transcript_state(g) for g in p1.ds.genes]
noise = NoiseModel(Modality(:metabolite, :lognormal, panel; floor = OBS_FLOOR),
                   Modality(:transcript, :poisson, tx))
observed, osp = observe_latent(latent, p1.ds.species, noise;
                               scales = Dict(:sigma_metabolite => SIGMA_MET),
                               noise_seed = NOISE_SEED)
assert_no_circularity(d11_models(), osp)

mkpath(OUT)
tsv = joinpath(OUT, "observations.tsv")
rows = state_rows(osp, p1.ds.species)
open(tsv, "w") do io
    println(io, join(("cell_seed", "t_s", "observable", "latent", "observed"), '\t'))
    for (c, s) in enumerate(seeds), (k, t) in enumerate(times), (j, o) in enumerate(osp)
        println(io, join((s, Int(t), o, repr(latent[c, rows[j], k]), repr(observed[c, j, k])), '\t'))
    end
end
digest = bytes2hex(open(sha256, tsv))

io = IOBuffer()
p(args...) = println(io, args...)
p("# The 15.7 synthetic dataset", p1.smoke ? " — SMOKE RUN, NOT A RESULT" : "")
p()
p("Cells run at commit $(p1.commit), job $(p1.job), by `dev/scripts/corea_dataset_15.jl`; ",
  "merged at $(strip(read(`git rev-parse --short HEAD`, String))) by `corea_dataset_15_merge.jl`, ",
  "Julia $VERSION. The panel comes from the sensitivity run at commit $(sens.commit), job ",
  "$(sens.job), $(sens.nseeds) seeds, merged at $(sens.merged_at). `observations.tsv` is ",
  "regenerable and not tracked. Its SHA-256 is `$digest`; the SHA depends on this Julia ",
  "version's random streams. Regeneration order: the sensitivity run and its merge, then ",
  "`sbatch dev/scripts/corea_dataset_15.slurm`, then `corea_dataset_15_merge.jl`.")
p()
p("**Model.** The published model: clamped drain, fractional carry, with D11's six freed ",
  "(`d11_models()`). $(length(seeds)) cells, seeds $(first(seeds)) to $(last(seeds)), each a ",
  "$(Int(last(p1.ds.times))) s cycle. Observations are at the $(length(times)) save points ",
  "from $(Int(first(times))) s to $(Int(last(times))) s, every 60 s; t = 0 is the fixed initial ",
  "condition and is left out. Values are in particles, as are the latent states.")
p()
p("**Truth.** Drawn from the prior by `draw_truth(Xoshiro($(p1.truth_seed)), d11_models(); ",
  "purpose = :recovery)`, Haldane-consistently (14c.1). The six targets, and the reverse ",
  "constants they derive:")
p()
p("| parameter | truth |")
p("|---|---|")
for (n, v) in p1.ds.truth.values
    p("| `$n` | $(@sprintf("%.6g", v)) |")
end
p()
p("**Observation model.** $(length(panel)) metabolites (15.6's panel) observed lognormally at ",
  "σ = $SIGMA_MET, from `Xoshiro($NOISE_SEED)`, with the prediction floored at $(OBS_FLOOR) ",
  "particle: a pool below one particle is not observable as a concentration, and a clipped ",
  "drain can leave a published pool at or just below zero. σ = $SIGMA_MET is ours, and is ",
  "itself a truth for check 9. The 17 transcripts are exact counts. The `:poisson` modality ",
  "only labels them as counts: phase 16 scores them by D10's exact conditional, not as ",
  "Poisson draws.")
p()
p("Panel: ", join(("`$m`" for m in panel), ", "), ".")
p()
p("**T2 row for the truth rule** (`$(p1.ds.truth.label.category)`, `$(p1.ds.truth.label.subject)`): ",
  p1.ds.truth.label.description, ".")
p()
p("## reduction_report of the model that produced the data")
p()
p("```")
p(p1.ds.report)
p("```")
write(joinpath(OUT, "meta.md"), take!(io))
println("Wrote $OUT")
