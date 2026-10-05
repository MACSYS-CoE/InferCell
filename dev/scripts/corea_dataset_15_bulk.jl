# Parent task 15.7b (spec §12 2026-10-05): bulk metabolite observations for the
# 15.7 dataset.
#
# It reads the 15.7 run's latent record (job 17703841's task files, written by
# corea_dataset_15.jl) and checks it against a fresh simulation of two of its
# cells (30001 and the last). Then, for each panel pool at each save, it observes
# the 200-cell mean particle count lognormally at σ_b = 0.10 from a recorded seed
# (`observe_bulk`). It writes `bulk.tsv` (regenerable, ignored, SHA in the meta)
# and adds a bulk section to `meta.md`.
#
# Three idealisations are stated in the meta, not hidden:
# - the same 200 cells are averaged at every save;
# - the mean is of particle counts, not concentrations (the volume grows ~1.07×);
# - the transcripts stay exact per-cell counts every 60 s.
#
# Usage: julia --project dev/scripts/corea_dataset_15_bulk.jl <dir of the 15.7 run's task files>
# Regenerate: sbatch dev/scripts/corea_dataset_15_bulk.slurm <dir>

using InferCell
using Printf
using Serialization
using SHA

const SRC = ARGS[1]
const OUT = joinpath(@__DIR__, "..", "data", "corea_dataset_15")
const SIGMA_B = 0.10      # the bulk assay CV, ours (§12 2026-10-05)
const BULK_SEED = 15072
const NCELLS = 200

function main()
    files = sort(filter(f -> occursin(r"^task_\d+\.jls$", f), readdir(SRC)))
    parts = [deserialize(joinpath(SRC, f)) for f in files]
    p1 = first(parts)
    all(p -> p.ds.truth.values == p1.ds.truth.values && p.ds.species == p1.ds.species &&
             p.ds.times == p1.ds.times && p.commit == p1.commit && p.job == p1.job, parts) ||
        error("the task files disagree")
    p1.smoke && error("the task files are a smoke run")
    order = sortperm(reduce(vcat, [p.ds.seeds for p in parts]))
    seeds = reduce(vcat, [p.ds.seeds for p in parts])[order]
    latent_all = cat((p.ds.latent for p in parts)...; dims = 1)[order, :, :]
    length(seeds) == NCELLS && allunique(seeds) || error("expected $NCELLS distinct cells")

    # The record against a fresh simulation of two cells.
    for c in (1, NCELLS)
        fresh = generate_dataset(p1.ds.truth, [seeds[c]]; horizon = last(p1.ds.times))
        fresh.species == p1.ds.species || error("species order differs")
        fresh.latent[1, :, :] == latent_all[c, :, :] ||
            error("cell $(seeds[c]): the record is not what the model now produces")
        @printf("cell %d: the record equals a fresh simulation at every save\n", seeds[c])
    end

    keep = findall(>(0), p1.ds.times)
    latent = latent_all[:, :, keep]
    times = p1.ds.times[keep]
    bulk = observe_bulk(latent, p1.ds.species, M0_PANEL; sigma = SIGMA_B, noise_seed = BULK_SEED)
    bulk == observe_bulk(latent, p1.ds.species, M0_PANEL; sigma = SIGMA_B, noise_seed = BULK_SEED) ||
        error("two observations at one seed differ")
    rows = state_rows(M0_PANEL, p1.ds.species)
    xbar = dropdims(sum(latent[:, rows, :]; dims = 1); dims = 1) ./ NCELLS

    mkpath(OUT)
    tsv = joinpath(OUT, "bulk.tsv")
    open(tsv, "w") do io
        println(io, join(("t_s", "observable", "mean_latent", "observed"), '\t'))
        for (k, t) in enumerate(times), (j, o) in enumerate(M0_PANEL)
            println(io, join((Int(t), o, repr(xbar[j, k]), repr(bulk[j, k])), '\t'))
        end
    end
    digest = bytes2hex(open(sha256, tsv))
    head = strip(read(`git rev-parse --short HEAD`, String))

    meta = joinpath(OUT, "meta.md")
    text = read(meta, String)
    marker = "\n## Bulk observations (parent task 15.7b, §12 2026-10-05)\n"
    occursin(marker, text) && (text = first(split(text, marker)))
    io = IOBuffer()
    print(io, rstrip(text), "\n", marker, "\n")
    job = get(ENV, "SLURM_JOB_ID", "not under Slurm")
    println(io, "Written by `dev/scripts/corea_dataset_15_bulk.jl` at $head, job $job, Julia $VERSION, from the ",
            "latent record of the cells above (commit $(p1.commit), job $(p1.job)). Cells $(seeds[1]) ",
            "and $(seeds[end]) were resimulated, and each equals the record at every save.")
    println(io)
    println(io, "**Observation model.** For each of the $(length(M0_PANEL)) panel metabolites at each ",
            "of the $(length(times)) saves, one bulk measurement of the mean over the $NCELLS cells, ",
            "lognormal at σ_b = $SIGMA_B with the prediction floored at one particle, from ",
            "`Xoshiro($BULK_SEED)` (`observe_bulk`). σ_b is ours, a typical bulk assay CV, and is ",
            "the truth for σ_b, whose prior is `LogNormal(log 0.2, 1)`. The transcripts are ",
            "unchanged: exact per-cell counts. The per-cell panel above is the superseded ",
            "observation model, kept for reference.")
    println(io)
    println(io, "**Stated idealisations.**")
    println(io, "- The same $NCELLS cells are averaged at every save. A real bulk assay samples a fresh cohort.")
    println(io, "- The mean is of particle counts, not concentrations. The volume grows about 1.07× over ",
            "the cycle, so the two differ by under 7%.")
    println(io, "- The transcripts are exact per-cell counts every 60 s from the same living cell.")
    println(io)
    println(io, "`bulk.tsv` (t_s, observable, mean_latent, observed) is regenerable and not tracked. Its ",
            "SHA-256 is `$digest`; it depends on this Julia version's random streams.")
    write(meta, take!(io))
    @printf("wrote %s (SHA-256 %s) and the bulk section of meta.md\n", tsv, digest)
end

main()
