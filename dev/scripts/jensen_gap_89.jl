# Issue #89 B3: the Jensen gap. How far is the mean-field description of gene
# expression from the population mean of the full two-way model, in the bulk
# metabolites? Throwaway prototype; gates B4 (spec/phases/approximate-method.md
# §8 KG, AV5).
#
# The mean-field cell: the full composition with transcription, decay and
# translation silenced (`Silent`: every state, parameter, edge and rebuild kept,
# no reactions), and B2's `MeanField` mean advanced by this script over each 1 s
# jump step at the driver's current rate constants. The driver runs everything
# else as published: the catalytic edges read the mean protein counts, the 1 s
# drains debit the mean costs (with the clip at zero), and the 60 s rebuild reads
# this cell's own pools, which are the mean-field pools, at 60m − 1 s. The ODE
# pools are written back continuously (`rounding = :continuous`).
#
# The jump state is integer, so the mean is written into it with fractional
# carry for the counters (exact in the drained total to one particle), and
# rounded to the nearest count for the proteins (the largest relative rounding
# error is reported). The transcripts are carried in the script; nothing in the
# silenced composition reads them.
#
# Reference: the 200-cell mean of the full model at the same truth, under
# fractional carry, cells seeded 40000 + k.
# - M0 (600 s): three truths, V6's replicate 16085 and V7's replicates 13 and 2
#   (ptsG 2.75, 1.13 and 33.5), ensembles run here.
# - Reduced syn3A (6,300 s): #88 A3's truths 1 and 2, whose `live60` ensembles
#   are reused from A3's archive. Its `src` differed from this script's base
#   only by `override_rebuild_pools!`, unset in `live60`.
#
# Score (§8 KG), per pool and in total over the M0_PANEL, at σ_b = 10%:
# Σ_t (log x̄_MF − log x̄_SSA)² / (2σ_b²), x̄ floored at one particle, against
# the noise Σ_t (log x̄_h1 − log x̄_h2)² / 2 / (2σ_b²) of the ensemble's halves,
# A3's statistic.
#
# Usage: julia --project dev/scripts/jensen_gap_89.jl ensemble <truth> <task> <ntasks> <outdir>
#        julia --project dev/scripts/jensen_gap_89.jl meanfield <truth> <outdir>
#        julia --project dev/scripts/jensen_gap_89.jl merge <outdir>
# Regenerate: dev/scripts/jensen_gap_89.sh <outdir outside the worktrees>

using InferCell
using InferCell: AbstractSubModel
using Random
using Serialization
using Statistics
using Printf
using Dates

const SMOKE = get(ENV, "SMOKE", "0") == "1"
const COMMIT = strip(read(`git rev-parse --short HEAD`, String))
const A3_ARCHIVE = get(ENV, "A3_ARCHIVE", "/fred/oz022/tkimpson/infercell-88-a3-archive-2026-10-19")
const NENSEMBLE = 200
const SIGMA_B = 0.10
# Each truth: the model, how its truth is drawn, its horizon, and where its
# reference ensemble comes from.
const TRUTHS = Dict(
    "m0_v6"   => (model = :m0, seed = 16_085, horizon = M0_HORIZON_S, ref = :here),
    "m0_low"  => (model = :m0, seed = 16_200_013, horizon = M0_HORIZON_S, ref = :here),
    "m0_high" => (model = :m0, seed = 16_200_002, horizon = M0_HORIZON_S, ref = :here),
    "corea1"  => (model = :corea, seed = 1507, horizon = COREA_CYCLE_S, ref = :a3),
    "corea2"  => (model = :corea, seed = 8802, horizon = COREA_CYCLE_S, ref = :a3),
)
const ORDER = ["m0_v6", "m0_low", "m0_high", "corea1", "corea2"]
horizon(t) = SMOKE ? 120.0 : TRUTHS[t].horizon

full_models(t) = TRUTHS[t].model === :m0 ? m0_models() : d11_models()

# M0's truths are its datasets' (m0_dataset draws the truth first); the reduced
# model's are A3's.
function truth_of(t)
    spec = TRUTHS[t]
    spec.model === :m0 ?
        draw_truth(Xoshiro(spec.seed), m0_models(); names = M0_TARGETS, purpose = :calibration).values :
        draw_truth(Xoshiro(spec.seed), d11_models(); purpose = :recovery).values
end

# A gene-expression module with its reactions removed and everything else kept.
struct Silent{M <: AbstractSubModel} <: AbstractSubModel
    m::M
end
for f in (:states, :parameters, :formalism, :inference_mode, :coupling, :inputs,
          :written_states, :rebuilt_params, :module_id, :membrane_protein_states,
          :contributed_states)
    @eval InferCell.$f(s::Silent) = InferCell.$f(s.m)
end
InferCell.rate_constants(p, t, s::Silent, pools) = InferCell.rate_constants(p, t, s.m, pools)
# One reaction that never fires, since the SSA cannot step a block with none.
InferCell.reactions(::Silent) = [InferCell.Reaction((u, p, t, w) -> 0.0, (u, w) -> nothing)]
InferCell.reduction_notes(s::Silent) =
    vcat(InferCell.reduction_notes(s.m),
         [:formalism => "$(module_id(s.m))'s reactions are replaced by the mean field (#89 B3)."])

const CME = (CoreATranscription, CoreATranscriptDecay, CoreATranslation)

ensemble_dir(out, t) = joinpath(out, "ensemble_$t")
meanfield_file(out, t) = joinpath(out, "meanfield_$t.jls")
slurm_job() = get(ENV, "SLURM_ARRAY_JOB_ID", get(ENV, "SLURM_JOB_ID", "none"))

panel_rows(species) = [findfirst(==(o), species) for o in M0_PANEL]

function run_ensemble(t, task, ntask, out)
    TRUTHS[t].ref === :here || error("truth $t's reference is A3's archive")
    dir = ensemble_dir(out, t); mkpath(dir)
    ms = full_models(t)
    truth = truth_of(t)
    cells, seeds = Matrix{Float64}[], Int[]
    species, times = Symbol[], Float64[]
    for k in 1:(SMOKE ? 4 : NENSEMBLE)
        (k - 1) % ntask == task || continue
        seed = 40_000 + k
        println("$(now()) ensemble $t seed $seed"); flush(stdout)
        Random.seed!(seed)
        d = build_problem(ms; tspan = (0.0, horizon(t)), complete = true)
        set_parameters!(d, ms, truth)
        Random.seed!(seed)
        r = emit_observables!(d, ms; every = 60)
        rows = panel_rows(r.species)
        species, times = r.species[rows], r.times
        push!(cells, r.counts[rows, :]); push!(seeds, seed)
    end
    serialize(joinpath(dir, "task_$(task).jls"),
              (; truth = t, values = truth, commit = COMMIT, smoke = SMOKE, job = slurm_job(),
               species, times, seeds, cells))
end

# One deterministic mean-field cell at truth `t`.
function run_meanfield(t, out)
    mkpath(out)
    full = full_models(t)
    truth = truth_of(t)
    cme = [m for m in full if any(T -> m isa T, CME)]
    length(cme) == 3 || error("expected the three gene-expression modules, found $(length(cme))")
    ms = AbstractSubModel[any(T -> m isa T, CME) ? Silent(m) : m for m in full]
    mf = MeanField(cme; edges = :none)
    d = build_problem(ms; tspan = (0.0, horizon(t)), complete = true, rounding = :continuous)
    set_parameters!(d, ms, truth)

    jump_models = [m for m in ms if formalism(m) === :jump]
    jnames = InferCell._block_names(ms, :jump)
    pnames = InferCell._block_param_names(jump_models)
    sidx = [findfirst(==(s), jnames) for s in states(mf)]
    pidx = [findfirst(==(p.name), pnames) for p in InferCell.model_free_params(parameters(mf))]
    (any(isnothing, sidx) || any(isnothing, pidx)) &&
        error("a mean-field state or parameter is missing from the driver")
    counters = Set(e.counter for m in ms for e in coupling(m) if e isa DeferredCounterEdge)
    is_counter = [s in counters for s in states(mf)]
    is_protein = [startswith(String(s), "P") for s in states(mf)]

    # The mean starts where the driver does; transcripts from their initial condition.
    x = Float64.(d.jump.u[sidx])
    carry = zeros(length(x))
    worst_round, worst_round_100, worst_state = 0.0, 0.0, :none
    ode_names = InferCell._block_names(ms, :ode)
    rows = panel_rows(ode_names)
    times, counts = Float64[], Vector{Float64}[]
    record!() = (push!(times, d.ode.t);
                 push!(counts, collect(d.ode.u[rows]) .* d.factor .+ d.rounding.remainders[rows]))
    record!()
    nsteps = round(Int, horizon(t) / d.interval)
    for k in 1:nsteps
        # The exchange (drain, rebuild, catalytic write, ODE step), then the
        # jump step's second, here the mean's, at the constants now in force.
        handshake_step!(d)
        p = d.jump.p[pidx]
        x1 = rk4(x, p, d.jump.t - d.interval, d.interval, mf)
        for i in eachindex(x)
            if is_counter[i]
                carry[i] += x1[i] - x[i]
                w = floor(carry[i])
                d.jump.u[sidx[i]] += Int(w)
                carry[i] -= w
            else
                n = round(Int, x1[i])
                d.jump.u[sidx[i]] = n
                if is_protein[i]
                    e = abs(x1[i] - n) / max(x1[i], 1.0)
                    e > worst_round && (worst_round = e; worst_state = states(mf)[i])
                    x1[i] >= 100 && (worst_round_100 = max(worst_round_100, e))
                end
            end
        end
        x = x1
        k % 60 == 0 && record!()
    end
    census = clipping_census(d)
    serialize(meanfield_file(out, t),
              (; truth = t, values = truth, commit = COMMIT, smoke = SMOKE, job = slurm_job(),
               species = M0_PANEL, times, counts = reduce(hcat, counts),
               worst_round, worst_round_100, worst_state, clipped = census.clipped, drains = census.drains,
               labels = [l.category for l in driver_declarations(d)]))
    @printf("%s mean field %s done: protein rounding at most %.2g (%s; %.2g at ≥ 100 copies), %d of %d drains clipped\n",
            now(), t, worst_round, worst_state, worst_round_100, census.clipped, census.drains)
end

# Four RK4 substeps over `h` of the linear mean equations at fixed parameters.
function rk4(x, p, t0, h, mf; n = 4)
    f(u, t) = collect(InferCell.dynamics(u, p, t, mf))
    s = h / n
    for i in 0:(n - 1)
        t = t0 + i * s
        k1 = f(x, t); k2 = f(x .+ s / 2 .* k1, t + s / 2)
        k3 = f(x .+ s / 2 .* k2, t + s / 2); k4 = f(x .+ s .* k3, t + s)
        x = x .+ s / 6 .* (k1 .+ 2k2 .+ 2k3 .+ k4)
    end
    return x
end

function load_parts(dir)
    files = sort(filter(f -> occursin(r"^task_\d+\.jls$", f), readdir(dir)))
    isempty(files) && error("no files in $dir")
    return [deserialize(joinpath(dir, f)) for f in files]
end

# The reference ensemble's cells (panel × times, particles), its times and its
# provenance.
function reference(out, t)
    if TRUTHS[t].ref === :here
        ps = load_parts(ensemble_dir(out, t))
        cells = vcat((p.cells for p in ps)...)
        return cells, first(ps).times, (; commit = only(unique(p.commit for p in ps)),
                                       job = only(unique(p.job for p in ps)), values = first(ps).values,
                                       smoke = only(unique(p.smoke for p in ps)))
    end
    a3 = Dict("corea1" => 1, "corea2" => 2)[t]
    ps = load_parts(joinpath(A3_ARCHIVE, "ensemble_live60_truth$a3"))
    sp = first(ps).species
    rows = panel_rows(sp)
    cells = [c[rows, :] for p in ps for c in p.cells]
    return cells, first(ps).times, (; commit = only(unique(p.commit for p in ps)),
                                   job = only(unique(p.job for p in ps)), values = first(ps).values,
                                   smoke = only(unique(p.smoke for p in ps)))
end

logbulk(cells, j) = [log(max(mean(c[i, j] for c in cells), 1.0)) for i in eachindex(M0_PANEL)]

function merge(out)
    io = IOBuffer()
    p(args...) = (println(io, args...); println(args...))
    p("# Issue #89 B3: the Jensen gap", SMOKE ? " — SMOKE RUN, NOT A RESULT" : "")
    p()
    p("Merged at $COMMIT. σ_b = $(Int(100SIGMA_B))%. Scores in nats; they scale as 1/σ_b².")
    p()
    summary = Dict{String, Any}()
    for t in ORDER
        isfile(meanfield_file(out, t)) || (p("## $t: mean field missing\n"); continue)
        mfr = deserialize(meanfield_file(out, t))
        mfr.smoke == SMOKE || error("$t: the mean field's smoke flag differs from this merge's")
        cells, rtimes, prov = reference(out, t)
        prov.values == mfr.values || error("$t: the reference ensemble's truth differs from the mean field's")
        # Saves after t = 0 that both have.
        common = [τ for τ in mfr.times if τ > 0 && τ in rtimes]
        jr = [findfirst(==(τ), rtimes) for τ in common]
        jm = [findfirst(==(τ), mfr.times) for τ in common]
        h = length(cells) ÷ 2
        L = reduce(hcat, [logbulk(cells, j) for j in jr])
        H1 = reduce(hcat, [logbulk(cells[1:h], j) for j in jr])
        H2 = reduce(hcat, [logbulk(cells[(h + 1):(2h)], j) for j in jr])
        M = log.(max.(mfr.counts[:, jm], 1.0))
        shift = vec(sum((M .- L) .^ 2; dims = 2)) ./ (2SIGMA_B^2)
        noise = vec(sum((H1 .- H2) .^ 2; dims = 2)) ./ 2 ./ (2SIGMA_B^2)
        worst = vec(maximum(abs.(M .- L); dims = 2)) ./ SIGMA_B
        pass_total = sum(shift) <= sum(noise)
        bad = [M0_PANEL[i] for i in eachindex(M0_PANEL) if shift[i] > 3noise[i]]
        verdict = pass_total && isempty(bad)
        summary[t] = (; shift, noise, worst, pass_total, bad, verdict, ncells = length(cells),
                      nsaves = length(common), worst_round = mfr.worst_round,
                      worst_round_100 = mfr.worst_round_100, clipped = mfr.clipped, drains = mfr.drains)
        spec = TRUTHS[t]
        p("## $t ($(spec.model === :m0 ? "M0" : "reduced syn3A"), truth seed $(spec.seed))")
        p()
        p(join([@sprintf("%s %.4g", k, v) for (k, v) in mfr.values], ", "), ".")
        p()
        p("- Reference: $(length(cells)) cells at $(prov.commit) (job $(prov.job))",
          spec.ref === :a3 ? ", #88 A3's `live60`" : "", ". Mean field at $(mfr.commit) (job $(mfr.job)).")
        p("- $(length(common)) saves. Protein rounding at most ", @sprintf("%.2g", mfr.worst_round),
          " relative, in `$(mfr.worst_state)` (", @sprintf("%.2g", mfr.worst_round_100), " at ≥ 100 copies). ",
          "$(mfr.clipped) of $(mfr.drains) mean-field drains clipped.")
        p("- **Verdict (KG): ", verdict ? "pass" : "FAIL", "**. Total ",
          @sprintf("%.3g", sum(shift)), " against noise ", @sprintf("%.3g", sum(noise)),
          isempty(bad) ? "; no pool above 3× its noise." : "; above 3× their noise: " * join(string.("`", bad, "`"), ", ") * ".")
        p()
        p("| Pool | noise | shift | max \\|Δ log x̄\\| / σ_b |")
        p("|---|---|---|---|")
        for (i, o) in enumerate(M0_PANEL)
            p("| `$o` | ", @sprintf("%.3g", noise[i]), " | ",
              shift[i] > 3noise[i] ? @sprintf("**%.3g**", shift[i]) : @sprintf("%.3g", shift[i]),
              " | ", @sprintf("%.2f", worst[i]), " |")
        end
        p("| **total** | ", @sprintf("**%.3g**", sum(noise)), " | ", @sprintf("**%.3g**", sum(shift)), " | |")
        p()
    end
    write(joinpath(out, "merged.md"), take!(io))
    serialize(joinpath(out, "merged.jls"), summary)
    println("Wrote ", joinpath(out, "merged.md"))
end

function main()
    stage = ARGS[1]
    if stage == "ensemble"
        run_ensemble(ARGS[2], parse(Int, ARGS[3]), parse(Int, ARGS[4]), ARGS[5])
    elseif stage == "meanfield"
        run_meanfield(ARGS[2], ARGS[3])
    elseif stage == "merge"
        merge(ARGS[2])
    else
        error("stage must be ensemble, meanfield or merge")
    end
end

main()
