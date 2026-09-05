# Spec §11 task 3.8, re-run at task 5.7: measure what the 1 s handshake costs.
#
# Runs the phase-3 toy — one jump gene-expression module and one ODE metabolite
# module, exchanging state every simulated second — and reports seconds of
# wall-clock per second of simulated time, together with the extrapolation to a
# full 6,300 s cell cycle that kill criterion K1 is stated against.
#
# **The extrapolation is a floor, not an estimate of Core A′.** The toy is two
# modules, one gene and two metabolite states; Core A′ is seven modules,
# seventeen genes and thirty-two states. A pass here is necessary and not
# sufficient, and K1 is actually decided by task 13.7 at full scale. The report
# says so in its own text so the number cannot travel without the caveat.
#
# Usage: julia --project dev/scripts/bench_handshake.jl

using InferCell
using Random
using Printf
using Dates
using Statistics: median

include(joinpath(@__DIR__, "..", "..", "test", "contribution_test_models.jl"))
include(joinpath(@__DIR__, "..", "..", "test", "jump_test_models.jl"))
include(joinpath(@__DIR__, "..", "..", "test", "hybrid_test_models.jl"))

const CYCLE = 6300.0        # s, the published cell cycle
const BUDGET = 10.0         # s of wall-clock per trajectory, spec §8 K1

# Best of `reps` after one warm-up, matching dev/scripts/bench_rhs_containers.jl:
# the minimum is the honest estimator for a timing floor on a shared node.
#
# The *spread* is returned beside it, because without it a reader cannot tell a
# real per-handshake cost from run-to-run scatter, and at this scale the two are
# the same size. It is `(median - min) / min` rather than `(max - min) / min`:
# a single repetition stalled by another job on the node — the first version of
# this measurement recorded one at 577x the minimum — makes a max-based spread
# say nothing about the other six. The median is what a second run of this
# script would be expected to reproduce.
function time_run(build, n_steps; reps = 7)
    d = build()
    run_handshake!(d, min(n_steps, 10))          # warm-up: compile for these types
    times = Float64[]
    for _ in 1:reps
        d = build()
        push!(times, @elapsed run_handshake!(d, n_steps))
    end
    return minimum(times), median(times)
end

build_toy(; kcat = 0.0) = () -> begin
    Random.seed!(20260904)
    build_problem([ToyPool(kcat = kcat), ToyExpression()]; tspan = (0.0, CYCLE))
end

# Phase 5's volume chain adds work to every handshake — the membrane counts are
# summed, the geometry is recomputed, and a change in the factor rescales every
# ODE state. It is measured rather than assumed to be free, and against the same
# budget, because a per-handshake cost is exactly what K1 is a threshold on.
#
# It is measured against a **matched control**, not against the `live rate law`
# row above. That row starts at zero protein, and the protein count fills the
# ODE rate law through a CatalyticEdge, so comparing against it would attribute
# a different ODE trajectory to the volume chain. The control is the same gene
# at the same 831 initial copies with no flag and no edge, so the pair differs
# in the chain and in nothing else.
build_growing() = () -> begin
    Random.seed!(20260904)
    build_problem([ToyPool(kcat = 30.0), ToyGrowingExpression(k_tx = 2.0)];
                  tspan = (0.0, CYCLE))
end

build_control() = () -> begin
    Random.seed!(20260904)
    build_problem([ToyPool(kcat = 30.0), ToyExpression(protein0 = 831)];
                  tspan = (0.0, CYCLE))
end

# Provenance is captured *before* the measurements, not after. Reading
# `git rev-parse` at write time stamps whatever HEAD happens to be when the run
# finishes, which on a long job can be a commit made after the code was loaded —
# the artefact then names source the measurement never saw.
const COMMIT = strip(read(`git rev-parse --short HEAD`, String))
const DIRTY = !isempty(strip(read(`git status --porcelain -- src test Project.toml Manifest.toml dev/scripts/bench_handshake.jl`, String)))

rows = NamedTuple[]
for horizon in (60, 300, 600)
    configs = (("frozen pool (kcat = 0)", build_toy(kcat = 0.0)),
               ("live rate law", build_toy(kcat = 30.0)),
               ("831 copies, no growth (control)", build_control()),
               ("831 copies, growth live", build_growing()))
    for (label, build) in configs
        secs, med = time_run(build, horizon)
        per_sim_s = secs / horizon
        push!(rows, (horizon = horizon, label = label, wall = secs,
                     spread = (med - secs) / secs,
                     per_sim_s = per_sim_s, cycle = per_sim_s * CYCLE))
    end
end

# The verdict is taken on the slowest configuration measured, since that is the
# one a real module set most resembles.
worst = maximum(r -> r.cycle, rows)
verdict = worst <= BUDGET ? "PASS" : "FAIL"

commit, dirty = COMMIT, DIRTY
jobid = get(ENV, "SLURM_JOB_ID", "not under Slurm")

io = IOBuffer()
println(io, "# The 1 s handshake: wall-clock (spec task 3.8)")
println(io)
println(io, "Generated by `dev/scripts/bench_handshake.jl` at commit $commit",
        dirty ? " (WARNING: tree dirty)" : "",
        ", Slurm job $jobid, on $(gethostname()), $(Dates.now()), Julia $(VERSION).")
println(io)
println(io, "Regenerate with `sbatch dev/scripts/bench_handshake.slurm`.")
println(io)
println(io, "Toy: `ToyPool` (2 ODE states, stiff `Rodas5P`, abstol 1e-10, reltol 1e-8) ")
println(io, "composed with `ToyExpression` (3 jump states, 3 reactions), exchanging ")
println(io, "every 1.0 s under the `:fractional_carry` rounding policy. The last ")
println(io, "two rows are a matched pair: the same gene at the same 831 initial ")
println(io, "protein copies, with and without the membrane flag and volume edge, so ")
println(io, "the pair differs in phase 5's volume chain and in nothing else. Each ")
println(io, "row is the minimum of seven repetitions; `spread` is ")
println(io, "(median - min) / min over those, which is what says whether a gap ")
println(io, "between rows means anything. The median rather than the maximum, ")
println(io, "because one repetition stalled by another job on the node says ")
println(io, "nothing about the other six.")
println(io)
println(io, "| horizon (s) | configuration | wall-clock (s) | spread (median vs min, 7 reps) | s per simulated s | extrapolated to 6,300 s |")
println(io, "|---|---|---|---|---|---|")
for r in rows
    @printf(io, "| %d | %s | %.3e | %+.1f%% | %.3e | %.3f |\n",
            r.horizon, r.label, r.wall, 100 * r.spread, r.per_sim_s, r.cycle)
end

# The chain's own cost, taken only against its matched control, and reported
# beside the scatter that bounds how much of it is real.
println(io)
println(io, "The volume chain against its matched control, per horizon:")
println(io)
for horizon in (60, 300, 600)
    ctrl = only(filter(r -> r.horizon == horizon && r.label == "831 copies, no growth (control)", rows))
    grow = only(filter(r -> r.horizon == horizon && r.label == "831 copies, growth live", rows))
    @printf(io, "- %d s: %+.1f%%, against a within-configuration spread of %.1f%% (control) and %.1f%% (growth).\n",
            horizon, 100 * (grow.per_sim_s - ctrl.per_sim_s) / ctrl.per_sim_s,
            100 * ctrl.spread, 100 * grow.spread)
end
println(io)
@printf(io, "**Verdict against K1's %.0f s budget: %s** — worst extrapolation %.3f s per 6,300 s trajectory.\n",
        BUDGET, verdict, worst)
println(io)
println(io, """
**This number is a floor, not Core A′'s cost.** The toy carries one gene and two
metabolite states against Core A′'s seventeen and thirty-two, and its ODE block
is two equations rather than twenty-two reactions across four modules. The cost
of the real composition can only be larger. A pass here is therefore necessary
and not sufficient: it says the *architecture* does not carry a fatal per-
handshake overhead, which is what phase 3 exists to establish. K1 is decided by
task 13.7, at full scale.

What the measurement does bound honestly is the per-handshake overhead of the
exchange itself, which is the quantity the mechanism choice of task 3.1 was
made on: the frozen-pool row runs the same 1 s loop with a zero derivative, so
the gap between it and the live row is integration and the frozen row is very
nearly the handshake alone.

**Read the volume chain's cost only against its matched control, and only
beside the spread.** The control and growth rows are the same gene at the same
831 initial copies, differing in the flag and the edge alone; the `live rate
law` row starts at zero protein, and that count fills the ODE rate law through
a CatalyticEdge, so a gap against *it* would be a different ODE trajectory as
much as a volume chain. Even against the matched control the gap is small
enough that the within-configuration spread printed beside it is the thing to
check first: on a shared node these timings scatter by several percent, and a
best-of-five minimum does not remove that. The chain is measured on **two**
diluted ODE states against Core A′'s thirty-two, and its cost scales with that
count, so this bounds the mechanism rather than estimating the model.""")

s = String(take!(io))
print(s)
write(joinpath(@__DIR__, "bench_handshake_result.md"), s)
