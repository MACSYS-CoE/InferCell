# Shared helpers for the pilot and evidence scripts, from KIM-68 on. The M4 scripts
# that ran at `a742224` keep their own copies and are not edited. A script in
# dev/scripts includes it:
#   include(joinpath(@__DIR__, "pilot_common.jl"))
#
# - the run's identity: `COMMIT` (short) and `COMMIT_FULL`, `JOB`, and `SMOKE`
#   (`SMOKE=1` is a plumbing pass at tiny sizes, never a result);
# - `V6`, B7's held-out `v6` fit, the posterior every pilot draws from;
# - `say`, timestamped and flushed;
# - `save` and `load`: .jls files carrying the commit, job, stage, threads, GC threads,
#   time, and where they ran (`RAN_ON`: host, CPU, Julia version, array ids); a
#   warning at 16 threads or more without a concurrent GC sweep thread;
# - `moments` and `net`: μ̂, σ̂² and the count per pair, and the spread net of the
#   MC term;
# - `plain`: a `Prediction` as plain arrays, Dicts and NamedTuples. It reads back
#   at any commit, which a saved `Prediction` does not, since its fields change
#   (the M4 archives' open only at `a742224`). Save it beside or instead of one.

using InferCell
using Serialization
using Statistics: mean, var
using Dates

const COMMIT = strip(read(`git rev-parse --short HEAD`, String))
# The short form grows with the repo, so compare a saved COMMIT as its prefix.
const COMMIT_FULL = strip(read(`git rev-parse HEAD`, String))
const JOB = get(ENV, "SLURM_JOB_ID", "none")
const SMOKE = get(ENV, "SMOKE", "0") == "1"
const V6 = ["/fred/oz022/tkimpson/infercell-89-b7-archive-2026-10-27/run/v6/chain_$c.jls" for c in 1:2]
# GC's threads: `ngcthreads()` counts mark and sweep together (3 for both
# `--gcthreads=3` and `2,1`), so the sweep thread is read from `--gcthreads`, or
# without it from JULIA_NUM_GC_THREADS, which the options do not record.
# KIM-68: at 16 threads and more, run with one (`--gcthreads=N,1`).
function _sweep_threads()
    o = Base.JLOptions()
    o.nmarkthreads > 0 && return Int(o.nsweepthreads)
    s = split(get(ENV, "JULIA_NUM_GC_THREADS", ""), ',')
    return length(s) == 2 ? something(tryparse(Int, strip(s[2])), 0) : 0
end
const GC_THREADS = (; total = Threads.ngcthreads(), sweep = _sweep_threads())
!SMOKE && Threads.nthreads() >= 16 && GC_THREADS.sweep == 0 &&
    @warn "$(Threads.nthreads()) threads and no concurrent GC sweep thread: start Julia with --gcthreads=$(Threads.nthreads() ÷ 2),1 (KIM-68)"

say(x...) = (println(now(), " ", x...); flush(stdout))

# Where a record was made. The node's CPU decides what compiled code it loads
# (2026-10-23), so `cpu_name` is LLVM's name for it.
const RAN_ON = (; host = gethostname(), cpu = Sys.cpu_info()[1].model, cpu_name = Sys.CPU_NAME,
                julia = string(VERSION), array_job = get(ENV, "SLURM_ARRAY_JOB_ID", "none"),
                array_task = get(ENV, "SLURM_ARRAY_TASK_ID", "none"))

save(out, name, x; stage = "") =
    serialize(joinpath(out, "$name.jls"),
              (; result = x, commit = COMMIT, commit_full = COMMIT_FULL, job = JOB, stage,
               smoke = SMOKE, threads = Threads.nthreads(), gc_threads = GC_THREADS,
               ran_on = RAN_ON, when = now()))
load(out, name) = deserialize(joinpath(out, "$name.jls")).result

# μ̂, σ̂² and the cell count of a quantity read from each pair's finished cells.
function moments(p::Prediction, q::Symbol)
    xs = [[getproperty(c.accounts, q) for c in cs if c.failure === :none] for cs in vec(p.cells)]
    return [isempty(x) ? NaN : mean(x) for x in xs], [length(x) >= 2 ? var(x) : NaN for x in xs],
           length.(xs)
end

# The spread over the pairs, net of the MC term, for independent pairs.
net(mu, s2, n) = var(mu) - mean(s2 ./ n)

_sorted(d::AbstractDict) = sort!([k => v for (k, v) in d]; by = kv -> string(kv[1]))

# A cell as plain values; `minima` and `accounts` only when asked, since they are
# most of a cell's size.
plain(c::CellRun; full = false) =
    (; c.seed, c.q, c.failure, c.message, c.drains, n_clips = length(c.clips),
     minima = full ? _sorted(c.minima) : nothing,
     accounts = full && c.accounts !== nothing ? Dict(pairs(c.accounts)) : nothing)

"""
    plain(p::Prediction; cells = false, full = false) -> NamedTuple

`p`'s results and settings as plain values. The per-cell records are kept only
with `cells` (and their minima and carbon accounts only with `full`).
"""
function plain(p::Prediction; cells = false, full = false)
    pc(A) = cells ? [[plain(c; full) for c in cs] for cs in A] : nothing
    return (; p.forward, p.completion, p.estimand, p.quantity, p.horizon, p.draws, p.held,
            p.theta_H, p.seeds, p.viability_seeds, p.mu, p.sigma2, p.n_ok, p.effect,
            p.effect_sigma2, p.viable,
            decisions = [(; d.viable, d.fraction, counts = _sorted(d.counts)) for d in p.decisions],
            p.weights, p.conditioning, p.summary, p.effect_summary,
            statuses = _sorted(p.statuses), statuses_intervened = _sorted(p.statuses_intervened),
            delta = p.delta === nothing ? nothing :
                    (; p.delta.hybrid_mean, p.delta.hybrid_se, p.delta.n_ok, p.delta.delta,
                     p.delta.delta_intervened, mean_field_q = p.delta.mean_field.q),
            intervention = p.intervention === nothing ? nothing : string(p.intervention),
            viability = p.viability === nothing ? nothing : p.viability.label,
            settings = merge(p.provenance.settings,
                             (; intervention = string(p.provenance.settings.intervention))),
            fork_seed = p.provenance.fork_seed, version = string(p.provenance.version),
            threads = p.provenance.threads, cells = pc(p.cells),
            viability_cells = pc(p.viability_cells), cells_intervened = pc(p.cells_intervened),
            commit = COMMIT, job = JOB)
end
