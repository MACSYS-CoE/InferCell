# Shared helpers for the pilot and evidence scripts, from KIM-68 on. The M4 scripts
# that ran at `a742224` keep their own copies and are not edited. A script in
# dev/scripts includes it:
#   include(joinpath(@__DIR__, "pilot_common.jl"))
#
# - the run's identity: `COMMIT`, `JOB`, and `SMOKE` (`SMOKE=1` is a plumbing pass
#   at tiny sizes, never a result);
# - `V6`, B7's held-out `v6` fit, the posterior every pilot draws from;
# - `say`, timestamped and flushed;
# - `save` and `load`: .jls files carrying the commit, job, stage, threads and time;
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
const JOB = get(ENV, "SLURM_JOB_ID", "none")
const SMOKE = get(ENV, "SMOKE", "0") == "1"
const V6 = ["/fred/oz022/tkimpson/infercell-89-b7-archive-2026-10-27/run/v6/chain_$c.jls" for c in 1:2]

say(x...) = (println(now(), " ", x...); flush(stdout))

save(out, name, x; stage = "") =
    serialize(joinpath(out, "$name.jls"),
              (; result = x, commit = COMMIT, job = JOB, stage, smoke = SMOKE,
               threads = Threads.nthreads(), when = now()))
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
