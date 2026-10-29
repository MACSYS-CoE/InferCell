# Issue #89 B8.4: profile the mixed likelihood on Core A′ and project AV7's cost
# (spec/phases/approximate-method.md §11 B8.4, §12 2026-10-29 "B8 as planned").
#
# The pilot dataset is corea_dataset(17_200_000), outside AV7's 200 (r = 1–200), at
# full size. δ and S come from NREF reference cells at the dataset's truth, enough
# to give the bulk term its full 1785 × 1785 shape; the timings do not depend on
# their values. B7.2's snapshot is the comparison: 2.143 s per gradient on m0_v6.
#
# Usage: julia --project dev/scripts/profile_89_b8.jl data <out dir>
#        julia --project dev/scripts/profile_89_b8.jl snapshot <out dir> <file.jls>
#
# `data` writes <out>/pilot.jls: the dataset, the reference cells, δ, S and the
# seconds each took. `snapshot` writes the log density and its gradient at the
# truth and the prior median, the median of three timings of each, and the cost
# by part, as profile_89_b7.jl does.

using InferCell
using InferCell: _mixed_solve
using Serialization
using Statistics
using LinearAlgebra
using Printf
using ForwardDiff

const COMMIT = strip(read(`git rev-parse --short HEAD`, String))
const PILOT_SEED = 17_200_000
const NREF = parse(Int, get(ENV, "NREF", "50"))
job() = get(ENV, "SLURM_JOB_ID", "none")

corea_setup(ds; kw...) = mixed_setup(ds; models = d11_models(), targets = D11_TARGETS, kw...)

truth_x(ds, s) = (tv = Dict(ds.truth.values);
                  vcat([log(tv[n]) for n in s.names], log(ds.sigma)))

function run_data(out)
    mkpath(out)
    t_ds = @elapsed ds = corea_dataset(PILOT_SEED)
    s0 = corea_setup(ds)
    x = truth_x(ds, s0)
    values = mixed_theta_values(s0, x)
    m_ref = mean_field_log_bulk(s0, _mixed_solve(s0, exp.(x[1:length(s0.names)])).sol)
    seeds = reference_seeds(ds, NREF)
    t_ens = @elapsed cells = reference_ensemble(values, seeds; models = d11_models(),
                                                horizon = COREA_CYCLE_S)
    d = delta_and_S(cells, m_ref; C = size(ds.transcripts, 1))
    serialize(joinpath(out, "pilot.jls"),
              (; ds, x, values, m_ref, seeds, cells, delta = d.delta, S = d.S, commit = COMMIT,
               job = job(), seconds = (dataset = t_ds, ensemble = t_ens)))
    @printf("data at %s (job %s): %d cells in %.0f s (%.1f s per cell); %d reference cells in %.0f s (%.1f s per cell)\n",
            COMMIT, job(), length(ds.seeds), t_ds, t_ds / length(ds.seeds), NREF, t_ens, t_ens / NREF)
end

med_time(f, k = 3) = (f(); ts = [@timed(f()) for _ in 1:k];
                      (seconds = median(t.time for t in ts), bytes = median(t.bytes for t in ts)))

function breakdown(s, x)
    n = length(s.names)
    θ = exp.(x[1:n])
    tag = ForwardDiff.Tag(:profile, Float64)
    θd = [ForwardDiff.Dual{typeof(tag)}(θ[i], ntuple(j -> Float64(i == j), n + 1)) for i in 1:n]
    primal = _mixed_solve(s, θ)
    dual = _mixed_solve(s, θd)
    k = route_transcription_rates(s, primal.sol, primal.p)
    μ = θ[findfirst(==(:krnadeg), s.names)] ./ s.lengths
    kd = route_transcription_rates(s, dual.sol, dual.p)
    μd = θd[findfirst(==(:krnadeg), s.names)] ./ s.lengths
    m = mean_field_log_bulk(s, primal.sol)
    md = mean_field_log_bulk(s, dual.sol)
    σ = exp(x[end])
    σd = ForwardDiff.Dual{typeof(tag)}(σ, ntuple(j -> Float64(j == n + 1), n + 1))
    stats(sol) = (; nf = sol.stats.nf, njacs = sol.stats.njacs, nw = sol.stats.nw,
                  naccept = sol.stats.naccept, nreject = sol.stats.nreject)
    return (; solve_primal = med_time(() -> _mixed_solve(s, θ)),
            solve_dual = med_time(() -> _mixed_solve(s, θd)),
            stats_primal = stats(primal.sol), stats_dual = stats(dual.sol),
            tx_primal = med_time(() -> transcript_loglik(s.factors, k, μ)),
            tx_dual = med_time(() -> transcript_loglik(s.factors, kd, μd)),
            bulk_primal = med_time(() -> bulk_loglik(s, m, σ)),
            bulk_dual = med_time(() -> bulk_loglik(s, md, σd)))
end

function run_snapshot(out, file)
    p = deserialize(joinpath(out, "pilot.jls"))
    s = with_reference(corea_setup(p.ds), p.delta, p.S)
    xs = [p.x, mixed_prior_median(s)]
    ℓ(z) = mixed_logdensity(s, z)
    vals = [ℓ(x) for x in xs]
    grads = [ForwardDiff.gradient(ℓ, x) for x in xs]
    density = med_time(() -> ℓ(xs[1]))
    gradient = med_time(() -> ForwardDiff.gradient(ℓ, xs[1]))
    parts = breakdown(s, xs[1])
    serialize(file, (; commit = COMMIT, job = job(), xs, vals, grads, density, gradient, parts,
                     threads = Threads.nthreads(), data_seconds = p.seconds,
                     ncells = length(p.ds.seeds), nref = length(p.seeds)))
    @printf("snapshot at %s (job %s): values %s\n", COMMIT, job(), vals)
    @printf("density %.3f s (%.2f GiB), gradient %.3f s (%.2f GiB)\n", density.seconds,
            density.bytes / 2^30, gradient.seconds, gradient.bytes / 2^30)
    for (k, v) in pairs(parts)
        println("  ", k, ": ", v isa NamedTuple && haskey(v, :seconds) ?
                @sprintf("%.4f s, %.3f GiB", v.seconds, v.bytes / 2^30) : v)
    end
end

function main(args)
    if args[1] == "data"
        run_data(args[2])
    elseif args[1] == "snapshot"
        run_snapshot(args[2], args[3])
    else
        error("unknown stage $(args[1])")
    end
end

main(ARGS)
