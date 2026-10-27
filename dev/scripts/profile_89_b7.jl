# Issue #89 B7.2: profile the mixed likelihood's gradient on m0_v6, and check that
# a speedup leaves the numbers unchanged (spec/phases/approximate-method.md §11 B7.2,
# §12 2026-10-24 "B7 as planned").
#
# Reads B5's archive (m0_v6's dataset, θ_ref and the 800-cell ensemble), so the
# setup carries B5's δ and S.
#
# Usage: julia --project dev/scripts/profile_89_b7.jl snapshot <out.jls>
#        julia --project dev/scripts/profile_89_b7.jl compare <before.jls> <after.jls>
#
# `snapshot` writes the log density and its gradient at five fixed x (θ_ref, the
# prior median and three of B5's chain-1 draws), the median of three timings of
# each, and the cost broken down by part: route (ii)'s solve in Float64 and in
# six-partial duals, with their `sol.stats`, the transcript factors and the bulk
# term. `compare` holds the second snapshot's values to the first's at 1e-10
# relative and prints the timings side by side.

using InferCell
using Serialization
using Statistics
using LinearAlgebra
using Printf
using ForwardDiff
using MCMCChains

const B5 = get(ENV, "B5_ARCHIVE", "/fred/oz022/tkimpson/infercell-89-b5-archive-2026-10-23")
const COMMIT = strip(read(`git rev-parse --short HEAD`, String))
const RTOL = 1e-10

function b5_setup()
    ds = deserialize(joinpath(B5, "dataset.jls"))
    ref = deserialize(joinpath(B5, "theta_ref.jls"))
    chunks = [deserialize(joinpath(B5, "ensemble_$k.jls")) for k in 1:16]
    cells = cat((c.cells for c in chunks)...; dims = 1)
    d = delta_and_S(cells, ref.m_ref; C = size(ds.transcripts, 1))
    s = with_reference(mixed_setup(ds), d.delta, d.S)
    return s, ref
end

function points(s, ref)
    ch = deserialize(joinpath(B5, "chain_1.jls")).chain
    n = length(mixed_parameter_names(s))
    X = Array(ch[:, [Symbol("x[$i]") for i in 1:n], 1])
    return [ref.x, mixed_prior_median(s), X[1, :], X[150, :], X[300, :]]
end

med_time(f, k = 3) = (f(); ts = [@timed(f()) for _ in 1:k];
                      (seconds = median(t.time for t in ts), bytes = median(t.bytes for t in ts)))

function breakdown(s, x)
    n = length(s.names)
    θ = exp.(x[1:n])
    tag = ForwardDiff.Tag(:profile, Float64)
    θd = [ForwardDiff.Dual{typeof(tag)}(θ[i], ntuple(j -> Float64(i == j), n + 1)) for i in 1:n]
    primal = InferCell._mixed_solve(s, θ)
    dual = InferCell._mixed_solve(s, θd)
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
    return (; solve_primal = med_time(() -> InferCell._mixed_solve(s, θ)),
            solve_dual = med_time(() -> InferCell._mixed_solve(s, θd)),
            stats_primal = stats(primal.sol), stats_dual = stats(dual.sol),
            tx_primal = med_time(() -> transcript_loglik(s.factors, k, μ)),
            tx_dual = med_time(() -> transcript_loglik(s.factors, kd, μd)),
            bulk_primal = med_time(() -> bulk_loglik(s, m, σ)),
            bulk_dual = med_time(() -> bulk_loglik(s, md, σd)))
end

function snapshot(out)
    s, ref = b5_setup()
    xs = points(s, ref)
    ℓ(z) = mixed_logdensity(s, z)
    vals = [ℓ(x) for x in xs]
    grads = [ForwardDiff.gradient(ℓ, x) for x in xs]
    all(isfinite, vals) || error("a snapshot point has a non-finite density: $vals")
    density = med_time(() -> ℓ(xs[1]))
    gradient = med_time(() -> ForwardDiff.gradient(ℓ, xs[1]))
    parts = breakdown(s, xs[1])
    serialize(out, (; commit = COMMIT, xs, vals, grads, density, gradient, parts,
                    threads = Threads.nthreads(), job = get(ENV, "SLURM_JOB_ID", "none")))
    @printf("snapshot at %s (job %s): density %.3f s (%.2f GiB), gradient %.3f s (%.2f GiB)\n",
            COMMIT, get(ENV, "SLURM_JOB_ID", "none"), density.seconds, density.bytes / 2^30,
            gradient.seconds, gradient.bytes / 2^30)
    for (k, v) in pairs(parts)
        println("  ", k, ": ", v isa NamedTuple && haskey(v, :seconds) ?
                @sprintf("%.4f s, %.3f GiB", v.seconds, v.bytes / 2^30) : v)
    end
end

relerr(a, b) = abs(a - b) / max(abs(a), abs(b), eps())

function compare(a_file, b_file)
    a, b = deserialize(a_file), deserialize(b_file)
    a.xs == b.xs || error("the two snapshots are at different x")
    ev = maximum(relerr.(a.vals, b.vals))
    eg = maximum(maximum(relerr.(ga, gb)) for (ga, gb) in zip(a.grads, b.grads))
    @printf("values: max rel diff %.3g; gradients: max rel diff %.3g (bound %.0e): %s\n",
            ev, eg, RTOL, ev <= RTOL && eg <= RTOL ? "PASS" : "FAIL")
    @printf("%-14s %12s %12s\n", "", a.commit, b.commit)
    @printf("%-14s %10.3f s %10.3f s\n", "density", a.density.seconds, b.density.seconds)
    @printf("%-14s %10.3f s %10.3f s\n", "gradient", a.gradient.seconds, b.gradient.seconds)
    for k in keys(a.parts)
        va, vb = a.parts[k], getproperty(b.parts, k)
        if haskey(va, :seconds)
            @printf("%-14s %10.4f s %10.4f s\n", k, va.seconds, vb.seconds)
        else
            println(rpad(k, 14), " ", va, "  →  ", vb)
        end
    end
    ev <= RTOL && eg <= RTOL || exit(1)
end

function main(args)
    if args[1] == "snapshot"
        snapshot(args[2])
    elseif args[1] == "compare"
        compare(args[2], args[3])
    else
        error("unknown stage $(args[1])")
    end
end

main(ARGS)
