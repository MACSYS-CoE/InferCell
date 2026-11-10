# KIM-99: where route (ii)'s log-density and gradient spend their time on the reduced
# model, on KIM-92's pilot setup, and what each candidate speed-up changes.
#
# Usage: julia --project dev/scripts/route_ii_speed_kim99.jl <stage> <out>
# Submit through dev/scripts/route_ii_speed_kim99.sh, from a worktree at a pushed commit.
# The setup is KIM-92's profile's (job 20454167): mixed_setup on d11_models() over
# COREA_TARGETS, δ = S = 0, default tolerances, on the pilot's dataset (200 cells, 52
# saves to 3120 s). x is the truth or the prior median.
# Stages:
#   split <out>     one log-density at each x, split into the route (ii) solve, the
#                   transcription rates, the transcript factor and the bulk factor,
#                   with time and GiB; then the same under a dual x of ForwardDiff's
#                   gradient chunk (11 partials, the first of the gradient's 2 passes);
#   rhs <out>       one _RouteRHS call and each module's dynamics, time and bytes, at
#                   Float64 and at dual u and p;
#   variants <out>  one solve per candidate: (a) as built; (b) FullSpecialize; (c)
#                   FullSpecialize with Rodas5P's ForwardDiff chunk set to N; (d)
#                   AutoFiniteDiff (Float64 only), each with nf, njacs, time, GiB, and the
#                   log-density's difference from (a), at Float64 and at the dual x;
#   grad <out>      the log-density and its 21-coordinate ForwardDiff gradient at each x,
#                   median of 3, as KIM-92's profile stage times them, with one solve's
#                   step counts. The before (on main) and after jobs of KIM-99.
# Each stage writes <stage>.jls and <stage>.md. SMOKE=1: a 2-cell, 120 s dataset built
# here; a plumbing pass, never a result.

include(joinpath(@__DIR__, "pilot_common.jl"))
using InferCell: _mixed_solve, _mixed_writes, _primal, route_transcription_rates,
                 mean_field_log_bulk, bulk_loglik, transcript_loglik, mean_field_parameters,
                 dynamics, _RouteRHS
using ForwardDiff
using ForwardDiff: Dual, Partials
using OrdinaryDiffEq
using SciMLBase
using Printf
using Statistics: median

const ADT = SciMLBase.ADTypes
const DATASET = "/fred/oz022/tkimpson/infercell-kim92-pilot-2026-11-10/dataset.jls"
const SEED = 17_200_001                         # KIM-92's
const NREPEAT = SMOKE ? 1 : 3
const CHUNKS = (8, 11, 12)
const NPART = 11                                # ForwardDiff's chunk for 21 coordinates

struct K99Tag end

function setup()
    ds = SMOKE ? corea_dataset(SEED; ncells = 2, horizon = 120.0, targets = COREA_TARGETS) :
         deserialize(DATASET).result
    s = mixed_setup(ds; models = d11_models(), targets = COREA_TARGETS)
    tv = Dict(ds.truth.values)
    truth = vcat([log(tv[n]) for n in s.names], log(ds.sigma))
    return s, [:truth => truth, :prior_median => mixed_prior_median(s)]
end

# x with its first NPART coordinates seeded, as ForwardDiff's first chunk seeds them.
dualx(x) = [Dual{K99Tag}(x[i], Partials(ntuple(j -> Float64(i == j), NPART))) for i in eachindex(x)]

timedrun(f) = (st = @timed f(); (; value = st.value, seconds = st.time, gib = st.bytes / 2^30,
                                  gc = st.gctime))

med(rs, k) = median(getproperty.(rs, k))

# The part of mixed_loglik after the solve, as the source states it, at a given solution.
function loglik_after(s, x, sol, p)
    n = length(s.names)
    θ, σ = exp.(x[1:n]), exp(x[n + 1])
    k = route_transcription_rates(s, sol, p)
    μ = θ[findfirst(==(:krnadeg), s.names)] ./ s.lengths
    return transcript_loglik(s.factors, k, μ) + bulk_loglik(s, mean_field_log_bulk(s, sol), σ)
end

# ---------------------------------------------------------------- split

function run_split(out)
    s, xs = setup()
    rows = []
    for (label, x0) in xs, kind in (:float, :dual)
        x = kind === :float ? x0 : dualx(x0)
        n = length(s.names)
        θ, σ = exp.(x[1:n]), exp(x[n + 1])
        ld = mixed_logdensity(s, x)                       # warm
        reps = map(1:NREPEAT) do _
            t_solve = timedrun(() -> _mixed_solve(s, θ))
            r = t_solve.value
            t_rates = timedrun(() -> route_transcription_rates(s, r.sol, r.p))
            k = t_rates.value
            μ = θ[findfirst(==(:krnadeg), s.names)] ./ s.lengths
            t_tx = timedrun(() -> transcript_loglik(s.factors, k, μ))
            t_bulk = timedrun(() -> bulk_loglik(s, mean_field_log_bulk(s, r.sol), σ))
            t_all = timedrun(() -> mixed_logdensity(s, x))
            (; solve = t_solve, rates = t_rates, tx = t_tx, bulk = t_bulk, all = t_all,
             stats = r.sol.stats)
        end
        part(k) = (; seconds = median([getproperty(r, k).seconds for r in reps]),
                   gib = median([getproperty(r, k).gib for r in reps]))
        st = reps[1].stats
        row = (; label, kind, logdensity = _primal(ld), solve = part(:solve), rates = part(:rates),
               tx = part(:tx), bulk = part(:bulk), all = part(:all), nf = st.nf,
               njacs = st.njacs, naccept = st.naccept, nreject = st.nreject)
        say(@sprintf("split %s %s: all %.1f s; solve %.1f s, rates %.2f s, tx %.1f s, bulk %.2f s; nf %d, njacs %d",
                     label, kind, row.all.seconds, row.solve.seconds, row.rates.seconds,
                     row.tx.seconds, row.bulk.seconds, row.nf, row.njacs))
        push!(rows, row)
    end
    save(out, "split", rows; stage = "split")
    md(out, "split") do p
        p("| x | number | all s | solve s (GiB) | rates s | transcript s (GiB) | bulk s | nf | njacs | steps |")
        p("|---|---|---|---|---|---|---|---|---|---|")
        for r in rows
            p(@sprintf("| %s | %s | %.1f | %.1f (%.1f) | %.2f | %.1f (%.1f) | %.2f | %d | %d | %d + %d |",
                       r.label, r.kind, r.all.seconds, r.solve.seconds, r.solve.gib,
                       r.rates.seconds, r.tx.seconds, r.tx.gib, r.bulk.seconds, r.nf, r.njacs,
                       r.naccept, r.nreject))
        end
    end
end

# ---------------------------------------------------------------- rhs

# Seconds per call and bytes per call of `f()`, after one warm call.
function percall(f; n = SMOKE ? 50 : 2000)
    f()
    b = @allocated f()
    t = @elapsed for _ in 1:n
        f()
    end
    return (; us = 1e6 * t / n, bytes = b)
end

function run_rhs(out)
    s, xs = setup()
    rt = s.route
    f = rt.prob.f.f
    f isa _RouteRHS || error("prob.f.f is a $(typeof(f)), not the route's RHS")
    x = Dict(xs)[:truth]
    p = mean_field_parameters(rt, _mixed_writes(s, exp.(x[1:length(s.names)])))
    rows = []
    for kind in (:float, :dual)
        T = kind === :float ? Float64 : Dual{K99Tag, Float64, NPART}
        u = T.(rt.prob.u0)
        q = kind === :float ? p : [Dual{K99Tag}(v, Partials(ntuple(j -> Float64(j == mod1(i, NPART)), NPART)))
                                   for (i, v) in enumerate(p)]
        du = similar(u, promote_type(eltype(u), eltype(q)))
        r = percall(() -> f(du, u, q, 0.0))
        push!(rows, (; kind, what = "_RouteRHS", r...))
        geo = InferCell._route_geometry(rt, rt.n_idx == 0 ? zero(eltype(u)) : u[rt.n_idx])
        for (m, sl, sb) in zip(f.models, f.slots, f.subs)
            uL, ui = u[sl.sidx], u[sl.iidx]
            pL = InferCell._substitute(q[sl.pidx], sb, u, geo)
            r = percall(() -> dynamics(uL, pL, 0.0, m, ui))
            push!(rows, (; kind, what = string(nameof(typeof(m)), " ", InferCell.module_id(m)), r...))
        end
    end
    for r in rows
        say(@sprintf("rhs %s %s: %.2f µs, %d bytes", r.kind, r.what, r.us, r.bytes))
    end
    save(out, "rhs", rows; stage = "rhs")
    md(out, "rhs") do p
        p("| number | call | µs per call | bytes per call |")
        p("|---|---|---|---|")
        for r in rows
            p(@sprintf("| %s | %s | %.2f | %d |", r.kind, r.what, r.us, r.bytes))
        end
    end
end

# ---------------------------------------------------------------- variants

# solve_mean_field's solve, with the problem's specialisation and the solver chosen.
function solve_variant(s, θ; full::Bool, solver)
    rt = s.route
    q = copy(mean_field_parameters(rt, _mixed_writes(s, θ)))
    u0 = convert(Vector{promote_type(eltype(q), eltype(rt.prob.u0))}, rt.prob.u0)
    prob = full ?
        ODEProblem{true, SciMLBase.FullSpecialize}(rt.prob.f.f, u0, rt.prob.tspan, q; rt.prob.kwargs...) :
        remake(rt.prob; u0, p = q)
    saveat = sort(unique(vcat(s.switches, s.saves)))
    sol = solve(prob, solver; saveat, abstol = s.abstol, reltol = s.reltol, maxiters = s.maxiters,
                save_everystep = false)
    return (; sol, p = q)
end

function run_variants(out)
    s, xs = setup()
    variants = vcat([(:a_as_built, false, Rodas5P()), (:b_full_specialize, true, Rodas5P())],
                    [(Symbol(:c_full_chunk_, N), true, Rodas5P(autodiff = ADT.AutoForwardDiff(chunksize = N)))
                     for N in CHUNKS],
                    [(:d_finite_diff, false, Rodas5P(autodiff = ADT.AutoFiniteDiff()))])
    rows = []
    for (label, x0) in xs, kind in (:float, :dual)
        x = kind === :float ? x0 : dualx(x0)
        θ = exp.(x[1:length(s.names)])
        ref = nothing
        for (name, full, solver) in variants
            kind === :dual && name === :d_finite_diff && continue
            solve_variant(s, θ; full, solver)            # warm
            rs = [timedrun(() -> solve_variant(s, θ; full, solver)) for _ in 1:NREPEAT]
            r = rs[1].value
            ℓ = loglik_after(s, x, r.sol, r.p)
            ref === nothing && (ref = ℓ)
            dv = abs(_primal(ℓ) - _primal(ref))
            dp = kind === :dual ? maximum(abs.(ForwardDiff.partials(ℓ) .- ForwardDiff.partials(ref))) : 0.0
            st = r.sol.stats
            row = (; label, kind, variant = name, seconds = med(rs, :seconds), gib = med(rs, :gib),
                   nf = st.nf, njacs = st.njacs, naccept = st.naccept, nreject = st.nreject,
                   retcode = r.sol.retcode, loglik = _primal(ℓ), dvalue = dv, dpartials = dp)
            say(@sprintf("variant %s %s %s: %.2f s, %.1f GiB, nf %d, njacs %d, |Δℓ| %.3g, |Δ∂ℓ| %.3g, %s",
                         label, kind, name, row.seconds, row.gib, row.nf, row.njacs, dv, dp, row.retcode))
            push!(rows, row)
        end
    end
    save(out, "variants", rows; stage = "variants")
    md(out, "variants") do p
        p("Δ is against (a), the problem as built, at the same x and number type.\n")
        p("| x | number | variant | solve s | GiB | nf | njacs | steps | log-lik | abs Δℓ | max abs Δ∂ℓ |")
        p("|---|---|---|---|---|---|---|---|---|---|---|")
        for r in rows
            p(@sprintf("| %s | %s | %s | %.2f | %.1f | %d | %d | %d + %d | %.6f | %.3g | %.3g |",
                       r.label, r.kind, r.variant, r.seconds, r.gib, r.nf, r.njacs, r.naccept,
                       r.nreject, r.loglik, r.dvalue, r.dpartials))
        end
    end
end

# ---------------------------------------------------------------- grad

function run_grad(out)
    s, xs = setup()
    ℓ(z) = mixed_logdensity(s, z)
    t_compile = @elapsed ForwardDiff.gradient(ℓ, xs[1][2])
    rows = []
    for (label, x) in xs
        d = [timedrun(() -> ℓ(x)) for _ in 1:NREPEAT]
        g = [timedrun(() -> ForwardDiff.gradient(ℓ, x)) for _ in 1:NREPEAT]
        st = _mixed_solve(s, exp.(x[1:length(s.names)])).sol.stats
        row = (; label, logdensity = d[1].value, grad = g[1].value,
               density = [Base.structdiff(r, (; value = 0)) for r in d],
               gradient = [Base.structdiff(r, (; value = 0)) for r in g],
               nf = st.nf, njacs = st.njacs, naccept = st.naccept, nreject = st.nreject)
        say(@sprintf("grad %s: density %s s (ℓ %.6f); gradient %s s (|g| %.4g); nf %d, njacs %d",
                     label, join((@sprintf("%.1f", r.seconds) for r in d), ", "), row.logdensity,
                     join((@sprintf("%.1f", r.seconds) for r in g), ", "), sqrt(sum(abs2, row.grad)),
                     row.nf, row.njacs))
        push!(rows, row)
    end
    save(out, "grad", (; rows, compile_seconds = t_compile); stage = "grad")
    md(out, "grad") do p
        p(@sprintf("Compiled in %.0f s (the first gradient at the truth, untimed below).\n", t_compile))
        p("| x | call | seconds (each) | median s | GiB (median) | log density or abs g |")
        p("|---|---|---|---|---|---|")
        for r in rows
            p(@sprintf("| %s | log density | %s | %.1f | %.1f | %.6f |", r.label,
                       join((@sprintf("%.1f", q.seconds) for q in r.density), ", "),
                       med(r.density, :seconds), med(r.density, :gib), r.logdensity))
            p(@sprintf("| %s | gradient | %s | %.1f | %.1f | %.4g |", r.label,
                       join((@sprintf("%.1f", q.seconds) for q in r.gradient), ", "),
                       med(r.gradient, :seconds), med(r.gradient, :gib), sqrt(sum(abs2, r.grad))))
        end
        p("\n| x | nf | njacs | accepted | rejected |\n|---|---|---|---|---|")
        for r in rows
            p(@sprintf("| %s | %d | %d | %d | %d |", r.label, r.nf, r.njacs, r.naccept, r.nreject))
        end
    end
end

# ---------------------------------------------------------------- output

function md(body, out, stage)
    io = IOBuffer()
    p(xs...) = println(io, xs...)
    p("# KIM-99 `", stage, "`\n")
    p("Commit `", COMMIT, "`, job ", JOB, ", on ", RAN_ON.host, " (", RAN_ON.cpu, "), ",
      Threads.nthreads(), " thread(s). SMOKE = ", SMOKE, ".\n")
    body(p)
    path = joinpath(out, "$stage.md")
    write(path, String(take!(io)))
    print(read(path, String))
end

function main(args)
    length(args) == 2 || error("usage: route_ii_speed_kim99.jl <stage> <out>")
    stage, out = args
    mkpath(out)
    say("KIM-99 ", stage, " at ", COMMIT, ", job ", JOB, ", on ", RAN_ON.host, " (", RAN_ON.cpu, ")")
    stage == "split" ? run_split(out) :
    stage == "rhs" ? run_rhs(out) :
    stage == "variants" ? run_variants(out) :
    stage == "grad" ? run_grad(out) : error("unknown stage $stage")
end

main(ARGS)
