# Which integrator should the ODE block use? (dev/notes/2026-10-13-inference-paths.md §3)
#
# Block 3 spends 97% of its time in the ODE step: Rodas5P, out of place over a
# 32-state `SVector`, with a dense ForwardDiff Jacobian every step and a 32×32
# `lu` that falls back to LAPACK and allocates (profile_b3_16a9b_result.md).
# syn3a-ensemble moved its full ODE to FBDF with a coloured sparse Jacobian and
# KLU, and measured 3.2× over Rodas5P (its task 2b.2). That was 146 states, its
# own tolerances, and a fresh problem at every hook. Ours is 32 states, atol 1e-10
# and rtol 1e-8, and one persistent integrator whose state is marked modified at
# every 1 s handshake, which resets a multistep method to first order. So the
# ranking has to be measured here, in the driver, and it may invert.
#
# Each variant is the same driver with a different integrator, run over the
# **same jump path**. The path is recorded once, from the production build under
# a fixed seed, and every variant replays it (`replay!`, task 16a.1), so the
# stochastic block is identical and only the ODE side differs. That makes the
# endpoints comparable: each variant's pools, in particles, are compared at every
# handshake with a tight-tolerance Rodas5P on the same path. The production
# variant's own error against that reference is the yardstick for "agrees to
# tolerance". A different integrator changes roundoff, and the rounded handshake
# can turn that into whole particles, so the report counts those too.
#
# Per variant: wall per handshake (best and median of REPS replays from one
# snapshot), bytes allocated per handshake and the GC share, and the integrator's
# own counts per handshake (accepted and rejected steps, RHS calls, Jacobians,
# factorisations), which say why a wall time is what it is.
#
# Variants, a ladder:
#   rodas5p      production: out of place over SVector, dense ForwardDiff, LAPACK lu
#   fbdf         FBDF swapped in through `ode_solver`, nothing else changed
#   rodas5p_ip   in place over a Vector, dense Jacobian, cached LU
#   fbdf_ip      the same with FBDF
#   rodas5p_klu  in place, coloured sparse Jacobian (pattern + diagonal), KLU
#   fbdf_klu     the same with FBDF: syn3a-ensemble's configuration
# The in-place variants wrap the production right-hand side, `du .= f(SVector(u))`,
# and swap the driver's integrator after the build; `src/` is not changed. The
# handshake's writes go through `_set_ode_state!`, which already takes a Vector.
#
# THREADS > 1 adds lane 6a's question: THREADS independent copies of each variant
# replay the path at once (julia -t THREADS), and the report gives the wall per
# copy against the serial best, with the GC share. That is where allocation costs.
#
# Environment:
#   MODEL      corea | m0 | both (default both). corea is build_corea() over a
#              full cycle, the 13.7 anchor (5.11 ms per handshake, job 17298242);
#              m0 is build_m0() over 600 s, as gibbs_setup builds each cell.
#   NSTEPS     handshakes per replay (default: the model's full horizon)
#   VARIANTS   subset, separated by + (sbatch --export splits on commas; default: all)
#   REPS       timed replays per variant (default 3)
#   THREADS    copies for the concurrent stage (default 1: skipped)
#   RESULT     output path (default dev/scripts/bench_ode_solver_result.md)
#
# Usage: julia --project dev/scripts/bench_ode_solver.jl
# Regenerate: sbatch dev/scripts/bench_ode_solver.slurm

using InferCell
using InferCell: HandshakeDriver
using OrdinaryDiffEq: FBDF, ODEFunction, ODEProblem, init
using StaticArrays: SVector
using SparseArrays
using LinearAlgebra: I
using UUIDs: UUID
using Random
using Printf
using Dates
using Statistics: median
import ForwardDiff

# Indirect dependencies, all in the project's manifest already: loaded by id so
# the benchmark adds nothing to Project.toml.
_load(name, uuid) = Base.require(Base.PkgId(UUID(uuid), name))
const LS = _load("LinearSolve", "7ed4a6bd-45f5-4d41-b270-4a48e9bafcae")
const ADT = _load("ADTypes", "47edcb42-4c32-4615-8424-f2b9edc5f35b")
const SMC = _load("SparseMatrixColorings", "0a514795-09f3-496d-8182-132a7b665d35")
const SCT = _load("SparseConnectivityTracer", "9f842d2f-2579-4b1d-911e-f412cf18a3f5")

const COMMIT = strip(read(`git rev-parse --short HEAD`, String))
const DIRTY = !isempty(strip(read(`git status --porcelain -- src test dev/scripts/bench_ode_solver.jl`, String)))
const JOB = get(ENV, "SLURM_JOB_ID", "none")
const SEED = 1310
const ABSTOL, RELTOL = 1e-10, 1e-8           # spec §3, pinned
# The reference only, tightest first. The driver's build sets no `maxiters`, and
# at 1e-14/1e-12 the default ran out at handshake 3345 of the full cycle (job
# 18510849), so the reference gets its own integrator with a large budget and
# steps down this ladder if that still fails. The report names the rung used.
const REF_LADDER = [(1e-14, 1e-12), (1e-13, 1e-11), (1e-12, 1e-10)]
const REF_MAXITERS = 10^8
const REPS = parse(Int, get(ENV, "REPS", "3"))
const NTHREADS = parse(Int, get(ENV, "THREADS", "1"))
const ALL_VARIANTS = [:rodas5p, :fbdf, :rodas5p_ip, :fbdf_ip, :rodas5p_klu, :fbdf_klu]
const VARIANTS = haskey(ENV, "VARIANTS") ?
    Symbol.(split(ENV["VARIANTS"], r"[,+ ]"; keepempty = false)) : ALL_VARIANTS
const RESULT = get(ENV, "RESULT", joinpath(@__DIR__, "bench_ode_solver_result.md"))
# The two pools job 17298242 saw at zero, where atol matters most.
const FLOOR_POOLS = (:M_gtp_c, :M_trna_chg_c)

const out = IOBuffer()
p(args...) = (println(out, args...); println(args...); flush(stdout))

const MODELS = Dict(
    :corea => (name = "reduced syn3A, full cycle", build = build_corea,
               models = corea_models, horizon = round(Int, COREA_CYCLE_S)),
    :m0 => (name = "M0, the sampler's build", build = build_m0,
            models = m0_models, horizon = round(Int, M0_HORIZON_S)),
)

# --- building each variant -------------------------------------------------

function seeded(build, nsteps; kw...)
    Random.seed!(SEED)
    return build(; tspan = (0.0, Float64(nsteps)), kw...)
end

# Swap a built driver's ODE integrator for `integ`, keeping every other field.
# The new integrator shares the driver's parameter vector, which the handshake
# writes the catalytic counts and geometry into.
swap_ode(d::HandshakeDriver, integ) =
    HandshakeDriver(integ, (getfield(d, f) for f in fieldnames(HandshakeDriver)[2:end])...)

# The production right-hand side, in place over a Vector.
# The length is a type parameter: captured as an Int, `SVector{n}(u)` cannot be
# inferred and every call dispatches dynamically and allocates (job 18512048).
inplace_rhs(f, n::Int) = inplace_rhs(f, Val(n))
function inplace_rhs(f, ::Val{n}) where {n}
    return function (du, u, p, t)
        du .= f(SVector{n}(u), p, t)
        return nothing
    end
end

# The Jacobian's sparsity: the union of ForwardDiff Jacobians at the sampled
# states, the tracer's pattern where it can trace the RHS, and the diagonal,
# which W = I - γJ needs and which KLU needs present in the prototype. A union
# errs towards nonzeros, so a missing entry cannot silently drop a coupling.
function jacobian_pattern(f, n, samples)
    P = falses(n, n)
    for (u, pv, t) in samples
        J = ForwardDiff.jacobian(x -> Vector(f(SVector{n}(x), pv, t)), Vector(u))
        P .|= J .!= 0
    end
    nnum = count(P)
    traced = try
        u, pv, t = first(samples)
        rhs! = inplace_rhs(f, n)
        S = ADT.jacobian_sparsity((du, x) -> rhs!(du, x, pv, t), zeros(n), Vector(u),
                                  SCT.TracerSparsityDetector())
        P .|= Matrix(S .!= 0)
        "the global tracer agreed or added $(count(P) - nnum) entries"
    catch e
        "the global tracer failed ($(typeof(e))), so the pattern is the sampled union alone"
    end
    P .|= Matrix(I, n, n)
    return sparse(Float64.(P)), nnum, traced
end

function variant_driver(v::Symbol, spec, nsteps, pattern)
    if v === :rodas5p
        return seeded(spec.build, nsteps)
    elseif v === :fbdf
        return seeded(spec.build, nsteps; ode_solver = FBDF())
    end
    d = seeded(spec.build, nsteps)
    n = length(d.ode.u)
    rhs! = inplace_rhs(d.ode.f, n)
    sparse_ad = ADT.AutoSparse(ADT.AutoForwardDiff();
                               sparsity_detector = ADT.KnownJacobianSparsityDetector(pattern),
                               coloring_algorithm = SMC.GreedyColoringAlgorithm())
    fn, alg = if v === :rodas5p_ip
        ODEFunction{true}(rhs!), Rodas5P()
    elseif v === :fbdf_ip
        ODEFunction{true}(rhs!), FBDF()
    elseif v === :rodas5p_klu
        ODEFunction{true}(rhs!; jac_prototype = copy(pattern)),
            Rodas5P(autodiff = sparse_ad, linsolve = LS.KLUFactorization())
    elseif v === :fbdf_klu
        ODEFunction{true}(rhs!; jac_prototype = copy(pattern)),
            FBDF(autodiff = sparse_ad, linsolve = LS.KLUFactorization())
    else
        error("unknown variant $v")
    end
    prob = ODEProblem(fn, plain(d.ode.u), (d.ode.t, Float64(nsteps)), d.ode.p)
    integ = init(prob, alg; abstol = ABSTOL, reltol = RELTOL, save_everystep = false)
    integ.u isa Vector{Float64} || error("the in-place state is a $(typeof(integ.u)), not a Vector")
    return swap_ode(d, integ)
end

# The production integrator, out of place, at a reference tolerance.
function reference_driver(spec, nsteps, abstol, reltol)
    d = seeded(spec.build, nsteps)
    prob = ODEProblem{false}(d.ode.f, d.ode.u, (d.ode.t, Float64(nsteps)), d.ode.p)
    return swap_ode(d, init(prob, Rodas5P(); abstol, reltol, maxiters = REF_MAXITERS,
                            save_everystep = false))
end

function reference(spec, path, nsteps)
    for (a, r) in REF_LADDER
        try
            _, X, _ = replay_recorded(reference_driver(spec, nsteps, a, r), path, nsteps)
            return X, (a, r)
        catch e
            p(@sprintf("- The reference at abstol %.0e, reltol %.0e failed: %s", a, r,
                       first(split(sprint(showerror, e), '\n'))))
        end
    end
    error("no rung of the reference ladder completed")
end

# --- measuring --------------------------------------------------------------

# A plain Vector. `collect(Float64, ::SVector)` is a `SizedVector`, because
# `similar` keeps the static axes, and so is a comprehension over one. KLU and
# sparse AD refuse it (job 18511245).
plain(u) = copyto!(Vector{Float64}(undef, length(u)), u)

# Pools in particles, as the sampler reads them (csmc.jl's latents).
particles(d) = plain(d.ode.u) .* d.factor .+ d.rounding.remainders

statvec(d) = (s = d.ode.stats; [s.naccept, s.nreject, s.nf, s.njacs, s.nw])

# Replay the path from a copy of `base`, recording the pools at every handshake.
function replay_recorded(base, path, nsteps)
    d = restore(base)
    r = PathReplay(path, d)
    X = Matrix{Float64}(undef, length(d.ode.u), nsteps)
    lo = Inf
    for k in 1:nsteps
        replay_step!(d, r)
        X[:, k] = particles(d)
        lo = min(lo, minimum(d.ode.u))
    end
    return d, X, lo
end

function timed_replay(base, path, nsteps)
    d = restore(base)
    r = PathReplay(path, d)
    s0 = statvec(d)
    t = @timed for _ in 1:nsteps
        replay_step!(d, r)
    end
    return t.time, t.bytes, t.gctime, statvec(d) .- s0
end

function concurrent_replay(base, path, nsteps, T)
    copies = [restore(base) for _ in 1:T]
    t = @timed Threads.@threads :static for i in 1:T
        replay!(copies[i], path, nsteps)
    end
    return t.time, t.gctime
end

# --- one model --------------------------------------------------------------

function run_model(key)
    spec = MODELS[key]
    nsteps = parse(Int, get(ENV, "NSTEPS", string(spec.horizon)))
    p("## ", spec.name, " (`", key, "`, ", nsteps, " handshakes)")
    p()

    # The path, recorded from the production build in the closed loop. Its
    # timing is the anchor: the warm handshakes should match 13.7's 5.11 ms on
    # the full cycle, with the SSA included.
    d = seeded(spec.build, nsteps)
    names = reduce(vcat, [states(m) for m in spec.models() if formalism(m) === :ode])
    record_path!(d)
    samples = Tuple{Vector{Float64}, Vector{Float64}, Float64}[]
    every = max(1, nsteps ÷ 20)
    twarm = @elapsed handshake_step!(d)
    trun = @elapsed for k in 2:nsteps
        handshake_step!(d)
        k % every == 0 && push!(samples, (plain(d.ode.u), copy(d.ode.p), d.ode.t))
    end
    path = recorded_path(d)
    closed = particles(d)
    p(@sprintf("Path recorded at seed %d: %d firings. The closed loop with the SSA, Rodas5P: first handshake %.1f s (compiles), then **%.3f ms per handshake**.",
               SEED, length(path), twarm, 1000 * trun / max(1, nsteps - 1)))

    n = length(d.ode.u)
    pattern, nnum, traced = jacobian_pattern(d.ode.f, n, samples)
    p(@sprintf("Jacobian pattern: %d of %d entries nonzero at the %d sampled states; with the diagonal and the tracer, %d (%s).",
               nnum, n * n, length(samples), nnz(pattern), traced))
    p()

    # The reference: the production integrator at tight tolerance, same path.
    Xref, ref_tol = reference(spec, path, nsteps)
    p(@sprintf("Reference: Rodas5P at abstol %.0e, reltol %.0e, maxiters %.0e, on the same path.",
               ref_tol..., REF_MAXITERS))
    p()
    floor_rows = filter(!isnothing, [findfirst(==(s), names) for s in FLOOR_POOLS])
    length(names) == n || error("$(length(names)) ODE names for $n states")

    rows = NamedTuple[]
    bases = Dict{Symbol, HandshakeDriver}()
    for v in VARIANTS
        base = try
            variant_driver(v, spec, nsteps, pattern)
        catch e
            p("- `", v, "`: did not build: ", sprint(showerror, e)[1:min(end, 300)])
            continue
        end
        res = try
            replay!(restore(base), path, min(5, nsteps))      # compile
            dfin, X, lo = replay_recorded(base, path, nsteps)
            times = [timed_replay(base, path, nsteps) for _ in 1:REPS]
            (; dfin, X, lo, times)
        catch e
            p("- `", v, "`: failed in the replay: ", sprint(showerror, e)[1:min(end, 300)])
            continue
        end
        bases[v] = base
        best = argmin(first.(res.times))
        wall, bytes, gct, st = res.times[best]
        E = abs.(res.X .- Xref)
        perstep = vec(maximum(E; dims = 1))
        push!(rows, (v = v, ms = 1000 * wall / nsteps,
                     spread = (median(first.(res.times)) - wall) / wall,
                     kb = bytes / nsteps / 1024, gc = gct / wall,
                     st = st ./ nsteps, err = maximum(perstep), err_end = perstep[end],
                     whole = count(>(0.5), perstep),
                     floor_err = maximum(E[floor_rows, :]), lo = res.lo,
                     clipped = clipping_census(res.dfin).clipped,
                     bitwise = v === :rodas5p ? particles(res.dfin) == closed : missing))
    end

    base_ms = only(r.ms for r in rows if r.v === :rodas5p)
    p("### Cost per handshake, replaying the path (no SSA)")
    p()
    p("| variant | ms | vs rodas5p | spread (median vs best) | KiB allocated | GC share | accepted | rejected | RHS calls | Jacobians | factorisations |")
    p("|---|---|---|---|---|---|---|---|---|---|---|")
    for r in rows
        p(@sprintf("| `%s` | %.3f | %.2f× | %+.1f%% | %.1f | %.1f%% | %.2f | %.3f | %.1f | %.2f | %.2f |",
                   r.v, r.ms, base_ms / r.ms, 100 * r.spread, r.kb, 100 * r.gc, r.st...))
    end
    p()
    p(@sprintf("### Agreement with the reference (Rodas5P at abstol %.0e, reltol %.0e, same path)", ref_tol...))
    p()
    p("Errors are in particles, the maximum over the $n pools. `whole` counts handshakes at which some pool is off by more than half a particle, which is where the rounded exchange has turned roundoff into a different count. `floor pools` is the maximum over ",
      join(("`$s`" for s in FLOOR_POOLS), " and "), ".")
    p()
    p("| variant | max error | at the end | whole | floor pools | lowest state (mM) | clipped drains |")
    p("|---|---|---|---|---|---|---|")
    for r in rows
        p(@sprintf("| `%s` | %.3g | %.3g | %d of %d | %.3g | %.3g | %d |",
                   r.v, r.err, r.err_end, r.whole, nsteps, r.floor_err, r.lo, r.clipped))
    end
    p()
    bw = only(r.bitwise for r in rows if r.v === :rodas5p)
    p("The production replay ", bw ? "reproduces" : "**does not reproduce**",
      " the recorded closed loop bitwise at the last handshake", bw ? "." : "; the comparison above is suspect until that is explained.")
    p()

    if NTHREADS > 1
        p("### $NTHREADS copies at once (julia -t $(Threads.nthreads()))")
        p()
        p("| variant | wall per copy (ms per handshake) | vs its serial best | GC share |")
        p("|---|---|---|---|")
        for r in rows
            replay!(restore(bases[r.v]), path, min(5, nsteps))     # compile the threaded path
            wall, gct = concurrent_replay(bases[r.v], path, nsteps, NTHREADS)
            # Each copy runs on its own thread, so the wall is one copy's time.
            ms = 1000 * wall / nsteps
            p(@sprintf("| `%s` | %.3f | %.2f× slower | %.1f%% |", r.v, ms, ms / r.ms, 100 * gct / wall))
        end
        p()
    end
    return rows
end

function main()
    p("# The ODE integrator in the driver: wall, allocation and agreement")
    p()
    p("Generated by `dev/scripts/bench_ode_solver.jl` at commit $COMMIT",
      DIRTY ? " (tree dirty)" : "", ", Slurm job $JOB, on $(gethostname()), ",
      "$(now()), Julia $VERSION, $(Threads.nthreads()) thread(s).")
    p()
    p("Regenerate with `sbatch dev/scripts/bench_ode_solver.slurm`. The question and ",
      "its constraints are in `dev/notes/2026-10-13-inference-paths.md` §3. Every ",
      "variant runs at abstol $ABSTOL and reltol $RELTOL.")
    p()
    which = get(ENV, "MODEL", "both")
    keys = which == "both" ? [:corea, :m0] : [Symbol(which)]
    for k in keys
        run_model(k)
    end
    write(RESULT, take!(out))
    println("Wrote ", RESULT)
end

main()
