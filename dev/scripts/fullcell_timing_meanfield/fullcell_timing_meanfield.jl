# KIM-184: what one gradient of the full syn3A cell's first-order mean field costs, solved
# by a sparse stiff solver under ForwardDiff dual numbers.
#
# Usage: julia --project=dev/scripts/fullcell_timing_meanfield \
#            dev/scripts/fullcell_timing_meanfield/fullcell_timing_meanfield.jl <stage> <out>
# Submit through fullcell_timing_meanfield.sh, from a worktree at a pushed commit.
#
# The model. The CME table of the full cell (syn3a-ensemble data/model/cme, the minute-1
# rate constants of `reactions.csv`), first-order reactions only (LM type 1, propensity
# k·x). The mean field of a first-order network is linear, dx/dt = S·(k ∘ x[r]), so the
# Jacobian is S·diag(k)·E, written analytically into a sparse prototype with a full
# structural diagonal. The state is every species a first-order reaction reads or
# changes, in the table's order. The 363 second-order reactions, the 175 metabolic ODE
# fluxes and the minute-by-minute rate rebuilds are left out, so the timing is a lower
# bound on the real system's cost and says nothing about accuracy.
#
# The gradient. θ_b = log(k_b / k_b⁰) over a block of B reactions, the B largest by
# flux integrated over the fit window at θ = 0 (ties by reaction index), so the blocks
# of 10, 20 and 40 are nested. The scalar is a Gaussian misfit over every state on the
# reduced model's save grid (60 s to 3120 s, solved to COREA_FIT_S = 3150 s) against
# the solution at θ = 0, with variance max(d, 1), evaluated at θ*_b = 0.05·(−1)^b. The
# gradient is forward mode in passes of CHUNK partials each, as ForwardDiff.gradient
# takes them; a pass is one dual solve. Threaded, the passes run one per thread.
#
# Stages (each writes <out>/<stage>.md and <stage>.jls):
#   check     state, reaction and non-zero counts; the RHS at the initial state against
#             an independent BigFloat sum of the table's propensities (1e-10 relative);
#             J·x0 against the RHS; the ranked block; the dual gradient against central
#             differences on 3 parameters at tight tolerances (1e-6 relative, as AV3);
#             the pass-wise gradient against ForwardDiff.gradient (bitwise; both go
#             through DualSparseKLU, so this checks pass assembly only); which linear
#             cache the solver builds under a dual state (asserted); and DualSparseKLU
#             against LinearSolve's own dual path (DualLinearCache) on one pass of CHUNK
#             partials over LIB_T s (bitwise, asserted)
#   time      serial: the primal solve, and the gradient at B = 10, 20 and 40, each the
#             median of NREPEAT after a warm-up, with per-pass times and solver stats
#   threaded  the same, the passes on threads (needs ≥ 4 threads); also the threaded
#             gradient against the serial one (bitwise, asserted when time.jls is
#             in the outdir)
#   library   one dual pass (B = 20, the first 10 partials) over the whole window through
#             LinearSolve's own dual path (DualLinearCache), against DualSparseKLU's,
#             each after a warm-up on a 30 s window; one repeat of the library pass;
#             the two asserted bitwise equal
# Every solve's retcode is checked; a failed guard stops the stage with an error. The
# guards were added after the runs of record (74755da, 6a20e2c) and change no number.
# SMOKE=1: the window cut to 300 s and 2 repeats; a plumbing pass, never a result.

using OrdinaryDiffEq
using SciMLBase
using LinearSolve
using SparseArrays
using LinearAlgebra
using ForwardDiff
using ForwardDiff: Dual, Partials, value, partials
using Printf
using Serialization
using Statistics: median
using Dates
using InteractiveUtils: versioninfo

const SMOKE = get(ENV, "SMOKE", "0") == "1"
const MODEL = "/fred/oz022/tkimpson/syn3a-ensemble/data/model"
const SPECIES_CSV = joinpath(MODEL, "cme", "species.csv")
const REACTIONS_CSV = joinpath(MODEL, "cme", "reactions.csv")
# The reduced model's fit window and save grid (InferCell's COREA_FIT_S, and the 60 s
# grid KIM-92's dataset saves on; corea_dataset keeps saves at 0 < t ≤ fit_to).
const COREA_FIT_S = 3150.0
const T_END = SMOKE ? 300.0 : COREA_FIT_S
const SAVES = collect(60.0:60.0:(T_END - mod(T_END, 60.0)))
# Route (ii)'s pinned tolerances (InferCell solve_mean_field), and its solver family.
const ABSTOL = 1e-10
const RELTOL = 1e-8
const MAXITERS = 10^5
# Tight tolerances for the finite-difference check, as AV3's test uses.
const ABSTOL_TIGHT = 1e-14
const RELTOL_TIGHT = 1e-12
const FD_H = 1e-4
const FD_RTOL = 1e-6
const RHS_RTOL = 1e-10
const CHUNK = 10                       # ForwardDiff.pickchunksize(B) at B = 10, 20 and 40
const BLOCKS = (10, 20, 40)
const NREPEAT = SMOKE ? 2 : 5
const LIB_T = 30.0                     # the window of the check against LinearSolve's dual path

# ---------------------------------------------------------------------------------------
# The linear solve under duals
#
# A dual solve gives Rodas5P a W = I/γ − J with dual entries. LinearSolve 3.73.0's
# ForwardDiff extension would take it: it solves the value part, then each partial
# against the same factorisation, which is right. But it has only the generic dense
# `partials_to_list` / `update_partials_list!`: each partial of W becomes a dense n×n
# matrix, filled by sparse getindex over all n² entries, and is multiplied densely. One
# `rhs_cache_valid` flag is cleared by both `setA!` and `setb!`, so the dense rebuild
# (`update_partials_list!`) runs at every stage solve, 8 per Rodas5P step, not once
# per W. On this system that cost 3.85 s per step over the whole window (`library`
# stage, job 21667143: 554.1 s for one pass of 10 partials over 144 steps, against
# 0.645 s here). LinearSolve 3.87 has sparse methods for both and separate A/b flags;
# it is untimed here. DualSparseKLU does
# the same algebra on the sparse pattern: KLU on the value part, refactored (symbolic
# reused) only when the values change, the primal solved, then all the partials at once
# as an n×N right-hand side, R = ∂b − ∂W·x₀, formed from the non-zeros. The primal
# solve keeps LinearSolve's KLUFactorization.

# Three LinearSolve internals: the KLU module vendored in its SparseArrays extension,
# the ForwardDiff extension's DualAbstractLinearProblem, and LinearSolve.__init (below).
const KLU = Base.get_extension(LinearSolve, :LinearSolveSparseArraysExt).KLU
const LSFD = Base.get_extension(LinearSolve, :LinearSolveForwardDiffExt)

# Holds a factorisation and buffers: one instance per solve (`solver` makes a fresh one
# each call). Not thread-safe to share between solves.
mutable struct DualSparseKLU
    F::Any                       # the KLU factorisation of value(W), or nothing
    vals::Vector{Float64}        # the values it factored
    x0::Vector{Float64}
    R::Matrix{Float64}
    nfact::Int                   # numeric factorisations
    nsolve::Int                  # solves
end
DualSparseKLU() = DualSparseKLU(nothing, Float64[], Float64[], zeros(0, 0), 0, 0)

function (s::DualSparseKLU)(A, b, u, p, isfresh, Pl, Pr, cacheval; kwargs...)
    T = eltype(A)
    N = ForwardDiff.npartials(T)
    nz = nonzeros(A); rv = rowvals(A); n = size(A, 1)
    if s.F === nothing
        s.vals = value.(nz)
        s.F = KLU.klu(SparseMatrixCSC(n, n, copy(SparseArrays.getcolptr(A)), copy(rv), copy(s.vals)))
        s.x0 = zeros(n); s.R = zeros(n, N); s.nfact += 1
    else
        length(nz) == length(s.vals) || error("W's pattern changed")
        changed = false
        @inbounds for k in eachindex(nz)
            v = value(nz[k])
            if v != s.vals[k]
                s.vals[k] = v; changed = true
            end
        end
        changed && (KLU.klu!(s.F, s.vals); s.nfact += 1)
    end
    x0 = s.x0; R = s.R
    @inbounds for r in 1:n
        x0[r] = value(b[r])
        pb = partials(b[r])
        for i in 1:N
            R[r, i] = pb[i]
        end
    end
    ldiv!(s.F, x0)
    @inbounds for c in 1:n
        xc = x0[c]
        xc == 0 && continue
        for k in nzrange(A, c)
            r = rv[k]; pk = partials(nz[k])
            for i in 1:N
                R[r, i] -= pk[i] * xc
            end
        end
    end
    ldiv!(s.F, R)
    @inbounds for r in 1:n
        u[r] = T(x0[r], Partials(ntuple(i -> R[r, i], Val(N))))
    end
    s.nsolve += 1
    return u
end

const DualKLUAlg = LinearSolveFunction{DualSparseKLU}
# Keep the dual problem away from DualLinearCache: a plain LinearCache, as LinearSolve
# itself does for SparspakFactorization.
SciMLBase.init(prob::LSFD.DualAbstractLinearProblem, alg::DualKLUAlg, args...; kwargs...) =
    LinearSolve.__init(prob, alg, args...; kwargs...)

# A concrete W, not OrdinaryDiffEq's lazy WOperator, which it builds for a solver that
# does not say it needs one.
LinearSolve.needs_concrete_A(::DualKLUAlg) = true

solver(::Type{<:Dual}) = Rodas5P(linsolve = LinearSolveFunction(DualSparseKLU()))
solver(::Type) = Rodas5P(linsolve = KLUFactorization())
# LinearSolve's own dual path, for the check against DualSparseKLU.
solver_library() = Rodas5P(linsolve = KLUFactorization())

struct K184Tag end
const TAG = typeof(ForwardDiff.Tag(K184Tag(), Float64))

say(s) = (println(Dates.format(now(), "HH:MM:SS"), "  ", s); flush(stdout))

# ---------------------------------------------------------------------------------------
# The table

"""The CME table's rows as strings, header dropped. No field holds a comma."""
csvrows(path) = [split(l, ',') for l in eachline(path)][2:end]

struct FirstOrderModel
    n::Int                       # states
    names::Vector{String}        # state names
    cme_index::Vector{Int}       # state → CME species index
    x0::Vector{Float64}          # initial counts (procid 1)
    rxn::Vector{Int}             # first-order reaction → its row index in the table
    k0::Vector{Float64}          # minute-1 rate constants
    reactant::Vector{Int}        # reaction → its reactant's state
    ptr::Vector{Int}             # reaction j's stoichiometry is ptr[j]:ptr[j+1]-1
    row::Vector{Int}             # … the state it changes
    coef::Vector{Float64}        # … by this much
    jac_proto::SparseMatrixCSC{Float64, Int}
    nzpos::Vector{Int}           # stoichiometry entry → its slot in jac_proto.nzval
    nnz_table::Int               # distinct (row, reactant) pairs, before the diagonal
end

function load_model()
    sp = csvrows(SPECIES_CSV)
    cme_names = String[r[2] for r in sp]
    cme_x0 = Float64[parse(Float64, r[3]) for r in sp]
    @assert all(parse(Int, r[1]) == i for (i, r) in enumerate(sp))
    rows = csvrows(REACTIONS_CSV)
    first = [r for r in rows if r[2] == "1"]
    # The state: every species a first-order reaction reads or changes, in CME order.
    touched = Set{Int}()
    for r in first
        push!(touched, parse(Int, r[4]))
        for pr in split(r[5], ';')
            push!(touched, parse(Int, split(pr, ':')[1]))
        end
    end
    cme_index = sort!(collect(touched))
    state = Dict(c => i for (i, c) in enumerate(cme_index))
    n = length(cme_index)
    rxn = Int[]; k0 = Float64[]; reactant = Int[]; ptr = Int[1]; row = Int[]; coef = Float64[]
    for r in first
        push!(rxn, parse(Int, r[1]))
        push!(k0, parse(Float64, r[3]))
        push!(reactant, state[parse(Int, r[4])])
        for pr in split(r[5], ';')
            s, d = split(pr, ':')
            push!(row, state[parse(Int, s)])
            push!(coef, parse(Float64, d))
        end
        push!(ptr, length(row) + 1)
    end
    # The Jacobian's pattern: (row, reactant) for every stoichiometry entry, plus the
    # diagonal so W = I/γ − J has a structurally complete one.
    cols = Int[]
    for j in eachindex(rxn), e in ptr[j]:(ptr[j + 1] - 1)
        push!(cols, reactant[j])
    end
    nnz_table = length(unique(zip(row, cols)))
    I = vcat(row, 1:n); J = vcat(cols, 1:n)
    proto = sparse(I, J, zeros(length(I)), n, n)
    nzpos = map(eachindex(row)) do e
        c = cols[e]
        rng = nzrange(proto, c)
        rng[findfirst(==(row[e]), view(rowvals(proto), rng))]
    end
    return FirstOrderModel(n, cme_names[cme_index], cme_index, cme_x0[cme_index], rxn, k0,
                           reactant, ptr, row, coef, proto, nzpos, nnz_table)
end

# ---------------------------------------------------------------------------------------
# The right-hand side, its Jacobian, and the objective

"""θ's parameters: `slot[j]` is reaction j's position in θ, or 0 if it is fixed."""
struct RHS{M}
    m::M
    slot::Vector{Int}
end

@inline rate(f::RHS, j, θ) =
    (s = f.slot[j]; s > 0 ? f.m.k0[j] * exp(θ[s]) : f.m.k0[j] * one(eltype(θ)))

function (f::RHS)(du, u, θ, t)
    m = f.m
    fill!(du, zero(eltype(du)))
    @inbounds for j in eachindex(m.k0)
        a = rate(f, j, θ) * u[m.reactant[j]]
        for e in m.ptr[j]:(m.ptr[j + 1] - 1)
            du[m.row[e]] += m.coef[e] * a
        end
    end
    return nothing
end

struct Jac{M}
    f::RHS{M}
end

function (jf::Jac)(J, u, θ, t)
    m = jf.f.m
    nz = nonzeros(J)
    fill!(nz, zero(eltype(nz)))
    @inbounds for j in eachindex(m.k0)
        kj = rate(jf.f, j, θ)
        for e in m.ptr[j]:(m.ptr[j + 1] - 1)
            nz[m.nzpos[e]] += m.coef[e] * kj
        end
    end
    return nothing
end

tgrad!(dT, u, p, t) = (fill!(dT, zero(eltype(dT))); nothing)

function problem(m, slot, θ; tend = T_END)
    f = RHS(m, slot)
    fun = ODEFunction{true, SciMLBase.FullSpecialize}(f; jac = Jac(f),
                                                     jac_prototype = copy(m.jac_proto),
                                                     tgrad = tgrad!)
    u0 = convert(Vector{eltype(θ)}, m.x0)
    return ODEProblem(fun, u0, (0.0, tend), θ)
end

msolve(prob; abstol = ABSTOL, reltol = RELTOL, saveat = SAVES, alg = solver(eltype(prob.p))) =
    solve(prob, alg; abstol, reltol, saveat, maxiters = MAXITERS, save_everystep = false,
          save_start = false)

"""The data: the primal solution at θ = 0 on the save grid, and its variances."""
function make_data(m, slot, B; kw...)
    sol = msolve(problem(m, slot, zeros(B)); kw...)
    SciMLBase.successful_retcode(sol) || error("data solve: $(sol.retcode)")
    d = reduce(hcat, sol.u)
    return (d = d, w = 1 ./ max.(d, 1.0))
end

"""The misfit and the solve's stats, at θ (Float64 or dual)."""
function objective(m, slot, data, θ; kw...)
    sol = msolve(problem(m, slot, θ); kw...)
    SciMLBase.successful_retcode(sol) || error("solve: $(sol.retcode)")
    T = eltype(θ)
    ℓ = zero(T)
    for (i, u) in enumerate(sol.u)
        @inbounds for s in eachindex(u)
            r = u[s] - data.d[s, i]
            ℓ -= r * r * data.w[s, i] / 2
        end
    end
    return ℓ, sol.stats
end

statsnt(s) = (nf = s.nf, njacs = s.njacs, nw = s.nw, nsolve = s.nsolve,
              naccept = s.naccept, nreject = s.nreject)

# ---------------------------------------------------------------------------------------
# The gradient, pass by pass

passes(B) = [i:min(i + CHUNK - 1, B) for i in 1:CHUNK:B]

function seed(θ, idx)
    return [Dual{TAG}(θ[i], Partials(ntuple(c -> (c <= length(idx) && idx[c] == i) ? 1.0 : 0.0,
                                            Val(CHUNK)))) for i in eachindex(θ)]
end

"""One dual pass: the value, its partials over `idx`, the seconds and the solve's stats."""
function dual_pass(m, slot, data, θ, idx; kw...)
    t = @elapsed begin
        v, st = objective(m, slot, data, seed(θ, idx); kw...)
    end
    return (value = value(v), partials = collect(partials(v))[1:length(idx)], seconds = t,
            stats = statsnt(st))
end

"""The gradient, its passes serial or one per thread; the wall time and each pass."""
function gradient_passes(m, slot, data, θ; threaded = false, kw...)
    ps = passes(length(θ))
    out = Vector{Any}(undef, length(ps))
    wall = @elapsed begin
        if threaded
            Threads.@threads :static for c in eachindex(ps)
                out[c] = dual_pass(m, slot, data, θ, ps[c]; kw...)
            end
        else
            for c in eachindex(ps)
                out[c] = dual_pass(m, slot, data, θ, ps[c]; kw...)
            end
        end
    end
    g = zeros(length(θ))
    for (c, idx) in enumerate(ps)
        g[idx] .= out[c].partials
    end
    return (grad = g, wall = wall, passes = out)
end

# ---------------------------------------------------------------------------------------
# The block

"""Every first-order reaction's flux integrated over the window at θ = 0, and the
reactions ranked by it, ties by reaction index. The integral is 3-point Gauss–Legendre
on each of the solver's own steps, through its dense interpolant, so a fast drain in
the first step is integrated, not smeared over a coarse grid (the SMOKE run's 10 s
trapezoid ranked 22 tRNA-charging drains of one molecule first at 150 each)."""
function rank_reactions(m)
    sol = solve(problem(m, zeros(Int, length(m.k0)), Float64[]), solver(Float64);
                abstol = ABSTOL, reltol = RELTOL, maxiters = MAXITERS, dense = true)
    SciMLBase.successful_retcode(sol) || error("ranking solve: $(sol.retcode)")
    nodes = (-sqrt(3 / 5), 0.0, sqrt(3 / 5)); weights = (5 / 9, 8 / 9, 5 / 9)
    xint = zeros(m.n)
    for i in 1:(length(sol.t) - 1)
        a, b = sol.t[i], sol.t[i + 1]
        h = (b - a) / 2
        for (z, w) in zip(nodes, weights)
            xint .+= (w * h) .* sol(a + h * (1 + z))
        end
    end
    flux = [m.k0[j] * xint[m.reactant[j]] for j in eachindex(m.k0)]
    order = sortperm(collect(zip(-flux, m.rxn)))
    return flux, order
end

function block_slot(m, order, B)
    slot = zeros(Int, length(m.k0))
    for (b, j) in enumerate(order[1:B])
        slot[j] = b
    end
    return slot
end

θstar(B) = [0.05 * (-1)^b for b in 1:B]

# ---------------------------------------------------------------------------------------
# Environment

function environment()
    io = IOBuffer()
    versioninfo(io)
    cpu = try
        strip(split(first(filter(l -> startswith(l, "Model name"), readlines(`lscpu`))), ':')[2])
    catch
        "unknown"
    end
    return (host = gethostname(), cpu = cpu, julia = string(VERSION),
            threads = Threads.nthreads(), blas = BLAS.get_num_threads(),
            gcsweep = Base.JLOptions().nsweepthreads, ngc = Threads.ngcthreads(),
            job = get(ENV, "SLURM_JOB_ID", "none"), smoke = SMOKE, t_end = T_END,
            nsaves = length(SAVES), versioninfo = String(take!(io)))
end

function header(io, stage, env)
    println(io, "# KIM-184 `$stage`", SMOKE ? " (SMOKE: a plumbing pass, not a result)" : "")
    println(io)
    println(io, "- job $(env.job) on $(env.host), $(env.cpu); Julia $(env.julia)")
    println(io, "- threads $(env.threads), BLAS threads $(env.blas), GC threads $(env.ngc) (sweep $(env.gcsweep))")
    println(io, "- window 0–$(env.t_end) s, $(env.nsaves) saves; Rodas5P, KLUFactorization " *
                "(primal) and DualSparseKLU (dual, unless a row says DualLinearCache), " *
                "abstol $ABSTOL, reltol $RELTOL; chunk $CHUNK")
    println(io)
end

# ---------------------------------------------------------------------------------------
# Stages

function stage_check(out)
    env = environment()
    say("check on $(env.host) ($(env.cpu)), job $(env.job)")
    m = load_model()
    rows = csvrows(REACTIONS_CSV)
    n_first = count(r -> r[2] == "1", rows)
    n_second = count(r -> r[2] == "2", rows)
    n_other = length(rows) - n_first - n_second
    nsp = length(csvrows(SPECIES_CSV))
    say("states $(m.n) of $nsp species; first-order $n_first, second-order $n_second, other $n_other")
    say("Jacobian non-zeros: $(m.nnz_table) from the table, $(nnz(m.jac_proto)) with the diagonal")
    @assert n_first == length(m.k0)

    # The RHS at the initial state against the table's propensities, summed
    # independently: straight from the rows, by species name, in BigFloat.
    setprecision(BigFloat, 256)
    sp = csvrows(SPECIES_CSV)
    x0name = Dict(String(r[2]) => parse(BigFloat, r[3]) for r in sp)
    cmename = Dict(parse(Int, r[1]) => String(r[2]) for r in sp)
    ref = Dict{String, BigFloat}()
    absum = Dict{String, BigFloat}()
    for r in rows
        r[2] == "1" || continue
        a = parse(BigFloat, r[3]) * x0name[cmename[parse(Int, r[4])]]
        for pr in split(r[5], ';')
            s, d = split(pr, ':')
            nm = cmename[parse(Int, s)]
            term = parse(BigFloat, d) * a
            ref[nm] = get(ref, nm, big(0)) + term
            absum[nm] = get(absum, nm, big(0)) + abs(term)
        end
    end
    f = RHS(m, zeros(Int, length(m.k0)))
    du = zeros(m.n)
    f(du, m.x0, Float64[], 0.0)
    nonzero_ref = 0; zero_ref = 0; zero_mismatch = 0; maxrel = 0.0; maxrel_abs = 0.0
    for (i, nm) in enumerate(m.names)
        rv = get(ref, nm, big(0))
        if iszero(rv)
            zero_ref += 1
            du[i] == 0 || (zero_mismatch += 1)
        else
            nonzero_ref += 1
            maxrel = max(maxrel, Float64(abs(du[i] - rv) / abs(rv)))
            maxrel_abs = max(maxrel_abs, Float64(abs(du[i] - rv) / absum[nm]))
        end
    end
    say(@sprintf("RHS vs BigFloat propensity sums: %d non-zero, max rel %.3e (vs Σ|terms| %.3e); %d exact zeros, %d not zero in the RHS",
                 nonzero_ref, maxrel, maxrel_abs, zero_ref, zero_mismatch))
    J = copy(m.jac_proto)
    Jac(f)(J, m.x0, Float64[], 0.0)
    jx = J * m.x0
    jrel = maximum(abs(jx[i] - du[i]) / max(abs(du[i]), eps()) for i in eachindex(du) if du[i] != 0)
    jzero = all(jx[i] == 0 for i in eachindex(du) if du[i] == 0)
    say(@sprintf("J·x0 vs RHS: max rel %.3e; zeros kept: %s", jrel, jzero))
    rhs_pass = n_first == 2002 && maxrel <= RHS_RTOL && zero_mismatch == 0

    # The block.
    flux, order = rank_reactions(m)
    say("ranked: $(count(>(0), flux)) of $(length(flux)) reactions carry flux over the window")
    ranked = [(rank = b, reaction = m.rxn[j], reactant = m.names[m.reactant[j]], k0 = m.k0[j],
               flux = flux[j]) for (b, j) in enumerate(order[1:maximum(BLOCKS)])]

    # How the dual linear solve is done: the solver's linear cache under a dual u.
    B = 20
    slot = block_slot(m, order, B)
    linpath = try
        pd = problem(m, slot, seed(θstar(B), 1:CHUNK))
        integ = init(pd, solver(eltype(pd.p)); abstol = ABSTOL, reltol = RELTOL, saveat = SAVES)
        lc = integ.cache.linsolve
        string(nameof(typeof(lc)), " / ", nameof(typeof(lc.alg)), " / W ",
               typeof(integ.cache.W).name.name, "{", eltype(integ.cache.W) <: Dual ? "Dual" : eltype(integ.cache.W), "}")
    catch e
        "could not read: " * sprint(showerror, e)[1:min(end, 200)]
    end
    say("dual linear solve: $linpath")

    # DualSparseKLU against LinearSolve's own dual path (DualLinearCache), one pass of
    # CHUNK partials over a short window, LIB_T s, where the library's cost allows it.
    lib = let tend = LIB_T, sv = collect(10.0:10.0:LIB_T)
        θd = seed(θstar(B), 1:CHUNK)
        dsol(alg) = (t = @elapsed (s = solve(problem(m, slot, θd; tend), alg; abstol = ABSTOL,
                                          reltol = RELTOL, saveat = sv, save_start = false,
                                          maxiters = MAXITERS));
                     SciMLBase.successful_retcode(s) || error("dual solve: $(s.retcode)"); (s, t))
        dsol(solver(eltype(θd))); dsol(solver_library())           # compile both
        (s1, t1) = dsol(solver(eltype(θd)))
        (s2, t2) = dsol(solver_library())
        U1 = reduce(hcat, s1.u); U2 = reduce(hcat, s2.u)
        num = maximum(maximum(abs, partials(U1[k]) .- partials(U2[k])) for k in eachindex(U1))
        den = maximum(maximum(abs, partials(U2[k])) for k in eachindex(U2))
        vrel = maximum(abs(value(U1[k]) - value(U2[k])) / max(abs(value(U2[k])), 1.0) for k in eachindex(U1))
        (t_end = tend, ours_s = t1, library_s = t2, steps_ours = s1.stats.naccept + s1.stats.nreject,
         steps_library = s2.stats.naccept + s2.stats.nreject, partials_rel = num / den, value_rel = vrel)
    end
    say(@sprintf("dual solve to %.0f s: DualSparseKLU %.3f s (%d steps), DualLinearCache %.3f s (%d steps); partials max |Δ| / max |∂| %.3e, values %.3e",
                 lib.t_end, lib.ours_s, lib.steps_ours, lib.library_s, lib.steps_library,
                 lib.partials_rel, lib.value_rel))

    # AV3: the dual gradient against central differences, at tight tolerances.
    tight = (abstol = ABSTOL_TIGHT, reltol = RELTOL_TIGHT)
    data = make_data(m, slot, B; tight...)
    θ = θstar(B)
    say("tight gradient (B = $B)")
    gt = gradient_passes(m, slot, data, θ; tight...)
    fdidx = [1, 10, 20]
    fd = map(fdidx) do b
        e = zeros(B); e[b] = FD_H
        (first(objective(m, slot, data, θ .+ e; tight...)) -
         first(objective(m, slot, data, θ .- e; tight...))) / (2FD_H)
    end
    fdrel = [abs(gt.grad[b] - fd[i]) / abs(fd[i]) for (i, b) in enumerate(fdidx)]
    for (i, b) in enumerate(fdidx)
        say(@sprintf("AV3 θ_%d (reaction %d): dual %.12e, FD %.12e, rel %.3e", b,
                     ranked[b].reaction, gt.grad[b], fd[i], fdrel[i]))
    end
    fd_pass = all(<=(FD_RTOL), fdrel)

    # The pass-wise gradient against ForwardDiff.gradient, at the timing tolerances.
    data_t = make_data(m, slot, B)
    gp = gradient_passes(m, slot, data_t, θ)
    obj = z -> first(objective(m, slot, data_t, z))
    gfd = ForwardDiff.gradient(obj, θ, ForwardDiff.GradientConfig(obj, θ, ForwardDiff.Chunk{CHUNK}()))
    fdiff = maximum(abs, gp.grad .- gfd)
    say(@sprintf("pass-wise vs ForwardDiff.gradient (chunk %d): max |Δ| %.3e", CHUNK, fdiff))

    res = (env = env, n = m.n, nspecies = nsp, n_first, n_second, nnz_table = m.nnz_table,
           nnz_proto = nnz(m.jac_proto), rhs_maxrel = maxrel, rhs_maxrel_abs = maxrel_abs,
           rhs_nonzero = nonzero_ref, rhs_zero = zero_ref, rhs_zero_mismatch = zero_mismatch,
           jx_rel = jrel, jx_zero = jzero, rhs_pass, ranked, nflux = count(>(0), flux),
           linpath, lib, fdidx, grad_tight = gt.grad, fd, fdrel, fd_pass, passwise_vs_fd = fdiff,
           tight_stats = [p.stats for p in gt.passes])
    serialize(joinpath(out, "check.jls"), res)
    open(joinpath(out, "check.md"), "w") do io
        header(io, "check", env)
        println(io, "## The system\n")
        println(io, "| quantity | value |\n|---|---|")
        println(io, "| CME species in the table | $nsp |")
        println(io, "| states (species a first-order reaction reads or changes) | $(m.n) |")
        println(io, "| first-order reactions (LM type 1) | $n_first |")
        println(io, "| second-order reactions (type 2), left out | $n_second |")
        println(io, "| other types | $n_other |")
        println(io, "| Jacobian non-zeros from the table | $(m.nnz_table) |")
        println(io, "| … with the structural diagonal | $(nnz(m.jac_proto)) |")
        println(io, "| reactions with flux over the window | $(count(>(0), flux)) |")
        println(io)
        println(io, "## Checks\n")
        println(io, @sprintf("- RHS at x0 vs the table's propensities summed in BigFloat: %d non-zero entries, max relative error %.3e (bound %.0e); %d entries exactly zero in both (%d mismatches). **%s**",
                             nonzero_ref, maxrel, RHS_RTOL, zero_ref, zero_mismatch, rhs_pass ? "PASS" : "FAIL"))
        println(io, @sprintf("- J·x0 vs the RHS: max relative difference %.3e; zeros kept: %s", jrel, jzero))
        println(io, "- Dual linear solve (the solver's linear cache under a dual state): $linpath")
        println(io, @sprintf("- DualSparseKLU vs LinearSolve's DualLinearCache, one pass of %d partials solved to %.0f s: %.3f s (%d steps) against %.3f s (%d steps); partials differ by %.3e of their largest, values by %.3e relative",
                             CHUNK, lib.t_end, lib.ours_s, lib.steps_ours, lib.library_s, lib.steps_library, lib.partials_rel, lib.value_rel))
        println(io, @sprintf("- AV3, B = %d, abstol %.0e, reltol %.0e, h = %.0e:", B, ABSTOL_TIGHT, RELTOL_TIGHT, FD_H))
        println(io, "\n| θ | reaction | dual | central FD | relative |\n|---|---|---|---|---|")
        for (i, b) in enumerate(fdidx)
            println(io, @sprintf("| θ_%d | %d | %.12e | %.12e | %.3e |", b, ranked[b].reaction,
                                 gt.grad[b], fd[i], fdrel[i]))
        end
        println(io, "\n  **$(fd_pass ? "PASS" : "FAIL")** (bound $FD_RTOL)")
        println(io, @sprintf("- Pass-wise gradient vs `ForwardDiff.gradient` with Chunk{%d}, timing tolerances: max |Δ| = %.3e", CHUNK, fdiff))
        println(io, "\n## The block (ranked by flux integrated over the window at θ = 0)\n")
        println(io, "| rank | reaction | reactant | k⁰ (1/s) | ∫ flux dt |\n|---|---|---|---|---|")
        for r in ranked
            println(io, @sprintf("| %d | %d | %s | %.6g | %.6g |", r.rank, r.reaction, r.reactant, r.k0, r.flux))
        end
        println(io, "\n```\n", env.versioninfo, "```")
    end
    say("check: RHS $(rhs_pass ? "PASS" : "FAIL"), AV3 $(fd_pass ? "PASS" : "FAIL")")
    occursin("LinearSolveFunction", linpath) || error("the dual solve did not reach DualSparseKLU: $linpath")
    lib.partials_rel == 0 && lib.value_rel == 0 ||
        error("DualSparseKLU and DualLinearCache differ: $(lib.partials_rel), $(lib.value_rel)")
    fdiff == 0 || error("pass-wise gradient differs from ForwardDiff.gradient by $fdiff")
    return rhs_pass && fd_pass
end

function stage_time(out, threaded::Bool)
    stage = threaded ? "threaded" : "time"
    env = environment()
    threaded && Threads.nthreads() < length(passes(maximum(BLOCKS))) &&
        error("threaded needs ≥ $(length(passes(maximum(BLOCKS)))) threads, has $(Threads.nthreads())")
    say("$stage on $(env.host) ($(env.cpu)), job $(env.job), threads $(env.threads), BLAS $(env.blas)")
    m = load_model()
    flux, order = rank_reactions(m)
    results = []
    for B in BLOCKS
        slot = block_slot(m, order, B)
        data = make_data(m, slot, B)
        θ = θstar(B)
        # Warm-up: compile both the primal and the dual paths.
        tc = @elapsed begin
            objective(m, slot, data, θ)
            gradient_passes(m, slot, data, θ; threaded)
        end
        say(@sprintf("B = %d: warm-up %.1f s", B, tc))
        prim = map(1:NREPEAT) do r
            local v, st
            t = @elapsed ((v, st) = objective(m, slot, data, θ))
            say(@sprintf("B = %d primal %d: %.3f s, ℓ %.10e, %s", B, r, t, v, statsnt(st)))
            (seconds = t, value = v, stats = statsnt(st))
        end
        grads = map(1:NREPEAT) do r
            g = gradient_passes(m, slot, data, θ; threaded)
            say(@sprintf("B = %d gradient %d: %.3f s (passes %s), |g| %.6e", B, r, g.wall,
                         join((@sprintf("%.3f", p.seconds) for p in g.passes), ", "), norm(g.grad)))
            g
        end
        push!(results, (B = B, npasses = length(passes(B)), warmup = tc, primal = prim,
                        grads = [(wall = g.wall, pass_seconds = [p.seconds for p in g.passes],
                                  pass_stats = [p.stats for p in g.passes], grad = g.grad)
                                 for g in grads]))
    end
    serialize(joinpath(out, "$stage.jls"), (env = env, results = results))
    # Against the serial gradient, if that stage ran first in this outdir.
    vs_serial = nothing
    if threaded && isfile(joinpath(out, "time.jls"))
        ser = deserialize(joinpath(out, "time.jls"))
        vs_serial = [maximum(abs, r.grads[1].grad .- s.grads[1].grad)
                     for (r, s) in zip(results, ser.results)]
    end
    open(joinpath(out, "$stage.md"), "w") do io
        header(io, stage, env)
        println(io, "Medians of $NREPEAT after a warm-up. A pass is one dual solve of $CHUNK partials.\n")
        println(io, "| B | passes | primal (s) | gradient (s) | gradient / primal | per pass (s) | primal steps (acc+rej) | pass steps | pass njacs | pass nw | s per pass step |")
        println(io, "|---|---|---|---|---|---|---|---|---|---|---|")
        for r in results
            tp = median([p.seconds for p in r.primal])
            tg = median([g.wall for g in r.grads])
            tpp = median(vcat([g.pass_seconds for g in r.grads]...))
            ps = r.grads[1].pass_stats[1]; pr = r.primal[1].stats
            println(io, @sprintf("| %d | %d | %.3f | %.3f | %.2f | %.3f | %d+%d | %d+%d | %d | %d | %.4f |",
                                 r.B, r.npasses, tp, tg, tg / tp, tpp, pr.naccept, pr.nreject,
                                 ps.naccept, ps.nreject, ps.njacs, ps.nw, tpp / (ps.naccept + ps.nreject)))
        end
        println(io, "\nEvery repeat (seconds):\n")
        for r in results
            println(io, "- B = $(r.B): primal ", join((@sprintf("%.3f", p.seconds) for p in r.primal), ", "),
                    "; gradient ", join((@sprintf("%.3f", g.wall) for g in r.grads), ", "))
        end
        if vs_serial !== nothing
            println(io, "\nThreaded gradient vs serial (first repeat), max |Δ| per B: ",
                    join((@sprintf("%.3e", d) for d in vs_serial), ", "))
        end
        println(io, "\n```\n", env.versioninfo, "```")
    end
    say("$stage done")
    vs_serial === nothing || all(==(0), vs_serial) ||
        error("threaded gradient differs from serial: $vs_serial")
    return true
end

function stage_library(out)
    env = environment()
    say("library on $(env.host) ($(env.cpu)), job $(env.job)")
    m = load_model()
    flux, order = rank_reactions(m)
    B = 20
    slot = block_slot(m, order, B)
    θd = seed(θstar(B), 1:CHUNK)
    dsol(alg; tend = T_END, saveat = SAVES) =
        (t = @elapsed (s = solve(problem(m, slot, θd; tend), alg; abstol = ABSTOL, reltol = RELTOL,
                                 saveat, save_start = false, maxiters = MAXITERS));
         SciMLBase.successful_retcode(s) || error("dual solve: $(s.retcode)"); (s, t))
    for alg in (solver(eltype(θd)), solver_library())
        dsol(alg; tend = LIB_T, saveat = collect(10.0:10.0:LIB_T))
    end
    ours = [dsol(solver(eltype(θd))) for _ in 1:NREPEAT]
    say(@sprintf("DualSparseKLU, whole window: %s s", join((@sprintf("%.3f", o[2]) for o in ours), ", ")))
    (s2, t2) = dsol(solver_library())
    s1 = ours[1][1]
    say(@sprintf("DualLinearCache, whole window: %.1f s, %d steps", t2, s2.stats.naccept + s2.stats.nreject))
    U1 = reduce(hcat, s1.u); U2 = reduce(hcat, s2.u)
    num = maximum(maximum(abs, partials(U1[k]) .- partials(U2[k])) for k in eachindex(U1))
    den = maximum(maximum(abs, partials(U2[k])) for k in eachindex(U2))
    vrel = maximum(abs(value(U1[k]) - value(U2[k])) / max(abs(value(U2[k])), 1.0) for k in eachindex(U1))
    st1 = statsnt(s1.stats); st2 = statsnt(s2.stats)
    t1 = median([o[2] for o in ours])
    res = (env = env, ours_s = [o[2] for o in ours], library_s = t2, stats_ours = st1,
           stats_library = st2, partials_rel = num / den, value_rel = vrel)
    serialize(joinpath(out, "library.jls"), res)
    open(joinpath(out, "library.md"), "w") do io
        header(io, "library", env)
        println(io, "One dual solve of $CHUNK partials (B = 20, θ*, partials 1–$CHUNK) over the whole window, saved on the $(length(SAVES))-point grid; no misfit.\n")
        println(io, "| linear solve | seconds | steps (acc+rej) | njacs | nw | s per step |\n|---|---|---|---|---|---|")
        println(io, @sprintf("| DualSparseKLU (median of %d) | %.3f | %d+%d | %d | %d | %.4f |", NREPEAT, t1,
                             st1.naccept, st1.nreject, st1.njacs, st1.nw, t1 / (st1.naccept + st1.nreject)))
        println(io, @sprintf("| LinearSolve DualLinearCache (1 run) | %.1f | %d+%d | %d | %d | %.3f |", t2,
                             st2.naccept, st2.nreject, st2.njacs, st2.nw, t2 / (st2.naccept + st2.nreject)))
        println(io, @sprintf("\nRatio %.0f×. Partials differ by %.3e of their largest, values by %.3e relative.", t2 / t1, num / den, vrel))
        println(io, "\nDualSparseKLU runs: ", join((@sprintf("%.3f", o[2]) for o in ours), ", "))
        println(io, "\n```\n", env.versioninfo, "```")
    end
    say("library done")
    num == 0 && vrel == 0 || error("DualSparseKLU and DualLinearCache differ: $num, $vrel")
    return true
end

function main(args)
    length(args) == 2 || error("usage: fullcell_timing_meanfield.jl <check|time|threaded|library> <out>")
    stage, out = args
    mkpath(out)
    ok = stage == "check" ? stage_check(out) :
         stage == "time" ? stage_time(out, false) :
         stage == "threaded" ? stage_time(out, true) :
         stage == "library" ? stage_library(out) : error("unknown stage $stage")
    ok || exit(1)
end

main(ARGS)
