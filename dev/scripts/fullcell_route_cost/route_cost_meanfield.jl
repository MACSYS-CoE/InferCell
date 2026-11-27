# KIM-195 item 3: what one gradient of the full syn3A cell's mean field costs, on a
# realistic system: KIM-184's first-order surrogate extended with the second-order
# reactions (product of means), the 175 metabolic fluxes and the minute rebuilds
# (model.jl). Solved under ForwardDiff dual numbers by Rodas5P with a sparse KLU solve,
# with two dual linear solves: DualSparseKLU (KIM-184's, in this file) and LinearSolve
# 3.87's own dual path with KLUFactorization explicit.
#
# Usage: julia --project=dev/scripts/fullcell_route_cost \
#            dev/scripts/fullcell_route_cost/route_cost_meanfield.jl <stage> <out>
# Julia 1.11.9: LinearSolve 3.87.0 does not load on 1.10.5 (its extensions are
# triggered by SparseArrays, which it also lists as a strong dependency; 1.10's loader
# fails with KeyError "SparseArrays"). Submit through route_cost_meanfield.sh, from a
# detached worktree at a pushed commit.
#
# The gradient. θ_b, b = 1..B, are log multipliers on a block (frame.md §1.2): GE, the B
# CME reactions with the largest flux integrated over the window at θ = 0; MET, the B
# kcatF constants whose reactions carry the largest integrated |flux|; BOTH, the top
# B/2 of each. The blocks of 10, 20 and 40 are nested. The scalar is KIM-184's: a
# Gaussian misfit over every state at the 52 minute saves (60 s to 3,120 s; solved to
# COREA_FIT_S = 3,150 s), against the solution at θ = 0, variance max(d, 1), at
# θ*_b = 0.05·(−1)^b. Forward mode, CHUNK partials a pass, one dual solve a pass.
#
# Stages (each writes <out>/<stage>.md and <stage>.jls):
#   check     sizes; the rate rules, p, y0 and fluxes at x0 against SpeedyCell's own
#             values (<out>/reference/, from reference_speedycell.jl; 1e-12 relative);
#             the rules at x0 against reactions.csv's `rate` column (reported); the RHS
#             at x0 against an independent BigFloat sum of the table's propensities, the
#             counters paid, plus the ODE part from SpeedyCell's fluxes (1e-10
#             relative); J·v against central differences of the RHS (1e-6); the ranked
#             blocks; the dual gradient against central differences on 3 parameters per
#             placement at tight tolerances (1e-6 of max|g| over
#             the block; per-component errors at h = 1e-4 and 1e-3 reported); DualSparseKLU against
#             LinearSolve's dual path on one pass over the first minute (stated
#             tolerance 1e-8 relative)
#   time      serial: the primal solve, and the gradient for each placement and B, each
#             the median of NREPEAT after a warm-up, with per-pass times and step counts
#   threaded  the same, the passes on threads (needs ≥ 4 threads)
#   library   B = 20, one dual pass of CHUNK partials over the whole window through
#             LinearSolve 3.87's own dual path, against DualSparseKLU's, for GE and MET
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

include(joinpath(@__DIR__, "model.jl"))

const SMOKE = get(ENV, "SMOKE", "0") == "1"
const COREA_FIT_S = 3150.0
const T_END = SMOKE ? 300.0 : COREA_FIT_S
const ABSTOL = 1e-10
const RELTOL = 1e-8
const MAXITERS = 10^5
const ABSTOL_TIGHT = 1e-14
const RELTOL_TIGHT = 1e-12
const FD_H = 1e-4
const FD_RTOL = 1e-6                    # scale-relative: |dual − FD| / max|g| (see stage_check)
const FD_H_WIDE = 1e-3                 # a second step, reported: out of the objective's noise floor
const RHS_RTOL = 1e-10
const REF_RTOL = 1e-12
const JV_RTOL = 1e-6
const LIB_RTOL = 1e-8
const CHUNK = 10
const BLOCKS = (10, 20, 40)
const PLACEMENTS = (:GE, :MET, :BOTH)
const NREPEAT = SMOKE ? 2 : 3
const NREPEAT_THREADED = 2
const T_FD = SMOKE ? 120.0 : 300.0     # the window of the finite-difference check: 5 rebuilds
const T_WARM = 60.0                    # warm-ups: one segment, to compile
const LOOSE = (abstol = 1e-6, reltol = 1e-6)   # one row at looser tolerances (a lever for KIM-194)

"The segment ends: every minute to T_END, then T_END itself."
function segment_ends(tend)
    e = collect(60.0:60.0:tend)
    (isempty(e) || last(e) < tend) && push!(e, tend)
    return e
end
"Which segment ends are saves: the whole minutes."
issave(t) = t % 60 == 0

# ---------------------------------------------------------------------------------------
# The linear solve under duals: KIM-184's DualSparseKLU, unchanged but for its comments.
# KLU on the value part, refactored only when the values change, the primal solved,
# then all partials at once as an n×N right-hand side R = ∂b − ∂W·x₀, from the non-zeros.

const KLU = Base.get_extension(LinearSolve, :LinearSolveSparseArraysExt).KLU
const LSFD = Base.get_extension(LinearSolve, :LinearSolveForwardDiffExt)

mutable struct DualSparseKLU
    F::Any
    vals::Vector{Float64}
    x0::Vector{Float64}
    R::Matrix{Float64}
    nfact::Int
    nsolve::Int
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
SciMLBase.init(prob::LSFD.DualAbstractLinearProblem, alg::DualKLUAlg, args...; kwargs...) =
    LinearSolve.__init(prob, alg, args...; kwargs...)
LinearSolve.needs_concrete_A(::DualKLUAlg) = true

solver(::Type{<:Dual}) = Rodas5P(linsolve = LinearSolveFunction(DualSparseKLU()))
solver(::Type) = Rodas5P(linsolve = KLUFactorization())
solver_library() = Rodas5P(linsolve = KLUFactorization())

struct K195Tag end
const TAG = typeof(ForwardDiff.Tag(K195Tag(), Float64))

say(s) = (println(Dates.format(now(), "HH:MM:SS"), "  ", s); flush(stdout))

# ---------------------------------------------------------------------------------------
# The solve, one minute at a time

const Stats = Dict{Symbol, Int}
addstats!(acc, s) = for k in (:nf, :njacs, :nw, :nsolve, :naccept, :nreject)
    acc[k] = get(acc, k, 0) + getfield(s, k)
end

"""Solve the mean field from x0 at θ to `tend`: (saves, the final state, stats, ok).
`alg(T)` gives each segment's solver (a fresh linear solve per segment)."""
function simulate(m::Model, blk::Block, θ; tend = T_END, abstol = ABSTOL, reltol = RELTOL,
                  alg = solver, dense = false, smooth = true, onsegment = nothing)
    T = eltype(θ) <: Dual ? eltype(θ) : Float64
    x = convert(Vector{T}, copy(m.x0))
    kode = ode_constants(m, blk, θ)
    saves = Vector{Vector{T}}()
    st = Stats()
    t0 = 0.0
    for (s, t1) in enumerate(segment_ends(tend))
        mn = rebuild!(x, m, blk, θ, kode; first = s == 1, smooth)
        f = SegRHS(m, mn)
        fun = ODEFunction{true, SciMLBase.FullSpecialize}(f; jac = SegJac(f),
                                                         jac_prototype = copy(m.proto), tgrad = seg_tgrad!)
        prob = ODEProblem(fun, x, (t0, t1))
        sol = solve(prob, alg(T); abstol, reltol, maxiters = MAXITERS, dense,
                    save_everystep = dense, save_start = dense)
        SciMLBase.successful_retcode(sol) || error("segment $(s) ($(t0)–$(t1) s): $(sol.retcode)")
        addstats!(st, sol.stats)
        onsegment === nothing || onsegment(sol, mn)
        x = copy(sol.u[end])
        issave(t1) && push!(saves, copy(x))
        t0 = t1
    end
    return saves, x, st
end

"""The data: the primal saves at θ = 0, and their variances."""
function make_data(m, blk; kw...)
    saves, _, _ = simulate(m, blk, zeros(blk.B); kw...)
    d = reduce(hcat, saves)
    return (d = d, w = 1 ./ max.(d, 1.0))
end

function objective(m, blk, data, θ; kw...)
    saves, _, st = simulate(m, blk, θ; kw...)
    T = eltype(θ)
    ℓ = zero(T)
    for (i, u) in enumerate(saves)
        @inbounds for s in eachindex(u)
            r = u[s] - data.d[s, i]
            ℓ -= r * r * data.w[s, i] / 2
        end
    end
    return ℓ, st
end

passes(B) = [i:min(i + CHUNK - 1, B) for i in 1:CHUNK:B]
seed(θ, idx) = [Dual{TAG}(θ[i], Partials(ntuple(c -> (c <= length(idx) && idx[c] == i) ? 1.0 : 0.0,
                                                  Val(CHUNK)))) for i in eachindex(θ)]

function dual_pass(m, blk, data, θ, idx; kw...)
    t = @elapsed begin
        v, st = objective(m, blk, data, seed(θ, idx); kw...)
    end
    return (value = value(v), partials = collect(partials(v))[1:length(idx)], seconds = t, stats = st)
end

function gradient_passes(m, blk, data, θ; threaded = false, kw...)
    ps = passes(length(θ))
    out = Vector{Any}(undef, length(ps))
    wall = @elapsed begin
        if threaded
            Threads.@threads :static for c in eachindex(ps)
                out[c] = dual_pass(m, blk, data, θ, ps[c]; kw...)
            end
        else
            for c in eachindex(ps)
                out[c] = dual_pass(m, blk, data, θ, ps[c]; kw...)
            end
        end
    end
    g = zeros(length(θ))
    for (c, idx) in enumerate(ps)
        g[idx] .= out[c].partials
    end
    return (grad = g, wall = wall, passes = out)
end

θstar(B) = [0.05 * (-1)^b for b in 1:B]

# ---------------------------------------------------------------------------------------
# The blocks

"""Each CME reaction's flux and each ODE reaction's |flux| (counts), integrated over the
window at θ = 0 by 3-point Gauss–Legendre on the solver's own steps in every segment."""
function integrated_fluxes(m::Model)
    blk = Block(m)
    cme = zeros(m.nr); ode = zeros(size(m.ode.S, 2))
    nodes = (-sqrt(3 / 5), 0.0, sqrt(3 / 5)); weights = (5 / 9, 8 / 9, 5 / 9)
    V = zeros(size(m.ode.S, 2))
    function acc(sol, mn)
        for i in 1:(length(sol.t) - 1)
            a, b = sol.t[i], sol.t[i + 1]
            h = (b - a) / 2
            for (z, w) in zip(nodes, weights)
                x = sol(a + h * (1 + z))
                for j in 1:m.nr
                    cme[j] += w * h * (m.rb[j] > 0 ? mn.k[j] * x[m.ra[j]] * x[m.rb[j]] : mn.k[j] * x[m.ra[j]])
                end
                m.ode.flux!(V, [x[i] / mn.c for i in m.ode.cme], mn.p, mn.kode)
                ode .+= (w * h * mn.c) .* abs.(V)
            end
        end
    end
    simulate(m, blk, Float64[]; dense = true, onsegment = acc)
    return cme, ode
end

"""The ranked candidates: CME reactions by flux; kcatF constants by their reaction's
|flux|. Ties by index."""
function rank_blocks(m::Model)
    cme, ode = integrated_fluxes(m)
    ge = sortperm(collect(zip(-cme, 1:m.nr)))
    kcat = [i for (i, nm) in enumerate(m.ode.k_names) if endswith(nm, ".kcatF")]
    rxn_of = Dict(id => j for (j, id) in enumerate(m.ode.reaction_ids))
    kflux = [ode[rxn_of[split(m.ode.k_names[i], '.')[1]]] for i in kcat]
    met = kcat[sortperm(collect(zip(-kflux, kcat)))]
    return (ge = ge, met = met, cme_flux = cme, ode_flux = ode, kflux = Dict(zip(kcat, kflux)))
end

function make_block(m::Model, r, place::Symbol, B::Int)
    blk = Block(m)
    cs = zeros(Int, m.nr); ks = zeros(Int, length(m.ode.k_default))
    if place === :GE
        for (b, j) in enumerate(r.ge[1:B]); cs[j] = b; end
    elseif place === :MET
        for (b, i) in enumerate(r.met[1:B]); ks[i] = b; end
    else
        h = B ÷ 2
        for (b, j) in enumerate(r.ge[1:h]); cs[j] = b; end
        for (b, i) in enumerate(r.met[1:(B - h)]); ks[i] = h + b; end
    end
    return Block(B, cs, ks)
end

# ---------------------------------------------------------------------------------------
# Environment and output

function environment()
    io = IOBuffer(); versioninfo(io)
    cpu = try
        strip(split(first(filter(l -> startswith(l, "Model name"), readlines(`lscpu`))), ':')[2])
    catch
        "unknown"
    end
    return (host = gethostname(), cpu = cpu, julia = string(VERSION), threads = Threads.nthreads(),
            blas = BLAS.get_num_threads(), ngc = Threads.ngcthreads(), job = get(ENV, "SLURM_JOB_ID", "none"),
            smoke = SMOKE, t_end = T_END, linearsolve = string(pkgversion(LinearSolve)),
            versioninfo = String(take!(io)))
end

function header(io, stage, env)
    println(io, "# KIM-195 `$stage`", SMOKE ? " (SMOKE: a plumbing pass, not a result)" : "")
    println(io)
    println(io, "- job $(env.job) on $(env.host), $(env.cpu); Julia $(env.julia); LinearSolve $(env.linearsolve)")
    println(io, "- threads $(env.threads), BLAS threads $(env.blas), GC threads $(env.ngc)")
    println(io, "- window 0–$(env.t_end) s in minute segments, saves at the whole minutes; Rodas5P, " *
                "KLUFactorization (primal) and DualSparseKLU (dual, unless a row says the library path); " *
                "abstol $ABSTOL, reltol $RELTOL; chunk $CHUNK")
    println(io)
end

statline(st) = @sprintf("steps %d (rejected %d), f %d, J %d, W %d, solves %d",
                        get(st, :naccept, 0) + get(st, :nreject, 0), get(st, :nreject, 0), get(st, :nf, 0),
                        get(st, :njacs, 0), get(st, :nw, 0), get(st, :nsolve, 0))

readref(out, f) = parse.(Float64, readlines(joinpath(out, "reference", f)))
relerr(a, b) = abs(a - b) / max(abs(b), floatmin(Float64))

# ---------------------------------------------------------------------------------------
# Stages

function stage_check(out)
    env = environment()
    say("check on $(env.host) ($(env.cpu)), job $(env.job)")
    m = get_model()
    ntype1 = count(==(1), m.lm_type); ntype2 = count(==(2), m.lm_type)
    say("states $(m.n); CME reactions $(m.nr) (type 1 $ntype1, type 2 $ntype2); ODE species $(length(m.ode.ids)), " *
        "fluxes $(size(m.ode.S, 2)), enzymes $(length(m.ode.enzyme_nodes)), named constants $(length(m.ode.k_default))")
    say("DAG nodes needed $(length(m.dag.ids)) of $(length(m.dag.op)); Jacobian non-zeros $(nnz(m.proto))")
    res = Dict{Symbol, Any}()

    # 1. The generic code at x0 against SpeedyCell's own values.
    x = copy(m.x0)
    blk0 = Block(m)
    kode0 = ode_constants(m, blk0, Float64[])
    mn = rebuild!(x, m, blk0, Float64[], kode0; first = true, smooth = false)
    cellv_given = m.x0[m.i_cellv]; cellv_ours = x[m.i_cellv]
    rates_ref = readref(out, "rates.csv"); p_ref = readref(out, "p.csv"); y0_ref = readref(out, "y0.csv")
    flux_ref = readref(out, "flux.csv"); fluxp_ref = readref(out, "flux_params.csv")
    r_rates = maximum(relerr.(mn.k, rates_ref)); n_rates_exact = count(mn.k .== rates_ref)
    r_p = maximum(relerr.(mn.p, p_ref))
    y0 = [x[i] / mn.c for i in m.ode.cme]
    nzy = findall(!=(0), y0_ref)
    r_y0 = maximum(relerr.(y0[nzy], y0_ref[nzy])); y0_zeros = all(y0[i] == 0 for i in eachindex(y0) if y0_ref[i] == 0)
    V = zeros(size(m.ode.S, 2)); m.ode.flux!(V, y0_ref, p_ref, m.ode.k_default)
    nzf = findall(!=(0), flux_ref)
    r_flux = maximum(relerr.(V[nzf], flux_ref[nzf])); flux_zeros = all(V[i] == 0 for i in eachindex(V) if flux_ref[i] == 0)
    say(@sprintf("vs SpeedyCell at x0: rates max rel %.2e (%d of %d bitwise), p %.2e, y0 %.2e (zeros kept %s), fluxes %.2e (zeros kept %s); CellV given %g, ours %g",
                 r_rates, n_rates_exact, m.nr, r_p, r_y0, y0_zeros, r_flux, flux_zeros, cellv_given, cellv_ours))
    ref_pass = max(r_rates, r_p, r_y0, r_flux) <= REF_RTOL && y0_zeros && flux_zeros && cellv_ours == cellv_given
    res[:reference] = (rates = r_rates, rates_bitwise = n_rates_exact, p = r_p, y0 = r_y0, flux = r_flux,
                       fluxp_vs_flux_bitwise = fluxp_ref == flux_ref, pass = ref_pass)
    # The rules at x0 against reactions.csv's `rate` column (minute 1, a different script).
    d = relerr.(mn.k, m.k_rate_col)
    differ = findall(>(1e-12), d)
    say(@sprintf("rules at x0 vs reactions.csv `rate`: %d of %d differ by > 1e-12 relative (max %.3e)",
                 length(differ), m.nr, isempty(d) ? 0.0 : maximum(d)))
    res[:rate_column] = (ndiffer = length(differ), maxrel = maximum(d),
                         differ_types = isempty(differ) ? Int[] : unique(m.lm_type[differ]))

    # 2. The RHS at x0 against an independent BigFloat sum: the table's propensities with
    #    SpeedyCell's rates, by species name; then the counters paid (writeback's lists);
    #    then the ODE part from SpeedyCell's own fluxes.
    setprecision(BigFloat, 256)
    rows = csvrows(joinpath(MODEL, "cme", "reactions.csv"))
    sp = csvrows(joinpath(MODEL, "cme", "species.csv"))
    nm = Dict(parse(Int, r[1]) => String(r[2]) for r in sp)
    xb = Dict(String(r[2]) => parse(BigFloat, r[3]) for r in sp)
    ref = Dict{String, BigFloat}(); absum = Dict{String, BigFloat}()
    addt!(s, v) = (ref[s] = get(ref, s, big(0)) + v; absum[s] = get(absum, s, big(0)) + abs(v))
    for (j, r) in enumerate(rows)
        re = parse.(Int, split(r[4], ';'))
        a = BigFloat(rates_ref[j]) * prod(xb[nm[i]] for i in re)
        for pr in split(r[5], ';')
            s, dd = split(pr, ':')
            addt!(nm[parse(Int, s)], parse(BigFloat, dd) * a)
        end
    end
    pays = vcat([(c, occursin("translat", c) ? ["M_gtp_c" => -1, "M_gdp_c" => 1, "M_pi_c" => 1] :
                                               ["M_atp_c" => -1, "M_adp_c" => 1, "M_pi_c" => 1]) for c in ATP_HYDRO_COUNTERS],
                [(c, [mm => -1, "M_ppi_c" => 1]) for (c, mm) in vcat(NTP_COUNTERS, DNTP_COUNTERS)],
                [(c, [mm => 1]) for (c, mm) in NMP_COUNTERS],
                [(c, [ch => -1, un => 1]) for (c, ch, un) in AATRNA_COUNTERS])
    for (c, lst) in pays
        dc = get(ref, c, big(0))
        for (t, f) in lst
            addt!(t, f * dc)
        end
        ref[c] = big(0)
    end
    cb = BigFloat(mn.c)
    st = csvrows(joinpath(MODEL, "ode", "stoichiometry.csv"))
    for r in st
        addt!(r[1], cb * parse(BigFloat, r[3]) * BigFloat(flux_ref[parse(Int, r[2])]))
    end
    f = SegRHS(m, mn)
    du = zeros(m.n); f(du, x, nothing, 0.0)
    nonzero_ref = 0; zero_ref = 0; zero_mismatch = 0; maxrel = 0.0; maxrel_abs = 0.0
    for (i, s) in enumerate(m.names)
        rv = get(ref, s, big(0))
        if iszero(rv)
            zero_ref += 1
            abs(du[i]) <= 1e-12 * Float64(get(absum, s, big(1))) || (zero_mismatch += 1)
        else
            nonzero_ref += 1
            maxrel = max(maxrel, Float64(abs(du[i] - rv) / abs(rv)))
            maxrel_abs = max(maxrel_abs, Float64(abs(du[i] - rv) / absum[s]))
        end
    end
    say(@sprintf("RHS vs BigFloat: %d non-zero, max rel %.3e (vs Σ|terms| %.3e); %d zero in the reference, %d not ~0 in the RHS",
                 nonzero_ref, maxrel, maxrel_abs, zero_ref, zero_mismatch))
    rhs_pass = maxrel <= RHS_RTOL && zero_mismatch == 0
    res[:rhs] = (nonzero = nonzero_ref, maxrel = maxrel, maxrel_abs = maxrel_abs, zero = zero_ref,
                 zero_mismatch = zero_mismatch, pass = rhs_pass)

    # 3. The Jacobian against central differences of the RHS along a direction.
    xs = copy(x) .+ 1.0                       # off the zero boundary, so every term is live
    vdir = [((i * 7919) % 13 - 6) / 6 * max(xs[i], 1.0) * 1e-3 for i in 1:m.n]
    J = copy(m.proto); SegJac(f)(J, xs, nothing, 0.0)
    jv = J * vdir
    h = 1e-3
    dp = zeros(m.n); dm = zeros(m.n)
    f(dp, xs .+ h .* vdir, nothing, 0.0); f(dm, xs .- h .* vdir, nothing, 0.0)
    fd = (dp .- dm) ./ (2h)
    scale = maximum(abs, fd)
    jv_rel = maximum(abs.(jv .- fd)) / scale
    say(@sprintf("J·v vs central differences: max |Δ| / max |fd| = %.3e", jv_rel))
    res[:jv] = (rel = jv_rel, pass = jv_rel <= JV_RTOL)

    # 4. The blocks.
    tr = @elapsed (r = rank_blocks(m))
    say(@sprintf("ranked in %.1f s: %d of %d CME reactions and %d of %d kcatF carry flux",
                 tr, count(>(0), r.cme_flux), m.nr, count(>(0), values(r.kflux)), length(r.met)))
    res[:ranking] = (ge = [(rank = b, reaction = j, type = m.lm_type[j], reactants = (m.names[m.ra[j]],
                             m.rb[j] > 0 ? m.names[m.rb[j]] : ""), flux = r.cme_flux[j])
                            for (b, j) in enumerate(r.ge[1:maximum(BLOCKS)])],
                     met = [(rank = b, constant = m.ode.k_names[i], value = m.ode.k_default[i], flux = r.kflux[i])
                            for (b, i) in enumerate(r.met[1:maximum(BLOCKS)])])
    serialize(joinpath(out, "ranking.jls"), (ge = r.ge, met = r.met))

    # 5. The dual gradient against central differences, tight tolerances, B = 20, over
    #    the first T_FD s (the derivative through T_FD/60 rebuilds).
    tight = (abstol = ABSTOL_TIGHT, reltol = RELTOL_TIGHT, tend = T_FD)
    B = 20; θ = θstar(B); fdidx = [1, 10, 20]
    fdres = Dict{Symbol, Any}()
    for place in PLACEMENTS
        blk = make_block(m, r, place, B)
        data = make_data(m, blk; tight...)
        g = gradient_passes(m, blk, data, θ; tight...)
        cfd(b, h) = (e = zeros(B); e[b] = h;
                     (first(objective(m, blk, data, θ .+ e; tight...)) - first(objective(m, blk, data, θ .- e; tight...))) / (2h))
        fdv = [cfd(b, FD_H) for b in fdidx]
        fdv3 = [cfd(b, FD_H_WIDE) for b in fdidx]
        scale = maximum(abs, g.grad)
        rel = [relerr(g.grad[b], fdv[i]) for (i, b) in enumerate(fdidx)]
        rel3 = [relerr(g.grad[b], fdv3[i]) for (i, b) in enumerate(fdidx)]
        srel = [abs(g.grad[b] - fdv[i]) / scale for (i, b) in enumerate(fdidx)]
        for (i, b) in enumerate(fdidx)
            say(@sprintf("%s θ_%d: dual %.12e, FD(h=%g) %.12e rel %.2e, FD(h=%g) %.12e rel %.2e; |Δ|/max|g| %.2e (max|g| %.3e)",
                         place, b, g.grad[b], FD_H, fdv[i], rel[i], FD_H_WIDE, fdv3[i], rel3[i], srel[i], scale))
        end
        # The pass criterion is scale-relative (Tom, 2026-11-27, after smoke 2's check, job
        # 21724627, failed on one small component whose central differences scatter with
        # h): |dual − FD| ≤ FD_RTOL · max|g| over the block. The per-component errors are
        # reported beside it.
        fdres[place] = (idx = fdidx, dual = g.grad[fdidx], fd = fdv, fd_wide = fdv3, rel = rel, rel_wide = rel3,
                        scale = scale, scale_rel = srel, pass = maximum(srel) <= FD_RTOL)
    end
    res[:fd] = fdres

    # 6. DualSparseKLU against LinearSolve's own dual path, one pass, the first minute.
    libres = Dict{Symbol, Any}()
    for place in (:GE, :MET)
        blk = make_block(m, r, place, B)
        θd = seed(θstar(B), 1:CHUNK)
        run(alg) = (t = @elapsed (out_ = simulate(m, blk, θd; tend = 60.0, alg)); (out_, t))
        run(solver); run(_ -> solver_library())
        (s1, t1) = run(solver); (s2, t2) = run(_ -> solver_library())
        u1 = s1[2]; u2 = s2[2]
        num = maximum(maximum(abs, partials(u1[k]) .- partials(u2[k])) for k in eachindex(u1))
        den = maximum(maximum(abs, partials(u2[k])) for k in eachindex(u2))
        vrel = maximum(abs(value(u1[k]) - value(u2[k])) / max(abs(value(u2[k])), 1.0) for k in eachindex(u1))
        say(@sprintf("%s first minute: DualSparseKLU %.3f s (%s), library %.3f s (%s); partials %.2e, values %.2e",
                     place, t1, statline(s1[3]), t2, statline(s2[3]), num / den, vrel))
        libres[place] = (ours_s = t1, library_s = t2, ours = s1[3], library = s2[3], partials_rel = num / den,
                         value_rel = vrel, pass = num / den <= LIB_RTOL && vrel <= LIB_RTOL)
    end
    res[:library_check] = libres
    # Which linear cache the solver builds under a dual state.
    linpath = try
        blk = make_block(m, r, :GE, B); xx = copy(m.x0); θd = seed(θstar(B), 1:CHUNK)
        xd = convert(Vector{eltype(θd)}, xx)
        mnd = rebuild!(xd, m, blk, θd, ode_constants(m, blk, θd); first = true)
        fd_ = SegRHS(m, mnd)
        fun = ODEFunction{true, SciMLBase.FullSpecialize}(fd_; jac = SegJac(fd_), jac_prototype = copy(m.proto), tgrad = seg_tgrad!)
        integ = init(ODEProblem(fun, xd, (0.0, 60.0)), solver(eltype(θd)); abstol = ABSTOL, reltol = RELTOL)
        lc = integ.cache.linsolve
        string(nameof(typeof(lc)), " / ", nameof(typeof(lc.alg)))
    catch e
        "could not read: " * sprint(showerror, e)[1:min(end, 200)]
    end
    say("dual linear cache: $linpath")
    res[:linpath] = linpath

    allpass = ref_pass && rhs_pass && res[:jv].pass && all(v.pass for v in values(fdres)) &&
              all(v.pass for v in values(libres))
    serialize(joinpath(out, "check.jls"), (env = env, res = res, pass = allpass))
    open(joinpath(out, "check.md"), "w") do io
        header(io, "check", env)
        println(io, "- reference (SpeedyCell at x0): ", res[:reference])
        println(io, "- rate column: ", res[:rate_column])
        println(io, "- RHS vs BigFloat: ", res[:rhs])
        println(io, "- J·v vs FD: ", res[:jv])
        for (k, v) in fdres; println(io, "- FD gradient $k: ", v); end
        for (k, v) in libres; println(io, "- library check $k: ", v); end
        println(io, "- dual linear cache: ", linpath)
        println(io, "\nGE block (top $(maximum(BLOCKS))):"); foreach(t -> println(io, "  ", t), res[:ranking].ge)
        println(io, "\nMET block (top $(maximum(BLOCKS))):"); foreach(t -> println(io, "  ", t), res[:ranking].met)
        println(io, "\n**", allpass ? "PASS" : "FAIL", "**")
    end
    say(allpass ? "check PASS" : "check FAIL")
    allpass || error("check failed: see check.md")
end

function load_ranked(out, m)
    p = joinpath(out, "ranking.jls")
    isfile(p) && return deserialize(p)
    r = rank_blocks(m); serialize(p, (ge = r.ge, met = r.met)); return (ge = r.ge, met = r.met)
end

function stage_time(out; threaded = false)
    stage = threaded ? "threaded" : "time"
    env = environment()
    threaded && env.threads < 4 && error("threaded needs ≥ 4 threads")
    say("$stage on $(env.host), job $(env.job), $(env.threads) threads")
    m = get_model(); r = load_ranked(out, m)
    rows = Any[]
    # The primal solve.
    blk0 = make_block(m, r, :GE, 10)
    simulate(m, blk0, zeros(10); tend = T_WARM)
    tp = [(t = @elapsed (s = simulate(m, blk0, zeros(10))); (t, s[3])) for _ in 1:NREPEAT]
    primal = median(first.(tp))
    say(@sprintf("primal solve: median %.3f s, %s", primal, statline(tp[1][2])))
    nrep = threaded ? NREPEAT_THREADED : NREPEAT
    for place in PLACEMENTS, B in (threaded ? (20, 40) : BLOCKS)
        blk = make_block(m, r, place, B)
        data = make_data(m, blk)
        θ = θstar(B)
        gradient_passes(m, blk, make_data(m, blk; tend = T_WARM), θ; threaded, tend = T_WARM)   # compile
        reps = [gradient_passes(m, blk, data, θ; threaded) for _ in 1:nrep]
        walls = [x.wall for x in reps]
        pass_s = [p.seconds for x in reps for p in x.passes]
        st = reps[1].passes[1].stats
        say(@sprintf("%s B = %d: gradient median %.3f s (%s), per pass median %.3f s; pass 1 %s",
                     place, B, median(walls), join(@sprintf("%.3f", w) for w in walls), median(pass_s), statline(st)))
        push!(rows, (place = place, B = B, npass = length(passes(B)), walls = walls, median = median(walls),
                     pass_median = median(pass_s), pass_mean = sum(pass_s) / length(pass_s),
                     stats = [p.stats for p in reps[1].passes], grad = reps[1].grad))
    end
    # One row at looser tolerances, serial only: GE, B = 20, one repeat.
    loose = nothing
    if !threaded
        blk = make_block(m, r, :GE, 20); θ = θstar(20)
        tl = @elapsed (sl = simulate(m, blk, zeros(20); LOOSE...))
        data = make_data(m, blk; LOOSE...)
        g = gradient_passes(m, blk, data, θ; LOOSE...)
        loose = (primal = tl, primal_stats = sl[3], wall = g.wall, stats = [p.stats for p in g.passes])
        say(@sprintf("loose (abstol %g, reltol %g): primal %.3f s (%s); GE B = 20 gradient %.3f s, pass 1 %s",
                     LOOSE.abstol, LOOSE.reltol, tl, statline(sl[3]), g.wall, statline(g.passes[1].stats)))
    end
    serialize(joinpath(out, "$stage.jls"), (env = env, primal = primal, primal_stats = tp[1][2], rows = rows, loose = loose))
    open(joinpath(out, "$stage.md"), "w") do io
        header(io, stage, env)
        println(io, @sprintf("primal solve (GE B = 10 at θ = 0): median %.3f s; %s", primal, statline(tp[1][2])))
        println(io)
        println(io, "| placement | B | passes | gradient median (s) | repeats (s) | per pass median (s) | pass 1 stats |")
        println(io, "|---|---|---|---|---|---|---|")
        for x in rows
            println(io, @sprintf("| %s | %d | %d | %.3f | %s | %.3f | %s |", x.place, x.B, x.npass, x.median,
                                 join((@sprintf("%.3f", w) for w in x.walls), ", "), x.pass_median, statline(x.stats[1])))
        end
        if loose !== nothing
            println(io)
            println(io, @sprintf("Looser tolerances (abstol %g, reltol %g), GE B = 20, one repeat: primal %.3f s (%s); gradient %.3f s; pass 1 %s",
                                 LOOSE.abstol, LOOSE.reltol, loose.primal, statline(loose.primal_stats), loose.wall, statline(loose.stats[1])))
        end
    end
end

function stage_library(out)
    env = environment()
    say("library on $(env.host), job $(env.job)")
    m = get_model(); r = load_ranked(out, m)
    B = 20; rows = Any[]
    for place in (:GE, :MET)
        blk = make_block(m, r, place, B)
        data = make_data(m, blk)
        θd = seed(θstar(B), 1:CHUNK)
        simulate(m, blk, θd; tend = T_WARM, alg = _ -> solver_library())     # compile
        simulate(m, blk, θd; tend = T_WARM)
        tl = @elapsed (vl, sl) = objective(m, blk, data, θd; alg = _ -> solver_library())
        to = [@elapsed(objective(m, blk, data, θd)) for _ in 1:2]
        vo, so = objective(m, blk, data, θd)
        num = maximum(abs, partials(vl) .- partials(vo)) / maximum(abs, partials(vo))
        say(@sprintf("%s: library pass %.1f s (%s); DualSparseKLU median %.3f s (%s); partials rel %.2e",
                     place, tl, statline(sl), median(to), statline(so), num))
        push!(rows, (place = place, library_s = tl, library = sl, ours_s = median(to), ours_all = to, ours = so,
                     partials_rel = num))
    end
    serialize(joinpath(out, "library.jls"), (env = env, rows = rows))
    open(joinpath(out, "library.md"), "w") do io
        header(io, "library", env)
        println(io, "| placement | LinearSolve dual path, one pass (s) | its stats | DualSparseKLU, median (s) | its stats | partials rel |")
        println(io, "|---|---|---|---|---|---|")
        for x in rows
            println(io, @sprintf("| %s | %.1f | %s | %.3f | %s | %.2e |", x.place, x.library_s, statline(x.library),
                                 x.ours_s, statline(x.ours), x.partials_rel))
        end
    end
end

# The model's flux function is made by @eval in load_model, so a function compiled
# before it was made cannot call it directly (world age). main loads the model, then
# runs the stage through invokelatest.
const MODEL_REF = Ref{Any}(nothing)
get_model() = (MODEL_REF[] === nothing && (MODEL_REF[] = load_model()); MODEL_REF[]::Model)

function main()
    length(ARGS) == 2 || error("usage: route_cost_meanfield.jl <check|time|threaded|library> <out>")
    stage, out = ARGS[1], abspath(ARGS[2])
    mkpath(out)
    MODEL_REF[] = load_model()
    Base.invokelatest(run_stage, stage, out)
end

function run_stage(stage, out)
    stage == "check" ? stage_check(out) :
    stage == "time" ? stage_time(out) :
    stage == "threaded" ? stage_time(out; threaded = true) :
    stage == "library" ? stage_library(out) : error("unknown stage $stage")
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
