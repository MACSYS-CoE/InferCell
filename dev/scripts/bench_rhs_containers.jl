# Spec §11 task 1.6: re-take the state-container decision at 32 states.
#
# Three right-hand sides over the same 32-state, two-module composition
# (EnergyPools + CarbonBlock from test/contribution_test_models.jl):
#
#   legacy   the pre-phase-1 `_build_rhs`, copied verbatim below. It cannot
#            execute contributions, so its derivative omits the two cross-module
#            terms; it is here to show what the framework cost before.
#   static   the phase-1 `_build_rhs`: out-of-place over an SVector, tuple of
#            modules, static index vectors, contributions folded in.
#   inplace  a hand-written in-place RHS over a Vector, the same equations and
#            the same contribution terms, as the alternative container.
#
# Reports bytes and nanoseconds per call, and wall-clock for one solve with a
# non-stiff and a stiff solver. Writes a markdown table to stdout and to
# dev/scripts/bench_rhs_containers_result.md.
#
# Usage: julia --project dev/scripts/bench_rhs_containers.jl

using InferCell
using StaticArrays
using OrdinaryDiffEq: ODEProblem, solve, Tsit5, Rodas5P
using Printf
using Dates

include(joinpath(@__DIR__, "..", "..", "test", "contribution_test_models.jl"))

# ---- the pre-phase-1 `_build_rhs`, verbatim (src/orchestrator.jl at f5fb809) ----
function legacy_build_rhs(models::Vector{<:AbstractSubModel},
                          contexts::Vector{InferCell.SubModelContext})
    model_inputs = [inputs(m) for m in models]
    function rhs(u, p, t)
        du_parts = map(models, contexts, model_inputs) do m, ctx, inp_syms
            u_local = u[ctx.state_idxs]
            p_local = p[ctx.param_idxs]
            u_inputs = [u[ctx.input_map[s]] for s in inp_syms]
            dynamics(u_local, p_local, t, m, u_inputs)
        end
        return vcat(du_parts...)
    end
    return rhs
end

# ---- the in-place alternative, same equations as the doubles ----
struct InplaceRHS
    n_energy::Int
    n_carbon::Int
    i_atp::Int
    i_pi::Int
    i_g6p::Int
end
function (r::InplaceRHS)(du, u, p, t)
    gamma_o, k_prod, k_cons, gamma_c = p[1], p[2], p[3], p[4]
    @inbounds for i in 1:r.n_energy
        du[i] = -gamma_o * u[i]
    end
    @inbounds for i in (r.n_energy + 1):(r.n_energy + r.n_carbon)
        du[i] = -gamma_c * u[i]
    end
    @inbounds du[r.i_atp] += k_prod * u[r.i_g6p]
    @inbounds du[r.i_pi] += -k_cons * u[r.i_pi]
    return nothing
end

# ---- measurement helpers: nothing captured, everything passed in ----
alloc_oop(f, u, p, t) = @allocated f(u, p, t)
alloc_ip(f, du, u, p, t) = @allocated f(du, u, p, t)
function time_oop(f, u, p, t, n)
    acc = zero(eltype(u))
    el = @elapsed for _ in 1:n
        acc += f(u, p, t)[1]
    end
    return el / n * 1e9, acc
end
function time_ip(f, du, u, p, t, n)
    acc = zero(eltype(u))
    el = @elapsed for _ in 1:n
        f(du, u, p, t); acc += du[1]
    end
    return el / n * 1e9, acc
end
function time_solve(prob, alg; reps = 5)
    solve(prob, alg; abstol = 1e-10, reltol = 1e-8, saveat = 60.0)   # warm-up
    best = Inf
    for _ in 1:reps
        best = min(best, @elapsed solve(prob, alg; abstol = 1e-10, reltol = 1e-8, saveat = 60.0))
    end
    return best
end

models = [EnergyPools(), CarbonBlock()]
prob_static = build_problem(models; tspan = (0.0, 6300.0))
u0, p0 = prob_static.u0, prob_static.p
N = length(u0)
@assert N == 32

contexts = InferCell._build_contexts(models)
InferCell._resolve_coupling(models, contexts)
rhs_legacy_raw = legacy_build_rhs(models, contexts)
# The legacy path returns whatever container the modules return. Its `u_local`
# is a UnitRange slice of the SVector, i.e. a plain Vector, so a dynamics that
# broadcasts over it returns a Vector and OrdinaryDiffEq refuses the solve
# ("non-constant types in an out-of-place ODE solve"). The five real models
# only ever worked because each returns SA[...] by hand. The conversion below
# is what a legacy module author had to do inside every dynamics.
rhs_legacy(u, p, t) = SVector{length(u)}(rhs_legacy_raw(u, p, t))
rhs_static = prob_static.f.f

gs(s) = findfirst(==(s), reduce(vcat, states.(models)))
rhs_inplace = InplaceRHS(length(ENERGY_SPECIES), length(CARBON_SPECIES),
                         gs(:M_atp_c), gs(:M_pi_c), gs(:M_g6p_c))

# The in-place RHS must agree with the static one exactly.
du_ip = zeros(N); rhs_inplace(du_ip, Vector(u0), p0, 0.0)
@assert du_ip == Vector(rhs_static(u0, p0, 0.0)) "in-place and static derivatives differ"

prob_legacy = ODEProblem{false}(rhs_legacy, u0, (0.0, 6300.0), p0)
prob_inplace = ODEProblem{true}(rhs_inplace, Vector(u0), (0.0, 6300.0), p0)

NCALLS = 200_000
# warm-ups
alloc_oop(rhs_legacy, u0, p0, 0.0); alloc_oop(rhs_static, u0, p0, 0.0); alloc_ip(rhs_inplace, du_ip, Vector(u0), p0, 0.0)
time_oop(rhs_legacy, u0, p0, 0.0, 10); time_oop(rhs_static, u0, p0, 0.0, 10); time_ip(rhs_inplace, du_ip, Vector(u0), p0, 0.0, 10)

rows = []
push!(rows, ("legacy out-of-place (no contributions; Vector->SVector conversion added)", alloc_oop(rhs_legacy, u0, p0, 0.0),
             time_oop(rhs_legacy, u0, p0, 0.0, NCALLS)[1],
             time_solve(prob_legacy, Tsit5()), time_solve(prob_legacy, Rodas5P())))
push!(rows, ("phase-1 static out-of-place", alloc_oop(rhs_static, u0, p0, 0.0),
             time_oop(rhs_static, u0, p0, 0.0, NCALLS)[1],
             time_solve(prob_static, Tsit5()), time_solve(prob_static, Rodas5P())))
uv = Vector(u0)
push!(rows, ("hand-written in-place Vector", alloc_ip(rhs_inplace, du_ip, uv, p0, 0.0),
             time_ip(rhs_inplace, du_ip, uv, p0, 0.0, NCALLS)[1],
             time_solve(prob_inplace, Tsit5()), time_solve(prob_inplace, Rodas5P())))

io = IOBuffer()
println(io, "# RHS container benchmark at 32 states (spec task 1.6)")
println(io)
println(io, "Generated by `dev/scripts/bench_rhs_containers.jl` on $(gethostname()), $(Dates.now()), Julia $(VERSION).")
println(io, "Composition: EnergyPools (13 states) + CarbonBlock (19 states), two contributions. ")
println(io, "Solves: tspan (0, 6300) s, abstol 1e-10, reltol 1e-8, saveat 60 s, best of 5.")
println(io)
println(io, "| RHS | bytes / call | ns / call | Tsit5 solve (s) | Rodas5P solve (s) |")
println(io, "|---|---|---|---|---|")
for (name, b, ns, t1, t2) in rows
    @printf(io, "| %s | %d | %.0f | %.4f | %.4f |\n", name, b, ns, t1, t2)
end
println(io)
println(io, "The legacy row omits the two cross-module terms because that RHS has no channel for them, ")
println(io, "and includes a Vector-to-SVector conversion of its result: its local state slice is a plain ")
println(io, "Vector, so a broadcasting dynamics returns one and OrdinaryDiffEq refuses the solve without it. ")
println(io, "The real models avoided that only by returning SA[...] explicitly.")
s = String(take!(io))
print(s)
write(joinpath(@__DIR__, "bench_rhs_containers_result.md"), s)
