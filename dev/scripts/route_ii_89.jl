# Issue #89 B4: route (ii) against the hybrid driver at a pinned mean-field path
# (spec/phases/approximate-method.md §3 "continuous drains", §12 2026-10-23 B4).
#
# Route (ii), `build_mean_field_problem`, runs the mean field and the metabolism as
# one ODE. B3's prototype (`jensen_gap_89.jl meanfield`, reproduced here as
# `prototype`) runs the same mean field inside the hybrid driver: 1 s deferred,
# clamped drains, rounded protein counts, the rebuild reading the pools at 60m and
# growth updated every handshake. Both are deterministic, so the check has no
# Monte Carlo noise. Rule, fixed in §12 before this ran: at every truth
#   - KG(route, prototype) = Σ_pools Σ_t (Δ log x̄)² / (2σ_b²), σ_b = 10%, is at
#     most 0.1 × B3's uncorrected Jensen-gap KG total at that truth;
#   - max |Δ log x̄| over pools and saves is at most 0.02.
# x̄ is in particles, floored at one, on M0_PANEL at the 60 s saves after t = 0.
#
# Diagnostic (information only): the prototype also runs at handshake intervals of
# 0.5 and 0.25 s. The splitting reductions (deferred drains, the snapshot shift,
# per-handshake growth) shrink with the interval; protein rounding does not.
#
# Usage: julia --project dev/scripts/route_ii_89.jl prototype <truth> <interval> <outdir>
#        julia --project dev/scripts/route_ii_89.jl route <truth> <outdir>
#        julia --project dev/scripts/route_ii_89.jl merge <outdir>
# Regenerate: dev/scripts/route_ii_89.sh <outdir outside the worktrees>

using InferCell
using InferCell: AbstractSubModel
using Random
using Serialization
using Statistics
using Printf
using Dates

const SMOKE = get(ENV, "SMOKE", "0") == "1"
const COMMIT = strip(read(`git rev-parse --short HEAD`, String))
const SIGMA_B = 0.10
const B3_MERGED = "/fred/oz022/tkimpson/infercell-89-b3-archive-2026-10-20/merged.jls"
const TOL_KG_FRACTION = 0.1
const TOL_MAX_LOG = 0.02
const INTERVALS = (1.0, 0.5, 0.25)
# B3's truths (jensen_gap_89.jl).
const TRUTHS = Dict(
    "m0_v6"   => (model = :m0, seed = 16_085, horizon = M0_HORIZON_S),
    "m0_low"  => (model = :m0, seed = 16_200_013, horizon = M0_HORIZON_S),
    "m0_high" => (model = :m0, seed = 16_200_002, horizon = M0_HORIZON_S),
    "corea1"  => (model = :corea, seed = 1507, horizon = COREA_CYCLE_S),
    "corea2"  => (model = :corea, seed = 8802, horizon = COREA_CYCLE_S),
)
const ORDER = ["m0_v6", "m0_low", "m0_high", "corea1", "corea2"]
horizon(t) = SMOKE ? 120.0 : TRUTHS[t].horizon
full_models(t) = TRUTHS[t].model === :m0 ? m0_models() : d11_models()
slurm_job() = get(ENV, "SLURM_ARRAY_JOB_ID", get(ENV, "SLURM_JOB_ID", "none"))

function truth_of(t)
    spec = TRUTHS[t]
    spec.model === :m0 ?
        draw_truth(Xoshiro(spec.seed), m0_models(); names = M0_TARGETS, purpose = :calibration).values :
        draw_truth(Xoshiro(spec.seed), d11_models(); purpose = :recovery).values
end

proto_file(out, t, h) = joinpath(out, "prototype_$(t)_$(h).jls")
route_file(out, t) = joinpath(out, "route_$t.jls")

# B3's prototype, unchanged but for the handshake interval.
struct Silent{M <: AbstractSubModel} <: AbstractSubModel
    m::M
end
for f in (:states, :parameters, :formalism, :inference_mode, :coupling, :inputs,
          :written_states, :rebuilt_params, :module_id, :membrane_protein_states,
          :contributed_states)
    @eval InferCell.$f(s::Silent) = InferCell.$f(s.m)
end
InferCell.rate_constants(p, t, s::Silent, pools) = InferCell.rate_constants(p, t, s.m, pools)
InferCell.reactions(::Silent) = [InferCell.Reaction((u, p, t, w) -> 0.0, (u, w) -> nothing)]
const CME = (CoreATranscription, CoreATranscriptDecay, CoreATranslation)

function rk4(x, p, t0, h, mf; n = 4)
    f(u, t) = collect(InferCell.dynamics(u, p, t, mf))
    s = h / n
    for i in 0:(n - 1)
        t = t0 + i * s
        k1 = f(x, t); k2 = f(x .+ s / 2 .* k1, t + s / 2)
        k3 = f(x .+ s / 2 .* k2, t + s / 2); k4 = f(x .+ s .* k3, t + s)
        x = x .+ s / 6 .* (k1 .+ 2k2 .+ 2k3 .+ k4)
    end
    return x
end

function run_prototype(t, h, out)
    mkpath(out)
    full = full_models(t)
    truth = truth_of(t)
    cme = [m for m in full if any(T -> m isa T, CME)]
    ms = AbstractSubModel[any(T -> m isa T, CME) ? Silent(m) : m for m in full]
    mf = MeanField(cme; edges = :none)
    d = build_problem(ms; tspan = (0.0, horizon(t)), complete = true, rounding = :continuous,
                      interval = h)
    # One integrator lives for the whole run, and the default cap of 1e5 steps
    # counts across handshakes: corea's 6300 s at 0.25 s reached it near 6040 s.
    d.ode.opts.maxiters = 10^9
    set_parameters!(d, ms, truth)
    jnames = InferCell._block_names(ms, :jump)
    pnames = InferCell._block_param_names([m for m in ms if formalism(m) === :jump])
    sidx = [findfirst(==(s), jnames) for s in states(mf)]
    pidx = [findfirst(==(p.name), pnames) for p in InferCell.model_free_params(parameters(mf))]
    counters = Set(e.counter for m in ms for e in coupling(m) if e isa DeferredCounterEdge)
    is_counter = [s in counters for s in states(mf)]
    x = Float64.(d.jump.u[sidx])
    carry = zeros(length(x))
    rows = [findfirst(==(o), InferCell._block_names(ms, :ode)) for o in M0_PANEL]
    times, counts = Float64[], Vector{Float64}[]
    record!() = (push!(times, d.ode.t);
                 push!(counts, collect(d.ode.u[rows]) .* d.factor .+ d.rounding.remainders[rows]))
    record!()
    nsteps = round(Int, horizon(t) / h)
    per_save = round(Int, 60 / h)
    wall = @elapsed for k in 1:nsteps
        handshake_step!(d)
        x1 = rk4(x, d.jump.p[pidx], d.jump.t - h, h, mf)
        for i in eachindex(x)
            if is_counter[i]
                carry[i] += x1[i] - x[i]
                w = floor(carry[i])
                d.jump.u[sidx[i]] += Int(w)
                carry[i] -= w
            else
                d.jump.u[sidx[i]] = round(Int, x1[i])
            end
        end
        x = x1
        k % per_save == 0 && record!()
    end
    census = clipping_census(d)
    serialize(proto_file(out, t, h),
              (; truth = t, interval = h, values = truth, commit = COMMIT, smoke = SMOKE,
               job = slurm_job(), times, counts = reduce(hcat, counts), wall,
               clipped = census.clipped, drains = census.drains))
    @printf("%s prototype %s at %.2f s: %.1f s wall, %d of %d drains clipped\n",
            now(), t, h, wall, census.clipped, census.drains)
end

function run_route(t, out)
    mkpath(out)
    truth = truth_of(t)
    rt = build_mean_field_problem(full_models(t); tspan = (0.0, horizon(t)), complete = true)
    p = mean_field_parameters(rt, truth)
    saves = collect(0.0:60.0:horizon(t))
    sol = solve_mean_field(rt; p, saveat = saves)          # compiles
    timing = [@timed(solve_mean_field(rt; p, saveat = saves)) for _ in 1:3]
    sol = timing[end].value
    string(sol.retcode) == "Success" || error("route (ii) at $t: retcode $(sol.retcode)")
    rows = [findfirst(==(o), rt.state_names) for o in M0_PANEL]
    factor(u) = InferCell._route_geometry(rt, rt.n_idx == 0 ? 0.0 : u[rt.n_idx]).factor
    counts = reduce(hcat, [u[rows] .* factor(u) for u in sol.u])
    cnames = unique(b.counter for b in rt.debits)
    census = [(c, sol.u[end][findfirst(==(c), rt.state_names)]) for c in cnames]
    vol = InferCell._route_geometry(rt, rt.n_idx == 0 ? 0.0 : sol.u[end][rt.n_idx]).volume
    serialize(route_file(out, t),
              (; truth = t, values = truth, commit = COMMIT, smoke = SMOKE, job = slurm_job(),
               times = sol.t, counts, wall = [x.time for x in timing],
               bytes = [x.bytes for x in timing], census,
               growth = vol / rt.initial_volume_litres,
               labels = [l.subject for l in route_declarations(rt)]))
    @printf("%s route %s: %.2f s per solve (median of 3), %.1f MiB, growth %.4f×\n", now(), t,
            median(x.time for x in timing), median(x.bytes for x in timing) / 2^20,
            vol / rt.initial_volume_litres)
end

logs(c) = log.(max.(c, 1.0))

# Route against a prototype run, on the saves after t = 0 both have.
function diff_of(r, q)
    common = [τ for τ in r.times if τ > 0 && any(s -> isapprox(s, τ; atol = 1e-6), q.times)]
    jr = [findfirst(==(τ), r.times) for τ in common]
    jq = [findfirst(s -> isapprox(s, τ; atol = 1e-6), q.times) for τ in common]
    return logs(r.counts[:, jr]) .- logs(q.counts[:, jq])
end

kg(D) = vec(sum(D .^ 2; dims = 2)) ./ (2SIGMA_B^2)

function merge(out)
    io = IOBuffer()
    p(args...) = (println(io, args...); println(args...))
    b3 = deserialize(B3_MERGED)
    p("# Issue #89 B4: route (ii) against the hybrid at a pinned mean-field path",
      SMOKE ? " — SMOKE RUN, NOT A RESULT" : "")
    p()
    p("Merged at $COMMIT (job $(slurm_job())). σ_b = $(Int(100SIGMA_B))%. Rule (§12 2026-10-23): ",
      "KG(route, prototype at 1 s) ≤ $(TOL_KG_FRACTION) × B3's uncorrected gap total, and ",
      "max |Δ log x̄| ≤ $(TOL_MAX_LOG), at every truth.")
    p()
    p("| Truth | KG(route, proto) | B3 gap total | bound | max \\|Δ log x̄\\| (pool) | verdict | ",
      "max \\|Δ\\| at 0.5 s | at 0.25 s | route s/solve | growth |")
    p("|---|---|---|---|---|---|---|---|---|---|")
    verdicts = Bool[]
    detail = Dict{String, Any}()
    for t in ORDER
        isfile(route_file(out, t)) && isfile(proto_file(out, t, 1.0)) ||
            (p("| $t | missing | | | | | | | | |"); push!(verdicts, false); continue)
        r = deserialize(route_file(out, t))
        q = Dict(h => deserialize(proto_file(out, t, h)) for h in INTERVALS if isfile(proto_file(out, t, h)))
        r.values == q[1.0].values || error("$t: route and prototype truths differ")
        D = Dict(h => diff_of(r, q[h]) for h in keys(q))
        k = kg(D[1.0])
        gap = sum(b3[t].shift)
        mx, imx = findmax(vec(maximum(abs.(D[1.0]); dims = 2)))
        ok = SMOKE ? true : (sum(k) <= TOL_KG_FRACTION * gap && mx <= TOL_MAX_LOG)
        push!(verdicts, ok)
        at(h) = haskey(D, h) ? @sprintf("%.4f", maximum(abs.(D[h]))) : "—"
        p(@sprintf("| %s | %.3g | %.3g | %.3g | %.4f (`%s`) | %s | %s | %s | %.2f | %.4f× |", t, sum(k),
                   gap, TOL_KG_FRACTION * gap, mx, M0_PANEL[imx], ok ? "pass" : "**FAIL**", at(0.5),
                   at(0.25), median(r.wall), r.growth))
        detail[t] = (; D, k, r, q)
    end
    p()
    p("**Verdict: ", all(verdicts) ? "pass" : "FAIL", "** at ", count(verdicts), " of ", length(ORDER), " truths.")
    p()
    for t in ORDER
        haskey(detail, t) || continue
        dt = detail[t]
        p("## $t")
        p()
        p("Route at $(dt.r.commit) (job $(dt.r.job)); prototype at $(dt.q[1.0].commit) (job $(dt.q[1.0].job)), ",
          "$(dt.q[1.0].clipped) of $(dt.q[1.0].drains) drains clipped.")
        p()
        p("| Pool | KG(route, proto) | max \\|Δ\\| at 1 s | 0.5 s | 0.25 s | mean Δ at 1 s |")
        p("|---|---|---|---|---|---|")
        for (i, o) in enumerate(M0_PANEL)
            at(h) = haskey(dt.D, h) ? @sprintf("%.4f", maximum(abs.(dt.D[h][i, :]))) : "—"
            p(@sprintf("| `%s` | %.3g | %s | %s | %s | %+.4f |", o, dt.k[i], at(1.0), at(0.5), at(0.25),
                       mean(dt.D[1.0][i, :])))
        end
        p()
        p("Counters at the end (expected particles): ",
          join([@sprintf("%s %.4g", c, v) for (c, v) in dt.r.census], ", "), ".")
        p()
    end
    write(joinpath(out, "merged.md"), take!(io))
    serialize(joinpath(out, "merged.jls"), detail)
    println("Wrote ", joinpath(out, "merged.md"))
end

function main()
    stage = ARGS[1]
    if stage == "prototype"
        run_prototype(ARGS[2], parse(Float64, ARGS[3]), ARGS[4])
    elseif stage == "route"
        run_route(ARGS[2], ARGS[3])
    elseif stage == "merge"
        merge(ARGS[2])
    else
        error("stage must be prototype, route or merge")
    end
end

main()
