# KIM-90 items 3, 4a and 4c: refit m0_v6 with the mixed route's bulk reference changed
# one way at a time, and measure each refit's shift against V6 and against an
# unchanged control.
#
# Every arm reads B7's archived prep for m0_v6 (dataset, θ_ref, the 800 reference
# seeds and cells, δ, S; job 19177809 at 7a102fc) and runs AV6's design: 2 NUTS
# chains of 200 + 200 at target acceptance 0.8, on v6's chain seeds
# (mixed_chain_seeds(ds, 2; arm = 0)) and jittered starts around θ_ref (0.25 prior
# SDs). Arms differ only in (δ, S):
#   control   the archived δ and S: m0_v6 rerun on this branch's code, a check that
#             it reproduces B7 (same seeds and starts, so not a noise floor);
#   n1600     δ and S from 1600 cells at θ_ref: the archived 800 plus the next 800
#             of reference_seeds(ds, 1600), whose first 800 are the archived ones;
#   nojensen  δ = 0, S held: the Jensen correction off;
#   postmean  δ and S from the archived 800 seeds at m0_v6's posterior mean instead
#             of θ_ref (item 4c, δ's θ-dependence; approved 2026-11-19, see KIM-90).
#
# Stages:
#   prep <out> <arm>        (δ, S) for the arm, to <out>/<arm>/prep.jls;
#   chain <out> <arm> <c>   one chain, to <out>/<arm>/chain_<c>.jls;
#   merge <out>             <out>/refit.md: each arm against V6, m0_v6 and control.
#
# Usage: julia --project dev/scripts/gap_diag_refit.jl <stage> <out> [<arm> [<c>]]
# Submit with dev/scripts/gap_diag_refit.sh <out outside the worktrees> <arms…>.
# SMOKE=1: 8 new reference cells, 10 + 10 NUTS steps (the dataset is the archive's).

using InferCell
using InferCell: _mixed_solve
using Random
using Serialization
using Statistics
using LinearAlgebra
using Printf
using Dates
using Turing: sample, NUTS
using MCMCChains

const SMOKE = get(ENV, "SMOKE", "0") == "1"
const COMMIT = strip(read(`git rev-parse --short HEAD`, String))
const NADAPT = SMOKE ? 10 : 200
const NDRAW = SMOKE ? 10 : 200
const NCHAIN = 2
const TARGET_ACCEPT = 0.8
const JITTER_SD = 0.25
const NCELLS = 50                                   # m0_v6's cells; δ_and_S's C
const NEXTRA = SMOKE ? 8 : 800
const ARMS = (:control, :n1600, :nojensen, :postmean)
const B7_RUN = "/fred/oz022/tkimpson/infercell-89-b7-archive-2026-10-27/run/v6"
const F10_DRAWS = joinpath(@__DIR__, "f10_16a10_draws.tsv")
slurm_job() = get(ENV, "SLURM_ARRAY_JOB_ID", get(ENV, "SLURM_JOB_ID", "none"))

armdir(out, arm) = joinpath(out, String(arm))
b7prep() = deserialize(joinpath(B7_RUN, "prep.jls"))
b7draws() = cat((Array(deserialize(joinpath(B7_RUN, "chain_$c.jls")).chain[:, [Symbol("x[$i]") for i in 1:6], 1])
                 for c in 1:2)...; dims = 3)              # draws × 6 × chains

function run_prep(out, arm)
    arm in ARMS || error("unknown arm $arm")
    mkpath(armdir(out, arm))
    p = b7prep()
    s0 = mixed_setup(p.ds)
    t = @elapsed if arm === :control
        δ, S, info = p.delta, p.S, "archived δ and S (800 cells at θ_ref)"
    elseif arm === :n1600
        seeds = reference_seeds(p.ds, 800 + NEXTRA)
        seeds[1:800] == p.seeds || error("reference_seeds(ds, 1600)'s first 800 are not the archive's")
        extra = reference_ensemble(p.values, seeds[801:end])
        cells = vcat(p.cells, extra)
        d = delta_and_S(cells, p.m_ref; C = NCELLS)
        δ, S, info = d.delta, d.S, "$(size(cells, 1)) cells at θ_ref"
    elseif arm === :nojensen
        δ, S, info = zero(p.delta), p.S, "δ = 0, archived S"
    else # :postmean
        x̄ = vec(mean(b7draws(); dims = (1, 3)))
        values = mixed_theta_values(s0, x̄)
        m_at = mean_field_log_bulk(s0, _mixed_solve(s0, exp.(x̄[1:length(s0.names)])).sol)
        cells = reference_ensemble(values, p.seeds)
        d = delta_and_S(cells, m_at; C = NCELLS)
        δ, S, info = d.delta, d.S, "800 archived seeds at the posterior mean $(round.(x̄; digits = 4))"
    end
    serialize(joinpath(armdir(out, arm), "prep.jls"),
              (; arm, delta = δ, S, info, commit = COMMIT, smoke = SMOKE, job = slurm_job(), seconds = t))
    @printf("%s prep %s: %s; ‖δ − δ_b7‖ = %.4g, tr S = %.4g (b7 %.4g), in %.0f s\n", now(), arm, info,
            norm(δ .- p.delta), tr(S), tr(p.S), t)
end

function run_chain(out, arm, c)
    p = b7prep()
    a = deserialize(joinpath(armdir(out, arm), "prep.jls"))
    s = with_reference(mixed_setup(p.ds), a.delta, a.S)
    seed = mixed_chain_seeds(p.ds, NCHAIN; arm = 0)[c]
    rng = Xoshiro(seed)
    st = mixed_jittered_start(rng, s, p.ref.x; scale = JITTER_SD)
    t = @elapsed chain = sample(rng, mixed_model(s), NUTS(NADAPT, TARGET_ACCEPT), NDRAW;
                                initial_params = st.x, progress = false)
    serialize(joinpath(armdir(out, arm), "chain_$c.jls"),
              (; arm, c, seed, start = st.x, chain, seconds = t, nadapt = NADAPT, ndraw = NDRAW,
               commit = COMMIT, job = slurm_job()))
    @printf("%s chain %s/%d: seed %d, %d + %d steps in %.0f s\n", now(), arm, c, seed, NADAPT, NDRAW, t)
end

function w1(a::AbstractVector, b::AbstractVector)
    xs = sort(vcat(a, b)); sa, sb = sort(a), sort(b); s = 0.0
    for i in 1:(length(xs) - 1)
        s += abs(searchsortedlast(sa, xs[i]) / length(sa) - searchsortedlast(sb, xs[i]) / length(sb)) *
             (xs[i + 1] - xs[i])
    end
    return s
end

function v6_draws()
    lines = readlines(F10_DRAWS)
    head = split(lines[1], '\t')
    names_ = ["S_JCVISYN3A_0607", "S_JCVISYN3A_0779", "krnadeg", "kcatF_R_ENO", "kcatF_R_FBA", "sigma_b"]
    cols = [findfirst(==(n), head) for n in names_]
    rows = [parse.(Float64, split(l, '\t')) for l in lines[2:end]]
    return reduce(hcat, [r[cols] for r in rows if r[1] in 1:4 && r[2] >= 101])'   # draws × 6
end

function run_merge(out)
    v6 = v6_draws(); sd = vec(std(v6; dims = 1)); μ6 = vec(mean(v6; dims = 1))
    b7 = b7draws()
    pooled(X) = reduce(vcat, [X[:, :, k] for k in axes(X, 3)])
    arms = [a for a in ARMS if isfile(joinpath(armdir(out, a), "chain_2.jls"))]
    X = Dict(a => cat((Array(deserialize(joinpath(armdir(out, a), "chain_$c.jls")).chain[:, [Symbol("x[$i]") for i in 1:6], 1])
                       for c in 1:2)...; dims = 3) for a in arms)
    io = IOBuffer(); say(s) = (println(s); println(io, s))
    say("# KIM-90: refits of m0_v6 with the bulk reference changed\n")
    for a in arms
        pr = deserialize(joinpath(armdir(out, a), "prep.jls"))
        jobs = [deserialize(joinpath(armdir(out, a), "chain_$c.jls")).job for c in 1:2]
        say("- **$a**: $(pr.info). prep job $(pr.job), chain jobs $(join(jobs, ", ")), commit $(pr.commit).")
    end
    say("\nUnits: V6 posterior SDs. shift = mean(arm) − mean(V6); Δ vs control = mean(arm) − mean(control).\n")
    pn = ["S_0607", "S_0779", "krnadeg", "ENO", "FBA", "σ_b"]
    say("| Arm | " * join(["$n W1 / shift / Δ ctrl" for n in pn], " | ") * " | max R̂ |")
    say("|---|" * repeat("---|", length(pn) + 1))
    ctrl = haskey(X, :control) ? pooled(X[:control]) : nothing
    for (lab, D) in vcat([("m0_v6 (B7 archive)", b7)], [(String(a), X[a]) for a in arms])
        P = pooled(D)
        cells = map(1:6) do q
            Δc = ctrl === nothing ? NaN : (mean(P[:, q]) - mean(ctrl[:, q])) / sd[q]
            @sprintf("%.3f / %+.3f / %+.3f", w1(P[:, q], v6[:, q]) / sd[q], (mean(P[:, q]) - μ6[q]) / sd[q], Δc)
        end
        rh = maximum(MCMCChains.rhat(Chains(permutedims(D, (1, 2, 3)))).nt.rhat)
        say("| $lab | " * join(cells, " | ") * @sprintf(" | %.3f |", rh))
    end
    write(joinpath(out, "refit.md"), take!(io))
end

if abspath(PROGRAM_FILE) == (@__FILE__)
    stage, out = ARGS[1], ARGS[2]
    stage == "prep" ? run_prep(out, Symbol(ARGS[3])) :
    stage == "chain" ? run_chain(out, Symbol(ARGS[3]), parse(Int, ARGS[4])) :
    stage == "merge" ? run_merge(out) : error("unknown stage $stage")
end
