# F6's second load-bearing cell: can transcripts alone shrink ENO? An upper bound
# (handoff 2026-10-03, after the bulk identifiability check left it unresolved).
#
# The transcripts are a function of a cell's full jump path X, so by the
# data-processing inequality the Fisher information they carry about the
# metabolic constants is at most the path's own. The path density p(X | θ) sees
# ENO and FBA through the rebuilt rate constants (transcription from ATP and
# GTP, translation from charged tRNA). So:
# - each task records one Core A′ cell at 15.7's truth (`Xoshiro(1507)`; seeds
#   30000 + task), and replays its path at ln θ ± δ for ENO and for FBA, at
#   δ ∈ {0.05, 0.2}, with each reverse constant derived through the nominal
#   equilibrium constant;
# - the merge forms the complete-data Fisher information of 200 cells,
#   I = 200 · E[s sᵀ] over the score s = ∂ log p(X | θ)/∂ ln θ by central
#   difference, and its observed counterpart from the curvature;
# - it reports posterior over prior SD for ENO and FBA under each target's
#   lognormal prior, with a bootstrap over cells.
# If even the full path cannot shrink ENO, transcripts cannot, and F6's cell is
# refuted. If it can, the cell remains possible, not shown.
#
# Usage: julia --project dev/scripts/eno_path_information_16a7.jl run <task>
#        julia --project dev/scripts/eno_path_information_16a7.jl merge
# Regenerate: sbatch --array=1-40 dev/scripts/eno_path_information_16a7.slurm run,
# then sbatch dev/scripts/eno_path_information_16a7.slurm merge.

using InferCell
using Random
using Serialization
using Statistics
using Printf
using LinearAlgebra
using Distributions: LogNormal, params

const DIR = joinpath(@__DIR__, "eno_path_information_16a7")
mkpath(DIR)
const FW = [:kcatF_R_ENO, :kcatF_R_FBA]
const DELTAS = [0.05, 0.2]
const NCELLS = 200
const HORIZON = COREA_CYCLE_S

function run_task(task)
    ms = d11_models()
    truth = draw_truth(Xoshiro(1507), ms; purpose = :recovery).values
    seed = 30_000 + task
    fresh() = (d = build_problem(ms; tspan = (0.0, HORIZON), complete = true);
               set_parameters!(d, ms, truth); d)
    Random.seed!(seed)
    base = snapshot(fresh())
    d = restore(base)
    Random.seed!(seed)
    record_path!(d)
    for _ in 1:round(Int, HORIZON)
        handshake_step!(d)
    end
    path = recorded_path(d)
    u0 = [log(Dict(truth)[f]) for f in FW]
    function lp(u)
        e = restore(base)
        set_parameters!(e, ms, derived_ode_values(ms, [f => exp(x) for (f, x) in zip(FW, u)]))
        return path_logdensity(replay!(e, path, round(Int, HORIZON); density = true))
    end
    lp0 = lp(u0)
    out = Dict{Tuple{Int, Float64, Int}, Float64}()
    for (i, f) in enumerate(FW), δ in DELTAS, s in (-1, 1)
        u = copy(u0)
        u[i] += s * δ
        out[(i, δ, s)] = lp(u)
        @printf("cell %d: %s %+g → Δ log p = %+.4f\n", seed, f, s * δ, out[(i, δ, s)] - lp0)
        flush(stdout)
    end
    serialize(joinpath(DIR, "cell_$(task).jls"), (; seed, lp0, out, events = length(path.times)))
end

function merge()
    files = sort(filter(f -> occursin(r"^cell_\d+\.jls$", f), readdir(DIR)))
    rs = [deserialize(joinpath(DIR, f)) for f in files]
    ms = d11_models()
    s0 = [params(InferCell._parameter_by_name(ms, f).prior)[2] for f in FW]
    @printf("# Complete-data information about ENO and FBA from one cell's full path\n\n")
    @printf("%d cells of 15.7's truth (seeds %d to %d), scaled to %d cells. Prior SD in ln θ: ENO %.3g, FBA %.3g.\n\n",
            length(rs), minimum(r.seed for r in rs), maximum(r.seed for r in rs), NCELLS, s0...)
    ratio(I) = sqrt.(diag(inv(I + Diagonal(1 ./ s0 .^ 2)))) ./ s0
    rng = Xoshiro(16_080)
    println("| δ | Score mean (ENO, FBA) | Expected-information ratio, ENO | FBA | 95% bootstrap, ENO | Observed-information ratio, ENO | FBA |")
    println("|---|---|---|---|---|---|---|")
    for δ in DELTAS
        S = reduce(hcat, [[(r.out[(i, δ, 1)] - r.out[(i, δ, -1)]) / 2δ for i in 1:2] for r in rs])'
        C = reduce(hcat, [[-(r.out[(i, δ, 1)] - 2r.lp0 + r.out[(i, δ, -1)]) / δ^2 for i in 1:2] for r in rs])'
        Iexp(rows) = NCELLS .* (S[rows, :]' * S[rows, :]) ./ length(rows)
        rE = ratio(Iexp(1:size(S, 1)))
        boot = [ratio(Iexp(rand(rng, 1:size(S, 1), size(S, 1))))[1] for _ in 1:2000]
        Iobs = Diagonal(NCELLS .* vec(mean(C; dims = 1)))
        rO = all(diag(Iobs) .> 0) ? ratio(Matrix(Iobs)) : [NaN, NaN]
        @printf("| %g | %+.3f, %+.3f | %.3f | %.3f | %.3f to %.3f | %.3f | %.3f |\n",
                δ, mean(S[:, 1]), mean(S[:, 2]), rE[1], rE[2], quantile(boot, 0.025),
                quantile(boot, 0.975), rO[1], rO[2])
    end
    println("\nA ratio near 1 is no shrinkage. The path's information bounds the transcripts' from above.")
end

if ARGS[1] == "run"
    run_task(parse(Int, ARGS[2]))
elseif ARGS[1] == "merge"
    merge()
else
    error("stage must be run or merge")
end
