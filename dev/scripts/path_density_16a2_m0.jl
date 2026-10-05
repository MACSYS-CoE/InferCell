# spec/phases/16-recovery.md task 16a.2, V1's importance identity on M0 at 10⁴
# paths: E_{X ~ p(·|θ0)}[p(X|θ1) / p(X|θ0)] = 1, for a CME-only shift (GAPD's
# promoter ×1.1) and an ODE-only shift (ENO's forward and reverse constants ×1.2,
# which keeps the equilibrium constant and moves the rate constants only through
# the replayed pools). Each task writes its weights; path_density_16a2_m0_merge.jl
# reports the means against 3 SE.
#
# Usage: julia --project dev/scripts/path_density_16a2_m0.jl <task> <ntasks>
# Regenerate: sbatch --array=0-99 dev/scripts/path_density_16a2_m0.slurm

using InferCell
using Random
using Serialization

const TASK = parse(Int, get(ARGS, 1, get(ENV, "SLURM_ARRAY_TASK_ID", "0")))
const NTASK = parse(Int, get(ARGS, 2, "100"))
const NPATHS = parse(Int, get(ENV, "NPATHS", "10000"))
const TRUTH_SEED = 1609
const DIR = joinpath(@__DIR__, "density_16a2_m0")
mkpath(DIR)

models = m0_models()
truth = draw_truth(Xoshiro(TRUTH_SEED), models; names = M0_TARGETS, purpose = :calibration)
θ0 = Dict(truth.values)
shifts = (cme = [:S_JCVISYN3A_0607 => 1.1 * θ0[:S_JCVISYN3A_0607]],
          ode = [:kcatF_R_ENO => 1.2 * θ0[:kcatF_R_ENO], :kcatR_R_ENO => 1.2 * θ0[:kcatR_R_ENO]])
function driver(values)
    d = build_m0()
    set_parameters!(d, models, values)
    return snapshot(d)
end
base = driver(truth.values)
alt = map(s -> driver(merge(θ0, Dict(s)) |> collect), shifts)
N = round(Int, M0_HORIZON_S)
logp(s, path) = path_logdensity(replay!(restore(s), path, N; density = true))

mine = [i for i in 1:NPATHS if (i - 1) % NTASK == TASK]
w = (cme = Float64[], ode = Float64[])
firings = Int[]
for i in mine
    d = restore(base)
    Random.seed!(160_900_000 + i)
    record_path!(d)
    for _ in 1:N
        handshake_step!(d)
    end
    path = recorded_path(d)
    push!(firings, length(path))
    l0 = logp(base, path)
    for k in keys(w)
        push!(w[k], exp(logp(alt[k], path) - l0))
    end
end
serialize(joinpath(DIR, "task_$(TASK).jls"),
          (commit = strip(read(`git rev-parse --short HEAD`, String)),
           job = get(ENV, "SLURM_ARRAY_JOB_ID", get(ENV, "SLURM_JOB_ID", "none")),
           truth = truth.values, paths = mine, w = w, firings = firings))
println("task $TASK: $(length(mine)) paths")
