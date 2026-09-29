# Spec §11 task 15.7: the synthetic dataset. One truth for D11's six, drawn from
# the prior under 14c.1's Haldane-consistent rule (`draw_truth`), and 200 cells
# of the published model (clamped drain, fractional carry) at that truth, each a
# full cycle recorded every 60 s. The 200 is §4 D11's derivation: 36% relative
# SE per cell on the weakest promoter needs about 51 cells for 5%, and the
# reverse channel's detection needs 11 to 850, so 200 by default with 850 as the
# contingency.
#
# This writes the latent record only. The observation noise is applied by
# corea_dataset_15_merge.jl once 15.6 fixes the metabolite panel.
# Task 0 also reruns its first cell and asserts it is identical.
#
# Usage: julia --project dev/scripts/corea_dataset_15.jl <task> <ntasks>
# Regenerate: sbatch dev/scripts/corea_dataset_15.slurm

using InferCell
using Random
using Serialization

const TASK = parse(Int, get(ARGS, 1, get(ENV, "SLURM_ARRAY_TASK_ID", "0")))
const NTASK = parse(Int, get(ARGS, 2, "1"))
const SMOKE = get(ENV, "SMOKE", "0") == "1"
const TRUTH_SEED = 1507
const CELLS = SMOKE ? (30_001:30_002) : (30_001:30_200)
const HORIZON = SMOKE ? 120.0 : COREA_CYCLE_S
const DIR = joinpath(@__DIR__, "dataset_15")
mkpath(DIR)

models = d11_models()
truth = draw_truth(Xoshiro(TRUTH_SEED), models; purpose = :recovery)
mine = [s for (j, s) in enumerate(CELLS) if (j - 1) % NTASK == TASK]
ds = generate_dataset(truth, mine; models, horizon = HORIZON)
if TASK == 0
    again = generate_dataset(truth, mine[1:1]; models, horizon = HORIZON)
    again.latent[1, :, :] == ds.latent[1, :, :] ||
        error("two runs at seed $(mine[1]) differ: the dataset is not reproducible")
    println("seed $(mine[1]) rerun: identical")
end
serialize(joinpath(DIR, "task_$(TASK).jls"),
          (commit = strip(read(`git rev-parse --short HEAD`, String)),
           job = get(ENV, "SLURM_ARRAY_JOB_ID", get(ENV, "SLURM_JOB_ID", "none")),
           smoke = SMOKE, truth_seed = TRUTH_SEED, ds = ds))
println("Done")
