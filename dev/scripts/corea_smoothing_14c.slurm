#!/bin/bash
#SBATCH --job-name=infercell_smoothing_14c
#SBATCH --partition=milan
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=02:00:00
#SBATCH --array=0-49
#SBATCH --output=dev/scripts/corea_smoothing_14c_%A_%a.log

# The sampler model against the published one (spec §11 tasks 14c.5 and 14c.7):
# 100 seeds over 50 tasks, three runs per seed (published, sampler, control);
# task 0 also runs the derivative scan. Merge afterwards with
#   julia --project dev/scripts/corea_smoothing_14c_merge.jl
module load julia/1.10.5
echo "Core A′ smoothing 14c, task ${SLURM_ARRAY_TASK_ID} at $(date) on $(hostname), HEAD $(git rev-parse --short HEAD)$(test -z "$(git status --porcelain -- src test dev/scripts/corea_smoothing_14c.jl)" || echo ' (tree dirty)')"
julia --project dev/scripts/corea_smoothing_14c.jl "${SLURM_ARRAY_TASK_ID}" 50
status=$?
echo "Exit code: $status  Finished at $(date)"
exit $status
