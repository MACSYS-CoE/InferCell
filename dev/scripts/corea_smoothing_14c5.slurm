#!/bin/bash
#SBATCH --job-name=infercell_smoothing_14c5
#SBATCH --partition=milan
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=02:00:00
#SBATCH --array=0-249
#SBATCH --output=dev/scripts/corea_smoothing_14c5_%A_%a.log

# The sampler model against the published one over 5,000 seeds (spec §11 task
# 14c.5, §12 2026-09-29): 20 seeds per task, two runs per seed. Merge afterwards
# with
#   julia --project dev/scripts/corea_smoothing_14c5_merge.jl
module load julia/1.10.5
echo "Core A′ smoothing 14c.5, task ${SLURM_ARRAY_TASK_ID} at $(date) on $(hostname), HEAD $(git rev-parse --short HEAD)$(test -z "$(git status --porcelain -- src test dev/scripts/corea_smoothing_14c5.jl)" || echo ' (tree dirty)')"
julia --project dev/scripts/corea_smoothing_14c5.jl "${SLURM_ARRAY_TASK_ID}" 250
status=$?
echo "Exit code: $status  Finished at $(date)"
exit $status
