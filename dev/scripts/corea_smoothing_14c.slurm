#!/bin/bash
#SBATCH --job-name=infercell_smoothing_14c
#SBATCH --partition=milan
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=04:00:00
#SBATCH --array=0-9
#SBATCH --output=dev/scripts/corea_smoothing_14c_%A_%a.log

# The smoothed model against the clipped one (spec §11 task 14c.5): one seed per
# task, both drains; task 0 also runs the derivative scan. Merge afterwards with
#   julia --project dev/scripts/corea_smoothing_14c_merge.jl
module load julia/1.10.5
echo "Core A′ smoothing 14c, task ${SLURM_ARRAY_TASK_ID} at $(date) on $(hostname), HEAD $(git rev-parse --short HEAD)$(test -z "$(git status --porcelain -- src test dev/scripts/corea_smoothing_14c.jl)" || echo ' (tree dirty)')"
julia --project dev/scripts/corea_smoothing_14c.jl "${SLURM_ARRAY_TASK_ID}" 10
status=$?
echo "Exit code: $status  Finished at $(date)"
exit $status
