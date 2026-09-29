#!/bin/bash
#SBATCH --job-name=infercell_dataset_15
#SBATCH --partition=milan
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=01:00:00
#SBATCH --array=0-19
#SBATCH --output=dev/scripts/corea_dataset_15_%A_%a.log

# The synthetic dataset of spec §11 task 15.7: 200 cells over 20 tasks, at a
# drawn truth. Merge afterwards with corea_dataset_15_merge.jl.
module load julia/1.10.5
echo "Core A′ dataset 15, task ${SLURM_ARRAY_TASK_ID} at $(date) on $(hostname), HEAD $(git rev-parse --short HEAD)$(test -z "$(git status --porcelain -- src dev/scripts/corea_dataset_15.jl)" || echo ' (tree dirty)')"
julia --project dev/scripts/corea_dataset_15.jl "${SLURM_ARRAY_TASK_ID}" 20
status=$?
echo "Exit code: $status  Finished at $(date)"
exit $status
