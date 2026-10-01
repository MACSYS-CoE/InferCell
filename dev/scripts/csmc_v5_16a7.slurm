#!/bin/bash
#SBATCH --job-name=infercell_v5
#SBATCH --partition=milan
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --output=dev/scripts/csmc_v5_16a7/slurm_%x_%A_%a.log

# One stage of V5 (spec/phases/16-recovery.md task 16a.7); csmc_v5_16a7.sh submits them.
module load julia/1.10.5
echo "16a.7 V5 $* task ${SLURM_ARRAY_TASK_ID:-none} at $(date) on $(hostname), HEAD $(git rev-parse --short HEAD)$(test -z "$(git status --porcelain -- src test dev/scripts/csmc_v5_16a7.jl)" || echo ' (tree dirty)')"
julia --project dev/scripts/csmc_v5_16a7.jl "$@" ${SLURM_ARRAY_TASK_ID:+"$SLURM_ARRAY_TASK_ID"}
status=$?
echo "Exit code: $status  Finished at $(date)"
exit $status
