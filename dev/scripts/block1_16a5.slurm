#!/bin/bash
#SBATCH --job-name=infercell_v4
#SBATCH --partition=milan
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --output=dev/scripts/block1_16a5/slurm_%x_%A_%a.log

# One stage of V4 (spec/phases/16-recovery.md task 16a.5); block1_16a5.sh submits them.
# Arguments: <stage> [<n>]; an array task's id is appended as the last argument.
module load julia/1.10.5
echo "16a.5 V4 $* task ${SLURM_ARRAY_TASK_ID:-none} at $(date) on $(hostname), HEAD $(git rev-parse --short HEAD)$(test -z "$(git status --porcelain -- src test dev/scripts/block1_16a5.jl)" || echo ' (tree dirty)')"
julia --project dev/scripts/block1_16a5.jl "$@" ${SLURM_ARRAY_TASK_ID:+"$SLURM_ARRAY_TASK_ID"}
status=$?
echo "Exit code: $status  Finished at $(date)"
exit $status
