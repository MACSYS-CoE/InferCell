#!/bin/bash
#SBATCH --job-name=v7_sbc
#SBATCH --partition=milan
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=5
#SBATCH --mem=8G
#SBATCH --output=dev/scripts/sbc_v7_16a9b/slurm_%x_%A_%a.log

# One stage of 16a.9b's V7 (spec/phases/16-recovery.md §12 2026-10-08):
# as an array task, replicate $SLURM_ARRAY_TASK_ID; otherwise the given stage
# (merge). sbc_v7_16a9b.sh submits them. 5 threads: a sweep costs 0.78 CPU-h
# here, against 1.22 at 20 (thread_scaling_16a9b_result.md).
module load julia/1.10.5
if [ -n "$SLURM_ARRAY_TASK_ID" ]; then set -- run "$SLURM_ARRAY_TASK_ID"; fi
echo "16a.9b V7 $* at $(date) on $(hostname), HEAD $(git rev-parse --short HEAD)$(test -z "$(git status --porcelain -- src dev/scripts/sbc_v7_16a9b.jl)" || echo ' (tree dirty)')"
julia --project --threads=$SLURM_CPUS_PER_TASK dev/scripts/sbc_v7_16a9b.jl "$@"
status=$?
echo "Exit code: $status  Finished at $(date)"
exit $status
