#!/bin/bash
#SBATCH --job-name=tscale
#SBATCH --partition=milan
#SBATCH --ntasks=1
#SBATCH --mem=16G
#SBATCH --output=dev/scripts/thread_scaling_16a9b/slurm_%x_%j.log

# §12 2026-10-07 option 2: a sweep's cost per thread count (thread_scaling_16a9b.jl).
# Run from a worktree at a pushed commit, once per count:
#   sbatch --cpus-per-task=<n> --time=<t> thread_scaling_16a9b.slurm <checkpoint> <nsweeps>
module load julia/1.10.5
echo "Thread scaling, $SLURM_CPUS_PER_TASK threads, on $1 at $(date) on $(hostname), HEAD $(git rev-parse --short HEAD)$(test -z "$(git status --porcelain -- src dev/scripts/thread_scaling_16a9b.jl)" || echo ' (tree dirty)')"
julia --project --threads=$SLURM_CPUS_PER_TASK dev/scripts/thread_scaling_16a9b.jl "$1" "$2"
status=$?
echo "Exit code: $status  Finished at $(date)"
exit $status
