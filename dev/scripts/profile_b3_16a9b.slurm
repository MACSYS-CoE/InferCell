#!/bin/bash
#SBATCH --job-name=profile_b3
#SBATCH --partition=milan
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --mem=16G
#SBATCH --time=01:30:00
#SBATCH --output=dev/scripts/profile_b3_16a9b/slurm_%x_%j.log

# §12 2026-10-07 option 2: where block 3's time goes (profile_b3_16a9b.jl).
# Run from a worktree at a pushed commit: sbatch profile_b3_16a9b.slurm <checkpoint>
# 20 threads, as in the pilot and V6 stage 1.
module load julia/1.10.5
echo "B3 profile on $1 at $(date) on $(hostname), HEAD $(git rev-parse --short HEAD)$(test -z "$(git status --porcelain -- src dev/scripts/profile_b3_16a9b.jl)" || echo ' (tree dirty)')"
julia --project --threads=$SLURM_CPUS_PER_TASK dev/scripts/profile_b3_16a9b.jl "$1"
status=$?
echo "Exit code: $status  Finished at $(date)"
exit $status
