#!/bin/bash
#SBATCH --job-name=infercell_bulk7c
#SBATCH --partition=milan
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --output=dev/scripts/csmc_bulk_16a7c/slurm_%x_%j.log

# One stage of 16a.7c's measurement (spec/phases/16-recovery.md §12 2026-10-05):
# run <N> <scans>, or merge. csmc_bulk_16a7c.sh submits them.
module load julia/1.10.5
echo "16a.7c bulk $* at $(date) on $(hostname), HEAD $(git rev-parse --short HEAD)$(test -z "$(git status --porcelain -- src dev/scripts/csmc_bulk_16a7c.jl)" || echo ' (tree dirty)')"
julia --project dev/scripts/csmc_bulk_16a7c.jl "$@"
status=$?
echo "Exit code: $status  Finished at $(date)"
exit $status
