#!/bin/bash
#SBATCH --job-name=infercell_bulkid
#SBATCH --partition=milan
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=01:00:00
#SBATCH --output=dev/scripts/bulk_identifiability_16a7_%j.log

# Bulk identifiability from 15.8's sensitivity run (spec/phases/16-recovery.md
# task 16a.7, handoff 2026-10-03). Argument: the directory of 15.8's task files.
module load julia/1.10.5
echo "16a.7 bulk identifiability $* at $(date) on $(hostname), HEAD $(git rev-parse --short HEAD)$(test -z "$(git status --porcelain -- src dev/scripts/bulk_identifiability_16a7.jl)" || echo ' (tree dirty)')"
julia --project dev/scripts/bulk_identifiability_16a7.jl "$@"
status=$?
echo "Exit code: $status  Finished at $(date)"
exit $status
