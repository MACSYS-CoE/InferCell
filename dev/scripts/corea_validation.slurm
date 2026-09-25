#!/bin/bash
#SBATCH --job-name=infercell_corea_validation
#SBATCH --partition=milan
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=03:00:00
#SBATCH --output=dev/scripts/corea_validation_%j.log

module load julia/1.10.5
echo "Core A′ phase 14a validation at $(date) on $(hostname), HEAD $(git rev-parse --short HEAD)$(test -z "$(git status --porcelain -- src test dev/scripts)" || echo ' (tree dirty)')"
echo "SMOKE=${SMOKE:-0}"
julia --project dev/scripts/corea_validation.jl
status=$?
echo "Exit code: $status  Finished at $(date)"
# Propagate it, so a driver that died on load does not look like a success.
exit $status
