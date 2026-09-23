#!/bin/bash
#SBATCH --job-name=infercell_trna_charging
#SBATCH --partition=milan
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=01:00:00
#SBATCH --output=dev/scripts/trna_charging_diagnostics_%j.log

module load julia/1.10.5
echo "tRNA charging diagnostics at $(date) on $(hostname), HEAD $(git rev-parse --short HEAD)$(test -z "$(git status --porcelain -- src test dev/scripts)" || echo ' (tree dirty)')"
julia --project dev/scripts/trna_charging_diagnostics.jl
status=$?
echo "Exit code: $status  Finished at $(date)"
# Propagate it, so a driver that died on load does not look like a success.
exit $status
