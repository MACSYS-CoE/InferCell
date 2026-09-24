#!/bin/bash
#SBATCH --job-name=infercell_translation
#SBATCH --partition=milan
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=02:00:00
#SBATCH --output=dev/scripts/translation_diagnostics_%j.log

module load julia/1.10.5
echo "Translation diagnostics at $(date) on $(hostname), HEAD $(git rev-parse --short HEAD)$(test -z "$(git status --porcelain -- src test dev/scripts)" || echo ' (tree dirty)')"
julia --project dev/scripts/translation_diagnostics.jl
status=$?
echo "Exit code: $status  Finished at $(date)"
# Propagate it, so a driver that died on load does not look like a success.
exit $status
