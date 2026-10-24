#!/bin/bash
#SBATCH --job-name=infercell_mixed89
#SBATCH --partition=milan
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=04:00:00
#SBATCH --output=dev/scripts/mixed_nuts_89_%j.log

# Issue #89 B5: one stage of the mixed likelihood's smoke. Submit through
# dev/scripts/mixed_nuts_89.sh. Arguments: the stage's own.
module load julia/1.10.5
echo "#89 B5 $* at $(date) on $(hostname), HEAD $(git rev-parse --short HEAD)$(test -z "$(git status --porcelain -- src dev/scripts/mixed_nuts_89.jl)" || echo ' (tree dirty)')"
julia --project dev/scripts/mixed_nuts_89.jl "$@"
status=$?
echo "Exit code: $status  Finished at $(date)"
exit $status
