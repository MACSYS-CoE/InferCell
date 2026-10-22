#!/bin/bash
#SBATCH --job-name=infercell_route89
#SBATCH --partition=milan
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=02:00:00
#SBATCH --output=dev/scripts/route_ii_89_%j.log

# Issue #89 B4: one stage of route (ii)'s check against the hybrid. Submit through
# dev/scripts/route_ii_89.sh. Arguments: the stage's own.
module load julia/1.10.5
echo "#89 B4 $* at $(date) on $(hostname), HEAD $(git rev-parse --short HEAD)$(test -z "$(git status --porcelain -- src dev/scripts/route_ii_89.jl)" || echo ' (tree dirty)')"
julia --project dev/scripts/route_ii_89.jl "$@"
status=$?
echo "Exit code: $status  Finished at $(date)"
exit $status
