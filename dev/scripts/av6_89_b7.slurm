#!/bin/bash
#SBATCH --job-name=infercell_av6
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G

# Issue #89 B7: one stage of AV6. Submit through dev/scripts/av6_89_b7.sh.
# Arguments: the stage's own.
module load julia/1.10.5
echo "#89 B7 $* at $(date) on $(hostname), HEAD $(git rev-parse --short HEAD)$(test -z "$(git status --porcelain -- src dev/scripts/av6_89_b7.jl)" || echo ' (tree dirty)')"
julia --project dev/scripts/av6_89_b7.jl "$@"
status=$?
echo "Exit code: $status  Finished at $(date)"
exit $status
