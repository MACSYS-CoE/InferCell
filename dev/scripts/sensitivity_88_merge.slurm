#!/bin/bash
#SBATCH --job-name=infercell_sens88_merge
#SBATCH --partition=milan
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=32G
#SBATCH --time=02:00:00
#SBATCH --output=dev/scripts/sensitivity_88_merge_%j.log

# Issue #88 A1/A2: merge a sensitivity run. Argument: its task-file directory.
module load julia/1.10.5
echo "#88 sensitivity merge $* at $(date) on $(hostname), HEAD $(git rev-parse --short HEAD)$(test -z "$(git status --porcelain -- src dev/scripts/sensitivity_88_merge.jl)" || echo ' (tree dirty)')"
julia --project dev/scripts/sensitivity_88_merge.jl "$@"
status=$?
echo "Exit code: $status  Finished at $(date)"
exit $status
