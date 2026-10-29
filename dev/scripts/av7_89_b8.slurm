#!/bin/bash
#SBATCH --job-name=infercell_av7
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G

# Issue #89 B8: one stage of AV7. Submit through dev/scripts/av7_89_b8.sh.
# Arguments: the stage's own; an argument @K is replaced by the array task id.
module load julia/1.10.5
args=()
for a in "$@"; do [ "$a" = @K ] && args+=("$SLURM_ARRAY_TASK_ID") || args+=("$a"); done
echo "#89 B8 ${args[*]} at $(date) on $(hostname), HEAD $(git rev-parse --short HEAD)$(test -z "$(git status --porcelain -- src dev/scripts/av7_89_b8.jl)" || echo ' (tree dirty)')"
julia --project dev/scripts/av7_89_b8.jl "${args[@]}"
status=$?
echo "Exit code: $status  Finished at $(date)"
exit $status
