#!/bin/bash
#SBATCH --job-name=infercell_replay_16a1
#SBATCH --partition=milan
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=01:00:00
#SBATCH --output=dev/scripts/path_replay_16a1_%j.log

# Path replay at full scale (spec/phases/16-recovery.md task 16a.1, V2).
module load julia/1.10.5
echo "16a.1 path replay at $(date) on $(hostname), HEAD $(git rev-parse --short HEAD)$(test -z "$(git status --porcelain -- src test dev/scripts/path_replay_16a1.jl)" || echo ' (tree dirty)')"
julia --project dev/scripts/path_replay_16a1.jl
status=$?
echo "Exit code: $status  Finished at $(date)"
exit $status
