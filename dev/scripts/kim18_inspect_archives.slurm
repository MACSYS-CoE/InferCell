#!/bin/bash
#SBATCH --job-name=kim18_inspect
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=02:00:00
#SBATCH --partition=milan
#SBATCH --exclude=dave[148-149]

# KIM-18: load the archived V7 and B7 posteriors. Submit from the worktree root:
#   sbatch --output=<log> dev/scripts/kim18_inspect_archives.slurm
module load julia/1.10.5
echo "KIM-18 inspect at $(date) on $(hostname), HEAD $(git rev-parse --short HEAD)$(test -z "$(git status --porcelain -- src dev/scripts/kim18_inspect_archives.jl)" || echo ' (tree dirty)')"
export JULIA_PKG_OFFLINE=1
julia --project -e 'using Pkg; Pkg.instantiate()'
julia --project dev/scripts/kim18_inspect_archives.jl
status=$?
echo "Exit code: $status  Finished at $(date)"
exit $status
