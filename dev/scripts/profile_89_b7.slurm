#!/bin/bash
#SBATCH --job-name=infercell_prof89
#SBATCH --partition=milan
#SBATCH --exclude=dave[148-149]
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=01:30:00
#SBATCH --output=dev/scripts/profile_89_b7_%j.log

# Issue #89 B7.2: one stage of the gradient profile. Arguments: the stage's own,
# e.g. `sbatch dev/scripts/profile_89_b7.slurm snapshot <out.jls>`.
module load julia/1.10.5
echo "#89 B7.2 $* at $(date) on $(hostname), HEAD $(git rev-parse --short HEAD)$(test -z "$(git status --porcelain -- src dev/scripts/profile_89_b7.jl)" || echo ' (tree dirty)')"
julia --project dev/scripts/profile_89_b7.jl "$@"
status=$?
echo "Exit code: $status  Finished at $(date)"
exit $status
