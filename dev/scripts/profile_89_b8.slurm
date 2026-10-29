#!/bin/bash
#SBATCH --job-name=infercell_prof89b8
#SBATCH --partition=milan
#SBATCH --exclude=dave[148-149]
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=32G
#SBATCH --time=08:00:00
#SBATCH --output=dev/scripts/profile_89_b8_%j.log

# Issue #89 B8.4: one stage of the Core A′ profile. Arguments: the stage's own,
# e.g. `sbatch dev/scripts/profile_89_b8.slurm data <out dir>`.
module load julia/1.10.5
echo "#89 B8.4 $* at $(date) on $(hostname), HEAD $(git rev-parse --short HEAD)$(test -z "$(git status --porcelain -- src dev/scripts/profile_89_b8.jl)" || echo ' (tree dirty)')"
julia --project dev/scripts/profile_89_b8.jl "$@"
status=$?
echo "Exit code: $status  Finished at $(date)"
exit $status
