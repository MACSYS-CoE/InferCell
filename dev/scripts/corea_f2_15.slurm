#!/bin/bash
#SBATCH --job-name=infercell_f2_15
#SBATCH --partition=milan
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=00:30:00
#SBATCH --output=dev/scripts/corea_f2_15_%j.log

# Output F2 and the metabolite panel (spec §11 task 15.5).
module load julia/1.10.5
echo "Core A′ F2 15.5 at $(date) on $(hostname), HEAD $(git rev-parse --short HEAD)$(test -z "$(git status --porcelain -- src test dev/scripts/corea_f2_15.jl)" || echo ' (tree dirty)')"
julia --project dev/scripts/corea_f2_15.jl
status=$?
echo "Exit code: $status  Finished at $(date)"
exit $status
