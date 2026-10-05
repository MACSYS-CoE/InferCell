#!/bin/bash
#SBATCH --job-name=infercell_kink_16a6
#SBATCH --partition=milan
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=03:00:00
#SBATCH --output=dev/scripts/fixed_path_kink_16a6_%j.log

# A clipped drain under a fixed path, on Core A′ (spec/phases/16-recovery.md task 16a.6).
module load julia/1.10.5
echo "16a.6 fixed-path kink at $(date) on $(hostname), HEAD $(git rev-parse --short HEAD)$(test -z "$(git status --porcelain -- src test dev/scripts/fixed_path_kink_16a6.jl)" || echo ' (tree dirty)')"
julia --project dev/scripts/fixed_path_kink_16a6.jl
status=$?
echo "Exit code: $status  Finished at $(date)"
exit $status
