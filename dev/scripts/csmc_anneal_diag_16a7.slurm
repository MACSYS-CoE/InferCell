#!/bin/bash
#SBATCH --job-name=infercell_adiag
#SBATCH --partition=milan
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=06:00:00
#SBATCH --output=dev/scripts/csmc_variants_16a7/slurm_%x_%j.log

# Inside one annealed window (spec/phases/16-recovery.md task 16a.7).
# Arguments: <cell> <windows> <chains> <K...>.
module load julia/1.10.5
echo "16a.7a anneal diagnostic $* at $(date) on $(hostname), HEAD $(git rev-parse --short HEAD)$(test -z "$(git status --porcelain -- src test dev/scripts/csmc_variants_16a7.jl dev/scripts/csmc_anneal_diag_16a7.jl)" || echo ' (tree dirty)')"
julia --project dev/scripts/csmc_anneal_diag_16a7.jl "$@"
status=$?
echo "Exit code: $status  Finished at $(date)"
exit $status
