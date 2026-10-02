#!/bin/bash
#SBATCH --job-name=infercell_rank
#SBATCH --partition=milan
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=04:00:00
#SBATCH --output=dev/scripts/csmc_variants_16a7/slurm_%x_%j.log

# Whether the true path is typical under the window proposal (spec/phases/16-recovery.md task 16a.7).
# Arguments: <scale> <cells> <proposals per window>.
module load julia/1.10.5
echo "16a.7 rank $* at $(date) on $(hostname), HEAD $(git rev-parse --short HEAD)$(test -z "$(git status --porcelain -- src test dev/scripts/csmc_variants_16a7.jl dev/scripts/csmc_rank_16a7.jl)" || echo ' (tree dirty)')"
julia --project dev/scripts/csmc_rank_16a7.jl "$@"
status=$?
echo "Exit code: $status  Finished at $(date)"
exit $status
