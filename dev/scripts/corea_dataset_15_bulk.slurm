#!/bin/bash
#SBATCH --job-name=infercell_ds15bulk
#SBATCH --partition=milan
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=01:00:00
#SBATCH --output=dev/scripts/corea_dataset_15_bulk_%j.log

# Parent task 15.7b: bulk observations for the 15.7 dataset (spec.md
# §12 2026-10-05). Argument: the directory of the 15.7 run's task files.
module load julia/1.10.5
echo "15.7b bulk dataset $* at $(date) on $(hostname), HEAD $(git rev-parse --short HEAD)$(test -z "$(git status --porcelain -- src dev/scripts/corea_dataset_15_bulk.jl)" || echo ' (tree dirty)')"
julia --project dev/scripts/corea_dataset_15_bulk.jl "$@"
status=$?
echo "Exit code: $status  Finished at $(date)"
exit $status
