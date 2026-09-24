#!/bin/bash
#SBATCH --job-name=infercell_drain_granularity
#SBATCH --partition=milan
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=08:00:00
#SBATCH --output=dev/scripts/drain_granularity_%j.log

module load julia/1.10.5
echo "Drain granularity at $(date) on $(hostname), HEAD $(git rev-parse --short HEAD)$(test -z "$(git status --porcelain -- src test dev/scripts)" || echo ' (tree dirty)')"
echo "SEEDS=${SEEDS:-50}"
julia --project dev/scripts/drain_granularity.jl
status=$?
echo "Exit code: $status  Finished at $(date)"
# Propagate it, so a driver that died on load does not look like a success.
exit $status
