#!/bin/bash
#SBATCH --job-name=infercell_full_cycle_recycling
#SBATCH --partition=milan
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=02:00:00
#SBATCH --output=dev/scripts/full_cycle_recycling_%j.log

module load julia/1.10.5
echo "Nucleotide recycling full-cycle checks at $(date) on $(hostname), HEAD $(git rev-parse --short HEAD)$(test -z "$(git status --porcelain -- src test dev/scripts)" || echo ' (tree dirty)')"
julia --project dev/scripts/full_cycle_recycling.jl
echo "Exit code: $?  Finished at $(date)"
