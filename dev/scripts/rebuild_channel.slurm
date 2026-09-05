#!/bin/bash
#SBATCH --job-name=infercell_rebuild_channel
#SBATCH --partition=milan
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=02:00:00
#SBATCH --output=dev/scripts/rebuild_channel_%j.log

module load julia/1.10.5
echo "Rebuild-channel measurements at $(date) on $(hostname), HEAD $(git rev-parse --short HEAD)$(test -z "$(git status --porcelain -- src test)" || echo ' (tree dirty)')"
julia --project dev/scripts/rebuild_channel.jl
echo "Exit code: $?  Finished at $(date)"
