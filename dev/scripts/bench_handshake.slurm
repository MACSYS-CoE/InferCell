#!/bin/bash
#SBATCH --job-name=infercell_bench_handshake
#SBATCH --partition=milan
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=01:00:00
#SBATCH --output=dev/scripts/bench_handshake_%j.log

module load julia/1.10.5
echo "Handshake wall-clock benchmark at $(date) on $(hostname), HEAD $(git rev-parse --short HEAD)$(test -z "$(git status --porcelain -- src test)" || echo ' (tree dirty)')"
julia --project dev/scripts/bench_handshake.jl
echo "Exit code: $?  Finished at $(date)"
