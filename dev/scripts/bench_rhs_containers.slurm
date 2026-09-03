#!/bin/bash
#SBATCH --job-name=infercell_bench_rhs
#SBATCH --partition=milan
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=00:30:00
#SBATCH --output=dev/scripts/bench_rhs_%j.log

module load julia/1.10.5
echo "RHS container benchmark at $(date) on $(hostname), HEAD $(git rev-parse --short HEAD)$(test -z "$(git status --porcelain -- src test)" || echo ' (tree dirty)')"
julia --project dev/scripts/bench_rhs_containers.jl
echo "Exit code: $?  Finished at $(date)"
