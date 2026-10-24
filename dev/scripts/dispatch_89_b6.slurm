#!/bin/bash
#SBATCH --job-name=infercell_b6_dispatch
#SBATCH --partition=milan
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=02:00:00
#SBATCH --output=test/dispatch_%j.log

module load julia/1.10.5

echo "B6 identity checks at $(date) on $(hostname)"
echo "HEAD: $(git rev-parse --short HEAD)$(test -z "$(git status --porcelain -- src test dev/scripts Project.toml Manifest.toml)" || echo ' (tree dirty)')"
INFERCELL_INTEGRATION_TESTS=true julia --project --threads=4 dev/scripts/dispatch_89_b6.jl
echo "Finished at $(date)"
