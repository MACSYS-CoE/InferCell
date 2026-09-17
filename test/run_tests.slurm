#!/bin/bash
#SBATCH --job-name=infercell_test
#SBATCH --partition=milan
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=00:30:00
#SBATCH --output=test/test_%j.log

module load julia/1.10.5

echo "Running InferCell tests at $(date)"
echo "Node: $(hostname)"
echo "HEAD: $(git rev-parse --short HEAD)$(test -z "$(git status --porcelain -- src test Project.toml Manifest.toml)" || echo ' (tree dirty)')"
echo ""

julia --project -e 'using Pkg; Pkg.test()'

echo ""
echo "Finished at $(date)"
