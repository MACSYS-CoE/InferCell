#!/bin/bash
#SBATCH --job-name=infercell_test
#SBATCH --partition=milan
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=01:00:00
#SBATCH --output=test/test_%j.log

module load julia/1.10.5

echo "Running InferCell tests at $(date)"
echo "Node: $(hostname)"
echo "HEAD: $(git rev-parse --short HEAD)$(test -z "$(git status --porcelain -- src test Project.toml Manifest.toml)" || echo ' (tree dirty)')"
echo ""

# Four threads, so the threaded-equals-serial tests (16a.9a) run threaded.
julia --project -e 'using Pkg; Pkg.test(; julia_args = ["--threads=4"])'

echo ""
echo "Finished at $(date)"
