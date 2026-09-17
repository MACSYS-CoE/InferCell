#!/bin/bash
#SBATCH --job-name=infercell_integration
#SBATCH --partition=milan
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=06:00:00
#SBATCH --output=test/integration_%j.log

# The full suite with the Bayesian integration tests switched on.
#
# These are PARKED — there is no CI workflow for them. This script is the only
# way they run, and it is kept so the work is recoverable, not because it
# currently succeeds: job 16617632 hit a 6-hour walltime on a dedicated core
# with no test output. Six hours is not enough. Raising the walltime below
# without first cutting the ABC-SMC configuration just buys a longer silence.
# See CONTRIBUTING.md, "The integration suite is parked".

module load julia/1.10.5

export INFERCELL_INTEGRATION_TESTS=true

echo "Running InferCell integration suite at $(date)"
echo "Node: $(hostname)"
echo "HEAD: $(git rev-parse --short HEAD)$(test -z "$(git status --porcelain -- src test Project.toml Manifest.toml)" || echo ' (tree dirty)')"
echo ""

julia --project -e 'using Pkg; Pkg.test()'

echo ""
echo "Finished at $(date)"
