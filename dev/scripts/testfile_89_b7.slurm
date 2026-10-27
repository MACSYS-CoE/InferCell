#!/bin/bash
#SBATCH --job-name=infercell_b7_tests
#SBATCH --partition=milan
#SBATCH --exclude=dave[148-149]
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=01:30:00
#SBATCH --output=dev/scripts/testfile_89_b7_%j.log

# Issue #89 B7: run named test files alone, a dev check before the suite.
# Usage: sbatch dev/scripts/testfile_89_b7.slurm test/test_mixed_likelihood.jl ...
module load julia/1.10.5
echo "#89 B7 tests $* at $(date) on $(hostname), HEAD $(git rev-parse --short HEAD)$(test -z "$(git status --porcelain -- src test)" || echo ' (tree dirty)')"
julia --project -e 'using Test; @testset "B7 dev check" begin; for f in ARGS; include(f); end; end' "$@"
status=$?
echo "Exit code: $status  Finished at $(date)"
exit $status
