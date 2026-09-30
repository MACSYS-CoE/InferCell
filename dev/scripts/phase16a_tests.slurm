#!/bin/bash
#SBATCH --job-name=infercell_16a_tests
#SBATCH --partition=milan
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=00:45:00
#SBATCH --output=dev/scripts/phase16a_tests_%j.log

# Only phase 16a's test files, with the doubles they need, to iterate before
# the full suite.
module load julia/1.10.5
echo "16a tests at $(date) on $(hostname), HEAD $(git rev-parse --short HEAD)$(test -z "$(git status --porcelain -- src test)" || echo ' (tree dirty)')"
julia --project -e 'using Test, InferCell
for f in ("corea_test_models.jl", "contribution_test_models.jl", "jump_test_models.jl",
          "hybrid_test_models.jl", "test_path_replay.jl", "test_path_density.jl")
    include(joinpath("test", f))
end'
status=$?
echo "Exit code: $status  Finished at $(date)"
exit $status
