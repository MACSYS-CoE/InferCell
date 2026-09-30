#!/bin/bash
#SBATCH --job-name=infercell_replay_tests
#SBATCH --partition=milan
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=00:45:00
#SBATCH --output=dev/scripts/path_replay_16a1_tests_%j.log

# Only test/test_path_replay.jl, to iterate on 16a.1 before the full suite.
module load julia/1.10.5
echo "16a.1 tests at $(date) on $(hostname), HEAD $(git rev-parse --short HEAD)$(test -z "$(git status --porcelain -- src test)" || echo ' (tree dirty)')"
julia --project -e 'using Test, InferCell; include("test/test_path_replay.jl")'
status=$?
echo "Exit code: $status  Finished at $(date)"
exit $status
