#!/bin/bash
#SBATCH --job-name=infercell_testfiles
#SBATCH --partition=milan
#SBATCH --exclude=dave[148-149]
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=01:30:00

# Named test files alone, at 4 threads (so the threaded-equals-serial tests run
# threaded), after four of the suite's six test doubles (corea, contribution, jump
# and hybrid; not nucleotide and trna, which test/runtests.jl loads later): a dev
# check before the suite (KIM-68). A file that needs the other two errors here.
# Run from the tree to test.
# Usage: sbatch --output=<log> dev/scripts/testfiles.slurm test/test_a.jl test/test_b.jl ...
module load julia/1.10.5
echo "test files $* at $(date) on $(hostname), HEAD $(git rev-parse --short HEAD)$(test -z "$(git status --porcelain -- src test)" || echo ' (tree dirty)')"
julia --project --threads=4 -e 'using Test, InferCell
for f in ("corea_test_models.jl", "contribution_test_models.jl", "jump_test_models.jl",
          "hybrid_test_models.jl")
    include(joinpath("test", f))
end
@testset "test files" begin
    for f in ARGS
        include(f)
    end
end' "$@"
status=$?
echo "Exit code: $status  Finished at $(date)"
exit $status
