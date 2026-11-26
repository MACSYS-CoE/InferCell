#!/bin/bash
#SBATCH --job-name=kim90_reads
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=03:00:00
#SBATCH --partition=milan

# KIM-90 items 1 and 2 in one job, so the note's numbers carry a job id:
# the dataset fingerprint at V6's commit (281d548), at B7's (7a102fc) and at
# this commit, the identity comparison, and the Monte Carlo error table.
#   sbatch --output=<out>/logs/%x_%j.log dev/scripts/gap_diag_reads.slurm <out>
# Run from a worktree at a pushed commit; <out> outside the worktrees.
set -uo pipefail
export EBROOTJULIA=/apps/modules/software/Julia/1.10.5
export PATH="$EBROOTJULIA/bin:$PATH"
export JULIA_PKG_OFFLINE=1
out=$1
mkdir -p "$out"
here=$(pwd)
wt=/fred/oz022/tkimpson/InferCell/.worktrees
echo "KIM-90 reads at $(date) on $(hostname), HEAD $(git rev-parse --short HEAD)"
s=0
for c in 281d548:v6-16a9b 7a102fc:89b7-7a102fc; do
    sha=${c%%:*}; dir=$wt/${c#*:}
    [ "$(git -C "$dir" rev-parse --short HEAD)" = "$sha" ] || { echo "$dir is not at $sha" >&2; exit 1; }
    julia --project="$dir" "$here/dev/scripts/gap_diag_identity.jl" fingerprint "$out/fp_$sha.tsv" || s=1
done
julia --project dev/scripts/gap_diag_identity.jl fingerprint "$out/fp_$(git rev-parse --short HEAD).tsv" || s=1
julia --project dev/scripts/gap_diag_identity.jl compare "$out"/fp_*.tsv | tee "$out/identity.txt"
s=$(( s | PIPESTATUS[0] ))
julia --project dev/scripts/gap_diag_mixing.jl "$out/mixing.md" || s=1
echo "Exit code: $s  Finished at $(date)"
exit $s
