#!/bin/bash
#SBATCH --job-name=infercell_predict_ref
#SBATCH --partition=milan
#SBATCH --exclude=dave[148-149]
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=02:00:00

# Runs dev/scripts/make_predict_reference.jl, from the tree sbatch was called in,
# against each <tree> in turn on this node at 4 threads, into
# <outdir>/predict_reference_<i>.jl. Every body (all but the '#' header and the line
# naming the tree's commit) must equal the first, byte for byte.
# - One tree, twice: the output is deterministic, and _1 is the fixture.
# - The frozen tree, then a later one: the bitwise old-against-new check of the M4
#   freeze's condition (KIM-68).
# Usage, from the generator's worktree at a pushed commit:
#   sbatch --output=<outdir>/%x_%j.log dev/scripts/make_predict_reference.slurm \
#       <outdir> <tree> [<tree>...]
module load julia/1.10.5
set -euo pipefail

outdir=$(realpath -m "$1"); shift
gen="$SLURM_SUBMIT_DIR/dev/scripts/make_predict_reference.jl"
mkdir -p "$outdir"
echo "node $(hostname), job ${SLURM_JOB_ID:-none}, generator at $(git -C "$SLURM_SUBMIT_DIR" rev-parse --short HEAD), $(date)"

i=0
for tree in "$@"; do
    i=$((i + 1))
    echo "== $i: $tree at $(git -C "$tree" rev-parse --short HEAD), $(date)"
    (cd "$tree" && julia --project=. --threads=4 "$gen" "$outdir/predict_reference_$i.jl")
done

# The body less its first line, which names the tree's own commit.
body() { grep -v -e '^#' -e '^const PREDICT_REFERENCE_COMMIT = ' "$1"; }
status=0
for j in $(seq 2 "$i"); do
    if cmp <(body "$outdir/predict_reference_1.jl") <(body "$outdir/predict_reference_$j.jl"); then
        echo "IDENTICAL: body $j equals body 1"
    else
        echo "DIFFERENT: body $j differs from body 1"
        status=2
    fi
done
sha256sum "$outdir"/predict_reference_*.jl
echo "finished at $(date)"
exit $status
