#!/bin/bash
# KIM-27's done-when: a smoke run of every stage of `predict`
# (dev/scripts/predict_smoke_kim27.jl), from a worktree at a pushed commit.
#   dev/scripts/predict_smoke_kim27.sh /fred/oz022/tkimpson/infercell-m3-smoke-<date>
# One job precompiles first, on milan; then one job per stage. Job ids go to
# <out>/jobs.tsv, logs to <out>/logs.
set -euo pipefail
outdir=$(realpath -m "$1")
case "$outdir" in "$(git rev-parse --show-toplevel)"/*|*/.worktrees/*|*/.claude/worktrees/*)
    echo "outdir $outdir is inside a worktree" >&2; exit 1;; esac
if [ -n "$(git status --porcelain -- src test dev/scripts Project.toml Manifest.toml)" ]; then
    echo "the tree is dirty; run from a clean worktree at a pushed commit" >&2; exit 1
fi
if [ -z "$(git branch -r --contains HEAD)" ]; then
    echo "HEAD $(git rev-parse --short HEAD) is on no remote branch; push it first" >&2; exit 1
fi
mkdir -p "$outdir/logs"
sub() { sbatch --parsable --partition=milan --exclude='dave[148-149]' --ntasks=1 \
               --cpus-per-task=1 --mem=16G --output="$outdir/logs/%x_%j.log" "$@"; }
: > "$outdir/jobs.tsv"
pre=$(sub --time=01:00:00 --job-name=kim27_pre --wrap 'module load julia/1.10.5; julia --project -e "using InferCell"')
printf 'precompile\t%s\n' "$pre" >> "$outdir/jobs.tsv"
for stage in quickstart ingest m0 intervention reduced wrapper; do
    j=$(sub --dependency=afterok:$pre --time=02:00:00 --job-name=kim27_$stage \
            --wrap "module load julia/1.10.5; julia --project dev/scripts/predict_smoke_kim27.jl $stage $outdir")
    printf '%s\t%s\n' "$stage" "$j" >> "$outdir/jobs.tsv"
done
echo "submitted at $(git rev-parse --short HEAD):"
cat "$outdir/jobs.tsv"
