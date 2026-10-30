#!/bin/bash
# KIM-66's done-when: the variance-split estimator checked by simulation on known
# splits (dev/scripts/variance_split_check_kim66.jl), from a worktree at a pushed
# commit.
#   dev/scripts/variance_split_check_kim66.sh /fred/oz022/tkimpson/infercell-m4-kim66-<date>
# One job; its id goes to <out>/jobs.tsv, its log to <out>/logs.
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
j=$(sbatch --parsable --partition=milan --exclude='dave[148-149]' --ntasks=1 \
           --cpus-per-task=1 --mem=8G --time=02:00:00 --job-name=kim66_check \
           --output="$outdir/logs/%x_%j.log" \
           --wrap "module load julia/1.10.5; julia --project dev/scripts/variance_split_check_kim66.jl $outdir")
printf 'check\t%s\n' "$j" > "$outdir/jobs.tsv"
echo "submitted at $(git rev-parse --short HEAD):"
cat "$outdir/jobs.tsv"
