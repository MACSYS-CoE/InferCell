#!/bin/bash
# Issue #89 B8.5: AV7's per-truth bound on the subset (spec/phases/approximate-method.md
# §12 2026-10-29, "B8 as planned").
#   dev/scripts/bound_89_b8.sh /fred/oz022/tkimpson/infercell-89-b8-bound-<date> [<r list>]
# Run it from a worktree at a pushed commit. Per injection (default 1–20): an array
# of NTASK tasks over the 1470 paired seeds (A2's split, 5 seeds a task), then the
# bound. SMOKE=1: 4 seeds, a 120 s horizon, 2 tasks. TASK_TIME sets the array's limit.
set -euo pipefail
outdir=$1; rs=${2:-$(seq -s ' ' 1 20)}
case "$outdir" in */.worktrees/*|*/.claude/worktrees/*) echo "outdir is inside a worktree" >&2; exit 1;; esac
mkdir -p "$outdir/logs"
smoke=${SMOKE:-0}
if [ "$smoke" = 1 ]; then ntask=2; task_time=01:00:00; else ntask=294; task_time=${TASK_TIME:-03:00:00}; fi
sub() { sbatch --parsable --partition=milan --exclude='dave[148-149]' --export=ALL,SMOKE="$smoke" \
               --output="$outdir/logs/%x_%A_%a.log" "$@"; }
job=dev/scripts/bound_89_b8.slurm
: >> "$outdir/jobs.tsv"
pre=$(sub --time=01:00:00 --mem=8G --job-name=bound_pre --wrap 'module load julia/1.10.5; julia --project -e "using InferCell"')
for r in $rs; do
    a=$(sub --dependency=afterok:$pre --array=0-$((ntask - 1)) --time=$task_time --job-name=bound_cells_r$r \
            $job cells "$outdir" "$r" @K "$ntask")
    b=$(sub --dependency=afterok:$a --time=02:00:00 --mem=32G --job-name=bound_r$r $job bound "$outdir" "$r")
    printf 'r%s\t%s\nr%s\t%s\n' "$r" "$a" "$r" "$b" >> "$outdir/jobs.tsv"
done
echo "precompile $pre; $(wc -l < "$outdir/jobs.tsv") jobs in $outdir/jobs.tsv"
