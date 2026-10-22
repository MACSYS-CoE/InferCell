#!/bin/bash
# Issue #89 B4: submit every stage of route (ii)'s check against the hybrid.
#   dev/scripts/route_ii_89.sh /fred/oz022/tkimpson/infercell-89-b4-archive-<date>
# SMOKE=1 runs a 120 s horizon. One job precompiles first, on milan, since a fresh
# worktree's cache built on mixed nodes is refused there (handoff 2026-10-23). The
# merge waits for everything (afterany). Prints the job ids.
set -euo pipefail
outdir=$1
case "$outdir" in */.worktrees/*|*/.claude/worktrees/*) echo "outdir is inside a worktree" >&2; exit 1;; esac
mkdir -p "$outdir/logs"
smoke=${SMOKE:-0}
sub() { sbatch --parsable --partition=milan --exclude='dave[148-149]' \
               --export=ALL,SMOKE="$smoke" --output="$outdir/logs/%x_%j.log" "$@"; }
job=dev/scripts/route_ii_89.slurm
pre=$(sub --time=01:00:00 --wrap 'module load julia/1.10.5; julia --project -e "using InferCell"')
all=()
for t in m0_v6 m0_low m0_high corea1 corea2; do
    all+=("$(sub --dependency=afterok:$pre $job route $t "$outdir")")
    for h in 1.0 0.5 0.25; do
        all+=("$(sub --dependency=afterok:$pre $job prototype $t $h "$outdir")")
    done
done
deps=$(IFS=:; echo "${all[*]}")
m=$(sub --time=00:30:00 --mem=8G --dependency=afterany:$deps $job merge "$outdir")
echo "precompile $pre; stages ${all[*]}; merge $m"
