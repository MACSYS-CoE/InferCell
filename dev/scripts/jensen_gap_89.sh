#!/bin/bash
# Issue #89 B3: submit every stage, chained.
#   dev/scripts/jensen_gap_89.sh /fred/oz022/tkimpson/infercell-89-b3-archive-<date>
# SMOKE=1 runs a 120 s horizon and 4 ensemble cells per M0 truth. The reduced
# model's references are #88 A3's live60 ensembles (A3_ARCHIVE), so only its mean
# field runs here. The merge waits for everything (afterany). Prints the job ids.
set -euo pipefail
outdir=$1
case "$outdir" in */.worktrees/*|*/.claude/worktrees/*) echo "outdir is inside a worktree" >&2; exit 1;; esac
mkdir -p "$outdir/logs"
smoke=${SMOKE:-0}
if [ "$smoke" = 1 ]; then nens=2; tens=00:30:00; tmf=00:30:00; else nens=10; tens=01:00:00; tmf=02:00:00; fi
sub() { sbatch --parsable --export=ALL,SMOKE="$smoke" --output="$outdir/logs/%x_%A_%a.log" "$@"; }
job=dev/scripts/jensen_gap_89.slurm
all=()
for t in m0_v6 m0_low m0_high; do
    all+=("$(sub --array=0-$((nens - 1)) --time=$tens $job ensemble $t $nens "$outdir")")
done
for t in m0_v6 m0_low m0_high corea1 corea2; do
    all+=("$(sub --time=$tmf $job meanfield $t "$outdir")")
done
deps=$(IFS=:; echo "${all[*]}")
m=$(sbatch --parsable --export=ALL,SMOKE="$smoke" --output="$outdir/logs/%x_merge_%j.log" \
    --time=00:30:00 --mem=8G --dependency=afterany:$deps $job merge "$outdir")
echo "stages ${all[*]}, merge $m"
