#!/bin/bash
# Issue #88 A3: submit every stage, chained.
#   dev/scripts/eno_feedback_88.sh /fred/oz022/tkimpson/infercell-88-a3-archive-<date>
# SMOKE=1 runs a 120 s horizon, 4 ensemble cells and 2 path cells per arm.
# Per truth: the live60, nominal and live6 ensembles and paths start at once;
# popmean's wait for live60's ensemble, whose mean pools they read. The merge
# waits for everything (afterany). TPATH6 sets the live6 path tasks' time limit.
# Prints the job ids.
set -euo pipefail
outdir=$1
case "$outdir" in */.worktrees/*|*/.claude/worktrees/*) echo "outdir is inside a worktree" >&2; exit 1;; esac
mkdir -p "$outdir"
smoke=${SMOKE:-0}
if [ "$smoke" = 1 ]; then
    nens=2; ncell=2; tens=00:30:00; tpath=00:30:00; tpath6=00:30:00
else
    nens=10; ncell=40; tens=01:00:00; tpath=01:00:00; tpath6=${TPATH6:-03:00:00}
fi
sub() { sbatch --parsable --export=ALL,SMOKE="$smoke" "$@"; }
job=dev/scripts/eno_feedback_88.slurm
all=()
for t in 1 2; do
    e60=$(sub --array=0-$((nens - 1)) --time=$tens $job ensemble live60 $t $nens "$outdir")
    all+=("$e60")
    for arm in nominal live6; do
        all+=("$(sub --array=0-$((nens - 1)) --time=$tens $job ensemble $arm $t $nens "$outdir")")
    done
    all+=("$(sub --array=0-$((nens - 1)) --time=$tens --dependency=afterok:$e60 $job ensemble popmean $t $nens "$outdir")")
    for arm in live60 nominal; do
        all+=("$(sub --array=1-$ncell --time=$tpath $job path $arm $t "$outdir")")
    done
    all+=("$(sub --array=1-$ncell --time=$tpath6 $job path live6 $t "$outdir")")
    all+=("$(sub --array=1-$ncell --time=$tpath --dependency=afterok:$e60 $job path popmean $t "$outdir")")
done
deps=$(IFS=:; echo "${all[*]}")
m=$(sub --time=00:30:00 --mem=8G --dependency=afterany:$deps $job merge "$outdir")
echo "stages ${all[*]}, merge $m"
