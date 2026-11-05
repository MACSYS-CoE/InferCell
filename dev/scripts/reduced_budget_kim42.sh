#!/bin/bash
# KIM-42: the reduced model's uncertainty budget, unscored
# (dev/scripts/reduced_budget_kim42.jl). It runs the tree it lives in: keep a detached
# worktree at a pushed commit for the run, call this launcher from it, and leave it
# untouched until the run ends.
#   dev/scripts/reduced_budget_kim42.sh <outdir>
#       the whole run: precompile, one job per width on its own node, then the report;
#   dev/scripts/reduced_budget_kim42.sh <outdir> missing
#       once none of the run's jobs is queued, from a tree at the run's commit and with
#       its SMOKE: precompiles, then reruns the widths with no output, then the report;
#   dev/scripts/reduced_budget_kim42.sh <outdir> progress
#       lists what is written, and the run's jobs in the queue;
#   dev/scripts/reduced_budget_kim42.sh <outdir> secheck
#       the SE method's check only: precompile and one single-core job.
# SMOKE=1 in the launcher's environment runs the smoke pass, and is passed to the jobs.
# THREADS (default 64, smoke 8), MEM (192G, smoke 32G) and WALLTIME (24:00:00, smoke
# 01:00:00) set each width's job. --gcthreads=N,1 is KIM-68's: one concurrent sweep thread.
# Job ids go to <outdir>/jobs.tsv, logs to <outdir>/logs.
set -euo pipefail
outdir=$(realpath -m "$1")
mode=${2:-run}
smoke=${SMOKE:-0}
if [ "$smoke" = 1 ]; then
    name=kim42smoke; threads=${THREADS:-8}; mem=${MEM:-32G}; wall=${WALLTIME:-01:00:00}
else
    name=kim42; threads=${THREADS:-64}; mem=${MEM:-192G}; wall=${WALLTIME:-24:00:00}
fi
jobs="${name}_pre,${name}_width,${name}_report,${name}_secheck"

if [ "$mode" = progress ]; then
    ls "$outdir"
    squeue -h -u "$USER" -n "$jobs" -o "%.14i %.18j %.8T %.10M %R" || true
    exit 0
fi
case "$mode" in run|missing|secheck) ;; *) echo "unknown mode $mode" >&2; exit 1;; esac

# The tree this launcher lives in, wherever it is called from.
cd "$(dirname "$(realpath "$0")")/../.."
source dev/scripts/launcher_guard.sh
script=dev/scripts/reduced_budget_kim42.jl
exp="--export=ALL,SMOKE=$smoke"
arr() { sbatch --parsable "${nodes[@]}" --ntasks=1 --output="$outdir/logs/%x_%A_%a.log" "$@"; }

head=$(git rev-parse HEAD)
widths=0-2
if [ "$mode" = missing ]; then
    queued=$(squeue -h -u "$USER" -n "$jobs") ||
        { echo "squeue failed, so whether the run has jobs queued is unknown" >&2; exit 1; }
    [ -z "$queued" ] || { echo "the run still has jobs in the queue; wait for them" >&2; exit 1; }
    ran() { awk -F'\t' -v k="$1" '$1 == k { print $2; exit }' "$outdir/jobs.tsv"; }
    [ "$(ran commit)" = "$head" ] ||
        { echo "the run is at '$(ran commit)', this tree at $head" >&2; exit 1; }
    [ "$(ran smoke)" = "$smoke" ] ||
        { echo "the run has SMOKE='$(ran smoke)', this launcher SMOKE=$smoke" >&2; exit 1; }
    widths=""
    for i in 0 1 2; do
        if [ ! -f "$outdir/prediction_w$i.jls" ]; then
            for f in w$i plain_w$i; do
                [ ! -f "$outdir/$f.jls" ] ||
                    { echo "$outdir/$f.jls is left by a task killed between its moves; move it aside" >&2; exit 1; }
            done
            widths="$widths${widths:+,}$i"
        fi
    done
    [ -n "$widths" ] || [ ! -f "$outdir/budget.md" ] || { echo "nothing is missing from $outdir"; exit 0; }
else
    [ -e "$outdir/jobs.tsv" ] && { echo "$outdir already holds a run; use missing" >&2; exit 1; }
fi

pre=$(sub $exp --cpus-per-task=1 --mem=16G --time=01:00:00 --job-name=${name}_pre \
          --wrap "$run -e 'using InferCell'")
if [ "$mode" != missing ]; then
    record commit "$head"
    record smoke "$smoke"
fi
record precompile "$pre"
dep="--dependency=afterok:$pre --kill-on-invalid-dep=yes"

if [ "$mode" = secheck ]; then
    record secheck "$(sub $exp $dep --cpus-per-task=1 --mem=16G --time=03:00:00 --job-name=${name}_secheck \
                          --wrap "$run $script secheck $outdir")"
else
    if [ -n "$widths" ]; then
        w=$(arr $exp $dep --array="$widths" --cpus-per-task="$threads" --mem="$mem" --time="$wall" \
                --job-name=${name}_width \
                --wrap "$run --threads=$threads --gcthreads=$((threads / 2)),1 $script width $outdir \$SLURM_ARRAY_TASK_ID")
        record width "$w"; dep="--dependency=afterany:$w"
    fi
    record report "$(sub $exp $dep --cpus-per-task=1 --mem=16G --time=01:00:00 --job-name=${name}_report \
                         --wrap "$run $script report $outdir")"
fi
echo "submitted at $(git rev-parse --short HEAD) ($mode):"
cat "$outdir/jobs.tsv"
