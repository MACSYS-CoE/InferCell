#!/bin/bash
# KIM-40 and KIM-41: M0's evidence run (dev/scripts/m0_evidence_kim40.jl). It runs
# the tree it lives in: keep a detached worktree at a pushed commit for the run,
# call this launcher from it, and leave it untouched until the run ends, since
# each task loads InferCell from it when the task starts.
#   dev/scripts/m0_evidence_kim40.sh <outdir>
#       the whole run: precompile, prep, the work array, finish, score;
#   dev/scripts/m0_evidence_kim40.sh <outdir> missing
#       once none of the run's jobs is queued: resubmits the work tasks that wrote
#       no output, then finish for each (r, p) with none, then score;
#   dev/scripts/m0_evidence_kim40.sh <outdir> progress
#       counts the outputs, and lists the run's jobs in the queue.
# SMOKE=1 in the launcher's environment runs the smoke pass, and is passed to the
# jobs. CONCURRENCY (default 150) caps the work array's running tasks, and MEM
# (default 6G) is each work task's memory.
# finish and score wait on `afterany`, so a failed task never leaves a job pending
# forever: they refuse on a missing input, and `missing` fills it in.
# Job ids go to <outdir>/jobs.tsv, logs to <outdir>/logs.
set -euo pipefail
outdir=$(realpath -m "$1")
mode=${2:-run}
smoke=${SMOKE:-0}
if [ "$smoke" = 1 ]; then name=kim40smoke; nwork=7; nfin=2; else name=kim40; nwork=1150; nfin=100; fi
jobs="${name}_pre,${name}_prep,${name}_work,${name}_finish,${name}_score"

# Indices 0…n−1 of <dir> with an output (`NNNN_*.jls`, zero-padded to <width>),
# or, for `absent`, without one, comma-separated.
has() { compgen -G "$outdir/$1/$(printf "%0${3}d" "$2")_*.jls" > /dev/null; }
count() { local n=0 i; for ((i = 0; i < $2; i++)); do has "$1" $i "$3" && n=$((n + 1)); done; echo $n; }
absent() {
    local out=() i
    for ((i = 0; i < $2; i++)); do has "$1" $i "$3" || out+=($i); done
    (IFS=,; echo "${out[*]}")
}

if [ "$mode" = progress ]; then
    echo "work $(count work $nwork 4) of $nwork; finish $(count finish $nfin 3) of $nfin;" \
         "score $([ -f "$outdir/score.jls" ] && echo done || echo 'not yet')"
    squeue -h -u "$USER" -n "$jobs" -o "%.14i %.18j %.8T %.10M %R" || true
    exit 0
fi
[ "$mode" = run ] || [ "$mode" = missing ] || { echo "unknown mode $mode" >&2; exit 1; }

# The tree this launcher lives in, wherever it is called from.
cd "$(dirname "$(realpath "$0")")/../.."
source dev/scripts/launcher_guard.sh
script=dev/scripts/m0_evidence_kim40.jl
exp="--export=ALL,SMOKE=$smoke"
# `sub` (launcher_guard.sh), with each array task's log named by its index.
arr() { sbatch --parsable --partition=milan --exclude='dave[148-149]' --ntasks=1 \
               --output="$outdir/logs/%x_%A_%a.log" "$@"; }

if [ "$mode" = run ]; then
    [ -e "$outdir/jobs.tsv" ] && { echo "$outdir already holds a run; use missing" >&2; exit 1; }
    pre=$(sub $exp --cpus-per-task=1 --mem=16G --time=01:00:00 --job-name=${name}_pre \
              --wrap "$run -e 'using InferCell'")
    record precompile "$pre"
    prep=$(sub $exp --dependency=afterok:$pre --kill-on-invalid-dep=yes --cpus-per-task=1 \
               --mem=16G --time=01:00:00 --job-name=${name}_prep \
               --wrap "$run $script prep $outdir")
    record prep "$prep"
    work_list="0-$((nwork - 1))"; wdep="--dependency=afterok:$prep --kill-on-invalid-dep=yes"
    fin_list="0-$((nfin - 1))"
else
    [ -f "$outdir/manifest.jls" ] || { echo "no $outdir/manifest.jls: start a run first" >&2; exit 1; }
    [ -z "$(squeue -h -u "$USER" -n "$jobs")" ] ||
        { echo "the run still has jobs in the queue; wait for them" >&2; exit 1; }
    work_list=$(absent work $nwork 4); wdep=""
    fin_list=$(absent finish $nfin 3)
fi

fdep=""
if [ -n "$work_list" ]; then
    work=$(arr $exp $wdep --array="$work_list%${CONCURRENCY:-150}" --cpus-per-task=1 \
               --mem=${MEM:-6G} --time=24:00:00 --job-name=${name}_work \
               --wrap "$run $script work $outdir \$SLURM_ARRAY_TASK_ID")
    record work "$work"; fdep="--dependency=afterany:$work"
fi
sdep=$fdep
if [ -n "$fin_list" ]; then
    fin=$(arr $exp $fdep --array="$fin_list" --cpus-per-task=1 --mem=16G --time=03:00:00 \
              --job-name=${name}_finish \
              --wrap "$run $script finish $outdir \$SLURM_ARRAY_TASK_ID")
    record finish "$fin"; sdep="--dependency=afterany:$fin"
fi
if [ ! -f "$outdir/score.jls" ]; then
    s=$(sub $exp $sdep --cpus-per-task=1 --mem=16G --time=01:00:00 --job-name=${name}_score \
            --wrap "$run $script score $outdir")
    record score "$s"
fi
echo "submitted at $(git rev-parse --short HEAD) ($mode):"
cat "$outdir/jobs.tsv"
