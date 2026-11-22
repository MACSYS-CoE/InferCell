#!/bin/bash
# KIM-95: submit the scored predictions' stages (dev/scripts/corea_predict_kim95.jl), from
# a worktree at a pushed commit.
#   dev/scripts/corea_predict_kim95.sh <outdir> sizing   precompile, prep, sizing (reads no truth)
#   dev/scripts/corea_predict_kim95.sh <outdir> runs     the truth, posterior and prior chunks as
#                                                        arrays, then collect and the scores (needs
#                                                        prep, and the frozen sizes in the driver)
#   dev/scripts/corea_predict_kim95.sh <outdir> resolve  the resolving pass: posterior only,
#                                                        SIZING_DRAWS (40) × SIZING_CELLS (100), to
#                                                        sizing_resolve.* (needs prep)
#   dev/scripts/corea_predict_kim95.sh <outdir> all      SMOKE only: every stage in one chain
# With SMOKE=1 in the environment it is a plumbing pass at tiny sizes. Job ids go to
# <outdir>/jobs.tsv (appended across calls), logs to <outdir>/logs. THREADS sets each
# prediction job's threads (default 32), RUN_TIME their wall limit.
set -euo pipefail
outdir=$(realpath -m "$1")
what=${2:?say what to submit: sizing, runs or all}
source "$(dirname "$0")/launcher_guard.sh"
smoke=${SMOKE:-0}
script=dev/scripts/corea_predict_kim95.jl
scorer=dev/scripts/corea_predict_kim95_score.jl
ex="--export=ALL,SMOKE=$smoke,SIZING_DRAWS=${SIZING_DRAWS:-},SIZING_CELLS=${SIZING_CELLS:-}"
if [ "$smoke" = 1 ]; then
    threads=2; mem=16G; t_prep=01:00:00; t_size=02:00:00; t_run=02:00:00
else
    # KIM-68: a reduced-model cell at 6,300 s allocates ~25 GiB; KIM-42's MaxRSS was
    # 43.4 GiB at 64 threads.
    threads=${THREADS:-32}; mem=128G; t_prep=02:00:00; t_size=12:00:00; t_run=${RUN_TIME:-2-00:00:00}
fi
[ "$what" = all ] && [ "$smoke" != 1 ] && { echo "'all' is for SMOKE=1 only" >&2; exit 1; }
need() { [ -f "$outdir/$1.jls" ] || [ "$what" = all ] || { echo "$outdir/$1.jls is missing" >&2; exit 1; }; }
after() { if [ -n "$1" ]; then echo "--dependency=afterok:$1"; fi; }
jl="$run --threads=$threads --gcthreads=$((threads > 1 ? threads / 2 : 1)),1"
# The chunk counts, read from the driver so that they cannot disagree with it.
chunks() {
    if [ "$smoke" = 1 ]; then echo 2; return; fi
    sed -n "s/^const CHUNKS = Dict(:posterior => SMOKE ? 2 : \([0-9]*\), :prior => SMOKE ? 2 : \([0-9]*\), :truth => SMOKE ? 2 : \([0-9]*\)).*/\1 \2 \3/p" "$script" |
        awk -v a="$1" '{ print (a == "posterior" ? $1 : a == "prior" ? $2 : $3) }'
}

resolve() {
    local j ex_r="--export=ALL,SMOKE=$smoke,SIZING_TAG=resolve,SIZING_ARMS=posterior"
    ex_r="$ex_r,SIZING_DRAWS=${SIZING_DRAWS:-40},SIZING_CELLS=${SIZING_CELLS:-100}"
    j=$(sub "$ex_r" --cpus-per-task=$threads --mem=$mem --time=${RESOLVE_TIME:-16:00:00} \
            --job-name=kim95_resolve --wrap "$jl $script sizing $outdir")
    record resolve "$j"
}
sizing() {   # -> $size_job
    local pre prep_job
    pre=$(sub "$ex" --cpus-per-task=1 --mem=16G --time=01:00:00 --job-name=kim95_pre \
              --wrap "$run -e 'using InferCell'")
    record precompile "$pre"
    prep_job=$(sub "$ex" --dependency=afterok:$pre --cpus-per-task=1 --mem=32G --time=$t_prep \
                   --job-name=kim95_prep --wrap "$run $script prep $outdir")
    record prep "$prep_job"
    size_job=$(sub "$ex" --dependency=afterok:$prep_job --cpus-per-task=$threads --mem=$mem \
                   --time=$t_size --job-name=kim95_sizing --wrap "$jl $script sizing $outdir")
    record sizing "$size_job"
}
runs() {   # $1: a job to wait for, or empty
    local arm k ids=() j
    if [ "$smoke" != 1 ] && ! grep -q '^const FROZEN = true' "$script"; then
        echo "the sizes are not frozen in $script" >&2; exit 1
    fi
    for arm in truth posterior prior; do
        k=$(chunks $arm)
        [ -n "$k" ] && [ "$k" -ge 1 ] || { echo "no chunk count for $arm" >&2; exit 1; }
        j=$(sub "$ex" $(after "$1") --array=1-$k --cpus-per-task=$threads --mem=$mem --time=$t_run \
                --job-name=kim95_$arm --output="$outdir/logs/%x_%A_%a.log" \
                --wrap "$jl $script run $outdir $arm \$SLURM_ARRAY_TASK_ID")
        record "run_$arm" "$j"
        ids+=("$j")
    done
    j=$(sub "$ex" --dependency=afterok:$(IFS=:; echo "${ids[*]}") --cpus-per-task=1 --mem=32G \
            --time=01:00:00 --job-name=kim95_collect \
            --wrap "$run $script collect $outdir && $run $scorer $outdir")
    record collect "$j"
}

case $what in
sizing) sizing ;;
runs) need prep; runs "" ;;
resolve) need prep; resolve ;;
all) sizing; runs "$size_job" ;;
*) echo "unknown: $what" >&2; exit 1 ;;
esac
echo "submitted $what at $(git rev-parse --short HEAD), SMOKE=$smoke:"
cat "$outdir/jobs.tsv"
