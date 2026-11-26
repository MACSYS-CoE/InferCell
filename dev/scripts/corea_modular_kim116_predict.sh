#!/bin/bash
# KIM-116 item 8: submit the plug-in's D3 scores (dev/scripts/corea_modular_kim116_predict.jl),
# from a clean worktree at a pushed commit.
#   dev/scripts/corea_modular_kim116_predict.sh <outdir> runs   precompile, prep, the 20
#                       posterior chunks as an array (KIM-95's sizes), then collect and
#                       KIM-95's scorer
# With SMOKE=1 in the environment it is a plumbing pass at KIM-95's smoke sizes, on the
# KIM-116 smoke pipeline's cut draws. MODULAR (optional) names the modular run's outdir.
# Job ids go to <outdir>/jobs.tsv, logs to <outdir>/logs. THREADS sets each chunk's threads
# (default 32), RUN_TIME their wall limit.
set -euo pipefail
outdir=$(realpath -m "$1")
what=${2:?say what to submit: runs}
source "$(dirname "$0")/launcher_guard.sh"
smoke=${SMOKE:-0}
script=dev/scripts/corea_modular_kim116_predict.jl
scorer=dev/scripts/corea_predict_kim95_score.jl
ex="--export=ALL,SMOKE=$smoke${MODULAR:+,MODULAR=$MODULAR}"
if [ "$smoke" = 1 ]; then
    threads=2; mem=16G; t_prep=01:00:00; t_run=02:00:00; chunks=2
else
    # KIM-95's: a reduced-model cell at 6,300 s allocates ~25 GiB; ~6 h a chunk on 32 threads.
    threads=${THREADS:-32}; mem=128G; t_prep=02:00:00; t_run=${RUN_TIME:-2-00:00:00}; chunks=20
fi
# The chunk count must be the driver's (KIM-95's CHUNKS, which it includes).
[ "$smoke" = 1 ] || grep -q '^const CHUNKS = Dict(:posterior => SMOKE ? 2 : 20,' dev/scripts/corea_predict_kim95.jl ||
    { echo "KIM-95's posterior chunk count is not 20" >&2; exit 1; }
jl="$run --threads=$threads --gcthreads=$((threads > 1 ? threads / 2 : 1)),1"

case $what in
runs)
    pre=$(sub "$ex" --cpus-per-task=1 --mem=16G --time=01:00:00 --job-name=kim116p_pre \
              --wrap "$run -e 'using InferCell'")
    record precompile "$pre"
    prep=$(sub "$ex" --dependency=afterok:$pre --cpus-per-task=1 --mem=32G --time=$t_prep \
               --job-name=kim116p_prep --wrap "$run $script prep $outdir")
    record prep "$prep"
    arr=$(sub "$ex" --dependency=afterok:$prep --array=1-$chunks --cpus-per-task=$threads --mem=$mem \
              --time=$t_run --job-name=kim116p_posterior --output="$outdir/logs/%x_%A_%a.log" \
              --wrap "$jl $script run $outdir \$SLURM_ARRAY_TASK_ID")
    record run_posterior "$arr"
    col=$(sub "$ex" --dependency=afterok:$arr --cpus-per-task=1 --mem=32G --time=01:00:00 \
              --job-name=kim116p_collect --wrap "$run $script collect $outdir && $run $scorer $outdir")
    record collect "$col" ;;
*) echo "unknown: $what" >&2; exit 1 ;;
esac
echo "submitted $what at $(git rev-parse --short HEAD), SMOKE=$smoke:"
cat "$outdir/jobs.tsv"
