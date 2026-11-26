#!/bin/bash
# KIM-178: submit dev/scripts/corea_replicates/corea_replicates.jl, from a worktree at a
# pushed commit. The design is frozen at prereg/kim-178.
#
#   corea_replicates.sh <outdir> kim92            prereg §5: the fitted MAP and Laplace on KIM-92's dataset
#   corea_replicates.sh <outdir> replicates       the ten (seeds 17_200_002–011), one array task each
#   corea_replicates.sh <outdir> spares <k>...    spare k of 1–5 (seed 17_200_011 + k), only after a
#                                                 recorded failure (prereg §3)
#   SMOKE=1 corea_replicates.sh <outdir> smoke    seed 17_200_199 at tiny sizes, then the kim92 stage on it
#
# Job ids are appended to <outdir>/jobs.tsv, and logs go to <outdir>/logs.
set -euo pipefail
outdir=$(realpath -m "$1")
what=${2:?say what to submit: kim92, replicates, spares or smoke}
source "$(dirname "$0")/../launcher_guard.sh"
smoke=${SMOKE:-0}
script=dev/scripts/corea_replicates/corea_replicates.jl
ex="--export=ALL,SMOKE=$smoke"
if [ "$smoke" = 1 ]; then
    t_rep=03:00:00; t_k92=03:00:00
else
    # KIM-92's prep took 7.3 h (job 20545098), and its metric 42 gradients in 2230 s (job 20651329).
    # The fitted MAP is a second optimisation of about the first's size.
    t_rep=2-00:00:00; t_k92=1-00:00:00
fi
case "$what" in smoke) [ "$smoke" = 1 ] || { echo "'smoke' needs SMOKE=1" >&2; exit 1; } ;;
                *) [ "$smoke" = 1 ] && { echo "SMOKE=1 is for 'smoke' only" >&2; exit 1; } ;; esac

pre=$(sub "$ex" --cpus-per-task=1 --mem=16G --time=01:00:00 --job-name=kim178_pre \
          --wrap "$run -e 'using InferCell'")
record precompile "$pre"
git rev-parse HEAD >> "$outdir/COMMIT"

case "$what" in
kim92)
    record kim92 "$(sub "$ex" --dependency=afterok:$pre --cpus-per-task=1 --mem=32G --time=$t_k92 \
                        --job-name=kim178_kim92 --wrap "$run $script kim92 $outdir")" ;;
replicates)
    # Array task i is seed 17_200_001 + i.
    record replicates "$(sub "$ex" --dependency=afterok:$pre --array=1-10 --cpus-per-task=1 --mem=32G \
                             --time=$t_rep --job-name=kim178_rep \
                             --output="$outdir/logs/%x_%A_%a.log" \
                             --wrap "$run $script replicate $outdir \$((17200001 + SLURM_ARRAY_TASK_ID))")" ;;
spares)
    shift 2
    [ $# -ge 1 ] || { echo "name the spares, 1 to 5" >&2; exit 1; }
    for k in "$@"; do
        [[ "$k" =~ ^[1-5]$ ]] || { echo "a spare is 1 to 5, not $k" >&2; exit 1; }
        record "spare_$k" "$(sub "$ex" --dependency=afterok:$pre --cpus-per-task=1 --mem=32G \
                                 --time=$t_rep --job-name=kim178_spare$k \
                                 --wrap "$run $script replicate $outdir $((17200011 + k))")"
    done ;;
smoke)
    r=$(sub "$ex" --dependency=afterok:$pre --cpus-per-task=1 --mem=16G --time=$t_rep \
            --job-name=kim178_smoke_rep --wrap "$run $script replicate $outdir 17200199")
    record smoke_replicate "$r"
    record smoke_kim92 "$(sub "$ex" --dependency=afterok:$r --cpus-per-task=1 --mem=16G --time=$t_k92 \
                              --job-name=kim178_smoke_kim92 \
                              --wrap "$run $script kim92 $outdir $outdir/seed_17200199")" ;;
*) echo "unknown: $what" >&2; exit 1 ;;
esac
cat "$outdir/jobs.tsv"
