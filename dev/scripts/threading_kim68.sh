#!/bin/bash
# KIM-68: the threading diagnosis (dev/scripts/threading_kim68.jl), from a worktree
# at a pushed commit.
#   dev/scripts/threading_kim68.sh <outdir>
# precompiles on milan, then submits every stage, and the report after them all.
# With SMOKE=1 in the environment it is a plumbing pass, at tiny sizes and few
# threads. Job ids go to <outdir>/jobs.tsv, logs to <outdir>/logs.
set -euo pipefail
outdir=$(realpath -m "$1")
source "$(dirname "$0")/launcher_guard.sh"
smoke=${SMOKE:-0}
script=dev/scripts/threading_kim68.jl
if [ "$smoke" = 1 ]; then big=4; mid=2; else big=64; mid=16; fi
ex() { echo "--export=ALL,SMOKE=$smoke${1:+,$1}"; }

pre=$(sub "$(ex)" --cpus-per-task=1 --mem=16G --time=01:00:00 --job-name=kim68_pre \
          --wrap "$run -e 'using InferCell'")
record precompile "$pre"
dep="--dependency=afterok:$pre"
all=""
add() { record "$1" "$2"; all="$all:$2"; }

add serial "$(sub "$(ex)" "$dep" --cpus-per-task=1 --mem=16G --time=01:00:00 --job-name=kim68_serial \
                  --wrap "$run $script serial $outdir")"
for n in $mid $big; do
    add "threaded$n" "$(sub "$(ex LOOPS=chunked+greedy,PREDICT=1)" "$dep" --cpus-per-task=$n \
                            --mem=$((n * 3))G --time=02:00:00 --job-name=kim68_threaded$n \
                            --wrap "$run --threads=$n $script threaded $outdir default")"
done
add "threaded${big}_gc_half_concurrent" \
    "$(sub "$(ex)" "$dep" --cpus-per-task=$big --mem=$((big * 3))G --time=01:00:00 \
           --job-name=kim68_gc_half_conc \
           --wrap "$run --threads=$big --gcthreads=$((big / 2)),1 $script threaded $outdir gc_half_concurrent")"
add "threaded${big}_gc_one" \
    "$(sub "$(ex)" "$dep" --cpus-per-task=$big --mem=$((big * 3))G --time=01:00:00 \
           --job-name=kim68_gc_one \
           --wrap "$run --threads=$big --gcthreads=1 $script threaded $outdir gc_one")"
add "threaded${big}_blas1" \
    "$(sub "$(ex OPENBLAS_NUM_THREADS=1)" "$dep" --cpus-per-task=$big --mem=$((big * 3))G \
           --time=01:00:00 --job-name=kim68_blas1 \
           --wrap "$run --threads=$big $script threaded $outdir blas1")"
# One node, 4 processes, each with its own heap and a quarter of the threads.
add "procs4x$mid" \
    "$(sub "$(ex)" "$dep" --nodes=1 --ntasks=4 --cpus-per-task=$mid --mem=$((4 * mid * 3))G \
           --time=01:00:00 --job-name=kim68_procs \
           --wrap "module load julia/1.10.5; srun --ntasks=4 --cpus-per-task=$mid julia --project --threads=$mid $script procs $outdir")"
add allocs "$(sub "$(ex)" "$dep" --cpus-per-task=1 --mem=16G --time=01:00:00 --job-name=kim68_allocs \
                  --wrap "$run $script allocs $outdir")"
record report "$(sub "$(ex)" --dependency=afterany$all --cpus-per-task=1 --mem=8G --time=00:20:00 \
                     --job-name=kim68_report --wrap "$run $script report $outdir")"
echo "submitted at $(git rev-parse --short HEAD), SMOKE=$smoke:"
cat "$outdir/jobs.tsv"
