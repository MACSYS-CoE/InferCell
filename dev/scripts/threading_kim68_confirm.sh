#!/bin/bash
# KIM-68, after the gate: confirms the diagnosis's one big gain, a concurrent GC
# sweep thread (`--gcthreads=N,1`), before KIM-42 runs on it. From a worktree at a
# pushed commit:
#   dev/scripts/threading_kim68_confirm.sh <outdir>
# runs dev/scripts/threading_kim68.jl's `threaded` stage (the new `_each`, and one
# full `predict` call) three times:
#   - 64 threads, `--gcthreads=32,1`, on a node other than the diagnosis's (dave13);
#   - 64 threads, `--gcthreads=8,1`: does the number of mark threads matter?
#   - 16 threads, `--gcthreads=8,1`: the default's 8 mark threads, plus the sweep.
# Then the report, which reads only <outdir>.
set -euo pipefail
outdir=$(realpath -m "$1")
source "$(dirname "$0")/launcher_guard.sh"
script=dev/scripts/threading_kim68.jl
ex="--export=ALL,SMOKE=0,LOOPS=greedy,PREDICT=1"

pre=$(sub "$ex" --cpus-per-task=1 --mem=16G --time=01:00:00 --job-name=kim68c_pre \
          --wrap "$run -e 'using InferCell'")
record precompile "$pre"
all=""
confirm() {     # label threads gcthreads [sbatch options]
    local id
    id=$(sub "$ex" --dependency=afterok:$pre --cpus-per-task=$2 --mem=$(($2 * 3))G \
             --time=01:00:00 --job-name=kim68c_$1_$2 "${@:4}" \
             --wrap "$run --threads=$2 --gcthreads=$3 $script threaded $outdir $1")
    record "threaded$2_$1" "$id"
    all="$all:$id"
}
confirm gc32c1 64 32,1 --exclude='dave[13,148-149]'
confirm gc8c1 64 8,1
confirm gc8c1 16 8,1
record report "$(sub "$ex" --dependency=afterany$all --cpus-per-task=1 --mem=8G --time=00:20:00 \
                     --job-name=kim68c_report --wrap "$run $script report $outdir")"
echo "submitted at $(git rev-parse --short HEAD):"
cat "$outdir/jobs.tsv"
