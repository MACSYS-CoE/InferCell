#!/bin/bash
# Probabilistic trajectories (dev/scripts/trajectories_predict.jl). It runs the tree it
# lives in: keep a detached worktree at a pushed commit for the run, call this launcher
# from it, and leave it untouched until the run ends.
#   dev/scripts/trajectories_predict.sh <outdir> run
#       precompile, the release's pre-check, the two replays, the dataset, then every check;
#   dev/scripts/trajectories_predict.sh <outdir> rest
#       the release, then the report. It refuses unless replay.md ends in ALL PASS (a smoke
#       pass goes on regardless), and refuses a second time. `rerun c` needs ALL PASS too;
#   dev/scripts/trajectories_predict.sh <outdir> bulk <a worktree at dfe570a>
#       compares the dataset's bulk data with V7's commit, in bulk_check.md;
#   dev/scripts/trajectories_predict.sh <outdir> rerun <stage> [<array spec, for a>]
#       one stage again, once its output has been moved aside. It waits only on precompile,
#       so rerun `replay` or `report` after the stages they read have finished;
#   dev/scripts/trajectories_predict.sh <outdir> progress
#       lists what is written, and the run's jobs in the queue.
# SMOKE=1 in the launcher's environment runs the smoke pass, and is passed to the jobs.
# AWALL, BWALL and CWALL set the stages' wall limits, so that raising one needs no commit.
# Job ids go to <outdir>/jobs.tsv, logs to <outdir>/logs.
set -euo pipefail
outdir=$(realpath -m "$1")
mode=${2:-run}
smoke=${SMOKE:-0}
if [ "$smoke" = 1 ]; then
    name=trajsmoke; big=8; bigmem=32G; bwall=01:00:00; cwall=01:00:00; pwall=01:00:00; pthreads=8; pmem=32G; awall=01:00:00
else
    name=traj; big=64; bigmem=192G; bwall=${BWALL:-08:00:00}; cwall=${CWALL:-24:00:00}; pwall=02:00:00; pthreads=16; pmem=48G; awall=${AWALL:-14:00:00}
fi
all="${name}_pre,${name}_a,${name}_b,${name}_cpre,${name}_c,${name}_data,${name}_replay,${name}_report,${name}_bulk"

if [ "$mode" = progress ]; then
    ls "$outdir"
    squeue -h -u "$USER" -n "$all" -o "%.14i %.18j %.8T %.10M %R" || true
    exit 0
fi
case "$mode" in run|rest|bulk|rerun) ;; *) echo "unknown mode $mode" >&2; exit 1;; esac

# The tree this launcher lives in, wherever it is called from.
cd "$(dirname "$(realpath "$0")")/../.."
source dev/scripts/launcher_guard.sh
script=dev/scripts/trajectories_predict.jl
exp="--export=ALL,SMOKE=$smoke"
head=$(git rev-parse HEAD)

if [ "$mode" = run ]; then
    [ ! -e "$outdir/jobs.tsv" ] || { echo "$outdir already holds a run" >&2; exit 1; }
    record commit "$head"
    record smoke "$smoke"
else
    ran() { awk -F'\t' -v k="$1" '$1 == k { print $2; exit }' "$outdir/jobs.tsv"; }
    [ "$(ran commit)" = "$head" ] ||
        { echo "the run is at '$(ran commit)', this tree at $head" >&2; exit 1; }
    [ "$(ran smoke)" = "$smoke" ] ||
        { echo "the run has SMOKE='$(ran smoke)', this launcher SMOKE=$smoke" >&2; exit 1; }
fi
[ "$mode" != rest ] || [ -z "$(ran c)" ] ||
    { echo "the release is already job $(ran c); to run it again, use rerun c" >&2; exit 1; }
if [ "$mode" = rest ] || { [ "$mode" = rerun ] && [ "${3:-}" = c ]; }; then
    if [ "$smoke" != 1 ] && [ "$(tail -n 1 "$outdir/replay.md" 2>/dev/null)" != "ALL PASS" ]; then
        echo "$outdir/replay.md does not end in ALL PASS, so the release does not go ahead" >&2
        exit 1
    fi
fi

pre=$(sub $exp --cpus-per-task=1 --mem=16G --time=01:00:00 --job-name=${name}_pre \
          --wrap "$run -e 'using InferCell'")
record precompile "$pre"
dep="--dependency=afterok:$pre --kill-on-invalid-dep=yes"

job_a() {
    sbatch --parsable "${nodes[@]}" --ntasks=1 --output="$outdir/logs/%x_%A_%a.log" $exp $dep \
        --array="$1" --cpus-per-task=1 --mem=6G --time=$awall --job-name=${name}_a \
        --wrap "$run $script a $outdir \$SLURM_ARRAY_TASK_ID"
}
job_big() {  # stage, job name, walltime
    sub $exp $dep --cpus-per-task=$big --mem=$bigmem --time=$3 --job-name=${name}_$2 \
        --wrap "$run --threads=$big --gcthreads=$((big / 2)),1 $script $1 $outdir"
}
job_pre() {
    sub $exp $dep --cpus-per-task=$pthreads --mem=$pmem --time=$pwall --job-name=${name}_cpre \
        --wrap "$run --threads=$pthreads --gcthreads=$((pthreads / 2)),1 $script c_precheck $outdir"
}
job_small() {  # stage, job name, memory, extra sbatch flags
    sub $exp ${4:-$dep} --cpus-per-task=1 --mem=$3 --time=01:00:00 --job-name=${name}_$2 \
        --wrap "$run $script $1 $outdir"
}

# `id=$(...)` and not `record x "$(...)"`: under `set -e` only the first stops on a failed sbatch.
case "$mode" in
run)
    if [ "$smoke" = 1 ]; then blocks=0-2; else blocks=0-10; fi
    a=$(job_a "$blocks"); record a "$a"
    b=$(job_big b b "$bwall"); record b "$b"
    p=$(job_pre); record c_precheck "$p"
    d=$(job_small data data 8G); record data "$d"
    id=$(job_small replay replay 16G "--dependency=afterok:$a:$b:$p:$d --kill-on-invalid-dep=yes")
    record replay "$id"
    ;;
rest)
    c=$(job_big c c "$cwall"); record c "$c"
    id=$(job_small report report 16G "--dependency=afterok:$c --kill-on-invalid-dep=yes")
    record report "$id"
    ;;
bulk)
    wt=$(realpath "${3:?give a worktree at dfe570a}")
    id=$(sub $exp $dep --chdir="$wt" --cpus-per-task=1 --mem=8G --time=02:00:00 --job-name=${name}_bulk \
             --wrap "$run $PWD/dev/scripts/trajectories_bulk_check.jl $outdir")
    record bulk "$id"
    ;;
rerun)
    stage=${3:?name a stage}
    case "$stage" in
    a) id=$(job_a "${4:?give the array spec}");;
    b) id=$(job_big b b "$bwall");;
    c) id=$(job_big c c "$cwall");;
    c_precheck) id=$(job_pre);;
    data) id=$(job_small data data 8G);;
    replay) id=$(job_small replay replay 16G);;
    report) id=$(job_small report report 16G);;
    *) echo "unknown stage $stage" >&2; exit 1;;
    esac
    record "$stage" "$id"
    ;;
esac
echo "submitted at $(git rev-parse --short HEAD) ($mode):"
cat "$outdir/jobs.tsv"
