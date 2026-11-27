#!/bin/bash
# KIM-195: submit the realistic mean-field timing (route_cost_meanfield.jl), from a
# detached worktree at a pushed commit.
#   dev/scripts/fullcell_route_cost/route_cost_meanfield.sh <outdir> [what]
# what = all (default): precompile both environments; the SpeedyCell reference; then
#        check (after the reference), time and library (after check, which writes the
#        ranking), threaded (after time), and clamp (SpeedyCell's dNTPs, then the mean
#        field's under both payment forms, in one job).
# what = reference | check | time | threaded | library | clamp: that stage alone, after its own
#        precompile (check needs <outdir>/reference; the others reuse ranking.jls).
# SMOKE=1 in the environment: a plumbing pass (300 s window, 2 repeats).
# Julia 1.11.9 (LinearSolve 3.87.0 does not load on 1.10.5). The reference runs in
# KIM-183's environment and depot, where SpeedyCell is developed by path.
# Job ids go to <outdir>/jobs.tsv, logs to <outdir>/logs.
set -euo pipefail
outdir=$(realpath -m "$1")
what=${2:-all}
case $what in all|reference|check|time|threaded|library|clamp) ;; *) echo "unknown: $what" >&2; exit 1 ;; esac
cd "$(dirname "$(realpath "$0")")/../../.."
source dev/scripts/launcher_guard.sh
smoke=${SMOKE:-0}
env=dev/scripts/fullcell_route_cost
script=$env/route_cost_meanfield.jl
nthr=4
if [ "$smoke" = 1 ]; then t=02:00:00; else t=48:00:00; fi
jl='export PATH=/fred/oz022/tkimpson/julia-1.11.9/bin:$PATH; export JULIA_PKG_OFFLINE=true; julia --project='"$env"
jlref='export PATH=/fred/oz022/tkimpson/julia-1.11.9/bin:$PATH; export JULIA_DEPOT_PATH=/fred/oz022/tkimpson/infercell-kim183-2026-11-26/depot:/fred/oz022/tkimpson/julia_depot; export JULIA_PKG_OFFLINE=true; julia --project=dev/scripts/fullcell_timing_cycle'
ex() { echo "--export=ALL,SMOKE=$smoke,OPENBLAS_NUM_THREADS=1${1:+,$1}"; }

pre=$(sub "$(ex)" --cpus-per-task=1 --mem=16G --time=01:00:00 --job-name=kim195_pre \
          --wrap "$jl -e 'using OrdinaryDiffEq, LinearSolve, SparseArrays, ForwardDiff, SciMLBase; include(\"$env/model.jl\"); load_model()' && $jlref -e 'using SpeedyCell'")
record precompile "$pre"
job() {   # $1 stage, $2 dependency list (afterok), $3 cpus
    local id cpus=${3:-1} wrap
    if [ "$1" = reference ]; then
        wrap="$jlref $env/reference_speedycell.jl $outdir"
    elif [ "$1" = clamp ]; then
        wrap="( $jlref $env/clamp_speedycell.jl $outdir ) && $jl $script clamp $outdir"
    elif [ "$cpus" -gt 1 ]; then
        wrap="$jl --threads=$cpus --gcthreads=$cpus,1 $script $1 $outdir"
    else
        wrap="$jl $script $1 $outdir"
    fi
    id=$(sub "$(ex)" --dependency=afterok:$pre${2:+:$2} --cpus-per-task=$cpus --mem=16G --time=$t \
             --job-name=kim195_$1 --wrap "$wrap")
    record "$1" "$id"
    echo "$id"
}
case $what in
all)
    rid=$(job reference)
    cid=$(job check "$rid")
    tid=$(job time "$cid")
    job library "$cid" >/dev/null
    job threaded "$cid:$tid" $nthr >/dev/null
    job clamp >/dev/null ;;
threaded) job threaded "" $nthr >/dev/null ;;
*) job "$what" >/dev/null ;;
esac
echo "submitted $what at $(git rev-parse --short HEAD), SMOKE=$smoke:"
cat "$outdir/jobs.tsv"
