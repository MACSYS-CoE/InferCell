#!/bin/bash
# KIM-116: submit dev/scripts/corea_modular_kim116.jl's stages to Slurm, from a clean
# worktree at a pushed commit. The design is frozen at prereg/kim-116.
#   dev/scripts/corea_modular_kim116.sh <mode> <outdir> [<variant>] [<k>]
# Modes (<variant> is cut or plugin):
#   profile          item 2: one job, each stage's density and gradient (profile.md)
#   stage1 v         v's warm metric, then its 4 chains once it is done
#   check v k        the stage-1 gate on k draws per chain (v/gate_<k>.md)
#   extend v         the one extension of v's 4 chains, after the gate failed at NDRAW
#   stage2 v         v's K = 200 inner chains, as an array of 20 tasks of 10
#   finish           each variant's PPC, then the merge (modular_summary.md)
#   all              the whole pipeline in one dependency chain, without the extension:
#                    stage1, check, stage2 and finish for both variants (for SMOKE=1)
# Job ids go to <outdir>/jobs.tsv, logs to <outdir>/logs, and each submission's commit
# to <outdir>/COMMIT. SMOKE=1 in the environment passes through to the jobs.
set -euo pipefail
mode=${1:?give the mode}
outdir=$(realpath -m "${2:?give the output directory}")
v=${3:-}
source "$(dirname "$0")/launcher_guard.sh"
script=dev/scripts/corea_modular_kim116.jl
git rev-parse HEAD >> "$outdir/COMMIT"
smoke=${SMOKE:-0}
if [ "$smoke" = 1 ]; then long=02:00:00; inner=1; per=4; else long=7-00:00:00; inner=20; per=10; fi
env=(--export=ALL,SMOKE="$smoke")
need_v() { case "$1" in cut|plugin) ;; *) echo "give the variant: cut or plugin" >&2; exit 1;; esac; }
after() { [ -n "$1" ] && echo "--dependency=afterok:$1" || true; }

submit_stage1() {  # variant [dependency] -> prints the chains' ids, colon-separated
    local m ids=() j
    m=$(sub --cpus-per-task=1 --mem=8G --time=12:00:00 --job-name="kim116_metric_$1" "${env[@]}" \
            $(after "${2:-}") --wrap "$run $script metric $outdir $1")
    record "metric_$1" "$m"
    for c in 1 2 3 4; do
        j=$(sub --cpus-per-task=1 --mem=8G --time="$long" --job-name="kim116_chain_$1_$c" "${env[@]}" \
                --dependency=afterok:"$m" --wrap "$run $script chain $outdir $1 $c")
        record "chain_$1" "$j"
        ids+=("$j")
    done
    (IFS=:; echo "${ids[*]}")
}
submit_check() {  # variant k dependency
    local j
    j=$(sub --cpus-per-task=1 --mem=8G --time=01:00:00 --job-name="kim116_check_$1" "${env[@]}" \
            $(after "${3:-}") --wrap "$run $script check $outdir $1 $2")
    record "check_$1" "$j"
    echo "$j"
}
submit_stage2() {  # variant dependency
    local j
    j=$(sub --cpus-per-task=1 --mem=8G --time=1-00:00:00 --job-name="kim116_stage2_$1" "${env[@]}" \
            --array=1-"$inner" $(after "${2:-}") \
            --wrap "a=\$(( (SLURM_ARRAY_TASK_ID - 1) * $per + 1 )); $run $script stage2 $outdir $1 \$a \$(( a + $per - 1 ))")
    record "stage2_$1" "$j"
    echo "$j"
}
submit_finish() {  # dependency
    local ids=() j
    for w in cut plugin; do
        j=$(sub --cpus-per-task=8 --mem=16G --time=04:00:00 --job-name="kim116_ppc_$w" \
                "${env[@]}" $(after "${1:-}") --wrap "$run --threads=8 $script ppc $outdir $w")
        record "ppc_$w" "$j"
        ids+=("$j")
    done
    j=$(sub --cpus-per-task=1 --mem=8G --time=01:00:00 --job-name=kim116_merge "${env[@]}" \
            --dependency=afterok:"${ids[0]}":"${ids[1]}" --wrap "$run $script merge $outdir")
    record merge "$j"
}

case "$mode" in
profile)
    [ -e "$outdir/profile.md" ] && { echo "$outdir/profile.md already exists" >&2; exit 1; }
    j=$(sub --cpus-per-task=1 --mem=16G --time=04:00:00 --job-name=kim116_profile "${env[@]}" \
            --wrap "$run $script profile $outdir")
    record profile "$j" ;;
stage1) need_v "$v"; submit_stage1 "$v" >/dev/null ;;
check) need_v "$v"; submit_check "$v" "${4:?give k}" "" >/dev/null ;;
extend)
    need_v "$v"
    for c in 1 2 3 4; do
        j=$(sub --cpus-per-task=1 --mem=8G --time="$long" --job-name="kim116_extend_${v}_$c" \
                "${env[@]}" --wrap "$run $script extend $outdir $v $c")
        record "extend_$v" "$j"
    done ;;
stage2) need_v "$v"; submit_stage2 "$v" "" >/dev/null ;;
finish) submit_finish "" ;;
all)
    k=$([ "$smoke" = 1 ] && echo 5 || echo 200)
    s2=()
    for w in cut plugin; do
        ch=$(submit_stage1 "$w" "")
        ck=$(submit_check "$w" "$k" "$ch")
        s2+=("$(submit_stage2 "$w" "$ck")")
    done
    submit_finish "${s2[0]}:${s2[1]}" ;;
*)
    echo "unknown mode $mode" >&2; exit 1 ;;
esac
echo "submitted at $(git rev-parse --short HEAD):"
cat "$outdir/jobs.tsv"
