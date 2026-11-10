#!/bin/bash
# KIM-92: submit the reduced-model pilot's stages (dev/scripts/corea_pilot_kim92.jl),
# from a worktree at a pushed commit. The design is frozen at prereg/kim-92.
#   dev/scripts/corea_pilot_kim92.sh <outdir> profile    precompile, dataset, profile (item 2)
#   dev/scripts/corea_pilot_kim92.sh <outdir> reprofile  precompile, profile, on the dataset
#                                                        already in <outdir> (new code, same cells)
#   dev/scripts/corea_pilot_kim92.sh <outdir> posterior  prep, 4 chains, the gate on their
#                                                        first NDRAW draws (needs the dataset)
#   dev/scripts/corea_pilot_kim92.sh <outdir> extend     the one extension of all 4 chains,
#                                                        then the gate on 2 NDRAW draws
#   dev/scripts/corea_pilot_kim92.sh <outdir> ppc        the PPC on the kept draws, then the merge
#   dev/scripts/corea_pilot_kim92.sh <outdir> all        SMOKE only: every stage in one chain
# With SMOKE=1 in the environment it is a plumbing pass at tiny sizes. Job ids go to
# <outdir>/jobs.tsv (appended across calls), logs to <outdir>/logs. CHAIN_TIME sets the
# chains' and the extension's wall limit (default 7 days, milan's maximum).
set -euo pipefail
outdir=$(realpath -m "$1")
what=${2:?say what to submit: profile, reprofile, posterior, extend, ppc or all}
source "$(dirname "$0")/launcher_guard.sh"
smoke=${SMOKE:-0}
script=dev/scripts/corea_pilot_kim92.jl
ex="--export=ALL,SMOKE=$smoke,NADAPT=${NADAPT:-},NDRAW=${NDRAW:-}"
if [ "$smoke" = 1 ]; then
    t_ds=01:00:00; t_prof=02:00:00; t_prep=02:00:00; t_chain=04:00:00; t_ppc=01:00:00; ppc_threads=2
    ndraw=${NDRAW:-5}
else
    t_ds=12:00:00; t_prof=2-00:00:00; t_prep=2-00:00:00; t_chain=${CHAIN_TIME:-7-00:00:00}
    t_ppc=2-00:00:00; ppc_threads=16; ndraw=${NDRAW:-200}
fi
[ "$what" = all ] && [ "$smoke" != 1 ] && { echo "'all' is for SMOKE=1 only" >&2; exit 1; }
need() { [ -f "$outdir/$1.jls" ] || [ "$what" = all ] || { echo "$outdir/$1.jls is missing" >&2; exit 1; }; }
after() { if [ -n "$1" ]; then echo "--dependency=afterok:$1"; fi; }

profile() {   # -> $ds_job
    local pre pr
    pre=$(sub "$ex" --cpus-per-task=1 --mem=16G --time=01:00:00 --job-name=kim92_pre \
              --wrap "$run -e 'using InferCell'")
    record precompile "$pre"
    ds_job=$(sub "$ex" --dependency=afterok:$pre --cpus-per-task=1 --mem=32G --time=$t_ds \
                 --job-name=kim92_dataset --wrap "$run $script dataset $outdir")
    record dataset "$ds_job"
    pr=$(sub "$ex" --dependency=afterok:$ds_job --cpus-per-task=1 --mem=32G --time=$t_prof \
             --job-name=kim92_profile --wrap "$run $script profile $outdir")
    record profile "$pr"
}
reprofile() {
    local pre
    pre=$(sub "$ex" --cpus-per-task=1 --mem=16G --time=01:00:00 --job-name=kim92_pre \
              --wrap "$run -e 'using InferCell'")
    record precompile "$pre"
    record profile "$(sub "$ex" --dependency=afterok:$pre --cpus-per-task=1 --mem=32G --time=$t_prof \
                          --job-name=kim92_profile --wrap "$run $script profile $outdir")"
}
posterior() {   # $1: a job to wait for, or empty -> $chains_job
    local p c ids=()
    p=$(sub "$ex" $(after "$1") --cpus-per-task=1 --mem=32G --time=$t_prep --job-name=kim92_prep \
            --wrap "$run $script prep $outdir")
    record prep "$p"
    for c in 1 2 3 4; do
        ids+=("$(sub "$ex" --dependency=afterok:$p --cpus-per-task=1 --mem=32G --time=$t_chain \
                     --job-name=kim92_chain$c --wrap "$run $script chain $outdir $c")")
        record chain "${ids[-1]}"
    done
    chains_job=$(IFS=:; echo "${ids[*]}")
    record check "$(sub "$ex" --dependency=afterok:$chains_job --cpus-per-task=1 --mem=16G \
                        --time=01:00:00 --job-name=kim92_check --wrap "$run $script check $outdir $ndraw")"
}
extend() {   # $1: a job to wait for, or empty -> $ext_job
    local c ids=()
    for c in 1 2 3 4; do
        ids+=("$(sub "$ex" $(after "$1") --cpus-per-task=1 --mem=32G --time=$t_chain \
                     --job-name=kim92_extend$c --wrap "$run $script extend $outdir $c")")
        record extend "${ids[-1]}"
    done
    ext_job=$(IFS=:; echo "${ids[*]}")
    record check "$(sub "$ex" --dependency=afterok:$ext_job --cpus-per-task=1 --mem=16G \
                        --time=01:00:00 --job-name=kim92_check2 \
                        --wrap "$run $script check $outdir $((2 * ndraw))")"
}
ppc() {   # $1: a job to wait for, or empty
    local p
    p=$(sub "$ex" $(after "$1") --cpus-per-task=$ppc_threads --mem=64G --time=$t_ppc \
            --job-name=kim92_ppc \
            --wrap "$run --threads=$ppc_threads --gcthreads=$((ppc_threads / 2)),1 $script ppc $outdir")
    record ppc "$p"
    record merge "$(sub "$ex" --dependency=afterok:$p --cpus-per-task=1 --mem=16G --time=01:00:00 \
                        --job-name=kim92_merge --wrap "$run $script merge $outdir")"
}

case $what in
profile) profile ;;
reprofile) need dataset; reprofile ;;
posterior) need dataset; posterior "" ;;
extend) need chain_4; extend "" ;;
ppc) need chain_4; ppc "" ;;
all) profile; posterior "$ds_job"; extend "$chains_job"; ppc "$ext_job" ;;
*) echo "unknown: $what" >&2; exit 1 ;;
esac
echo "submitted $what at $(git rev-parse --short HEAD), SMOKE=$smoke:"
cat "$outdir/jobs.tsv"
