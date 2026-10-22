#!/bin/bash
# Issue #89, after B3: submit every stage of the rethink, chained.
#   dev/scripts/jensen_rethink_89.sh /fred/oz022/tkimpson/infercell-89-rethink-archive-<date>
# Stage 1: at B3's three M0 truths, the full arm (with a clip census) and the
# drains_own, drains_mf, enzymes_mf and popmean arms. Stage 2A: full ensembles at
# F10's posterior mean (800 cells) and eight posterior draws, and at ±1 posterior
# SD around m0_low and m0_high (800-cell full arms at those truths give δ). Every
# point gets a mean field; the arms that read it wait for it. The score waits for
# everything (afterany). SMOKE=1 runs 4 cells, a 120 s horizon and a few points.
set -euo pipefail
outdir=$1
case "$outdir" in */.worktrees/*|*/.claude/worktrees/*) echo "outdir is inside a worktree" >&2; exit 1;; esac
mkdir -p "$outdir/logs"
smoke=${SMOKE:-0}
job=dev/scripts/jensen_rethink_89.slurm
sub() { sbatch --parsable --export=ALL,SMOKE="$smoke" --output="$outdir/logs/%x_%A_%a.log" "$@"; }
if [ "$smoke" = 1 ]; then
    stage1=(m0_v6); refs=(f10_ref); targets=(f10_d1 low_krnadeg_+); per=2; t=00:30:00
else
    stage1=(m0_v6 m0_low m0_high)
    refs=(f10_ref)
    targets=($(julia --project dev/scripts/jensen_rethink_89.jl points | awk '{print $1}' | grep -E '^(f10_d|low_|high_)'))
    per=20; t=01:30:00
fi
declare -A mfjob
all=()
for p in "${stage1[@]}" "${refs[@]}" "${targets[@]}"; do
    mfjob[$p]=$(sub --time=02:00:00 $job meanfield "$p" "$outdir"); all+=("${mfjob[$p]}")
done
# ensemble <point> <arm> <ncells> <seed0> [dependency]
ens() {
    local n=$3 tasks
    if [ "$smoke" = 1 ]; then tasks=2; else tasks=$(( n / per )); fi
    all+=("$(sub --array=0-$((tasks - 1)) --time=$t ${5:+--dependency=afterok:$5} $job ensemble "$1" "$2" "$tasks" "$n" "$4" "$outdir")")
}
for p in "${stage1[@]}"; do
    if [ "$p" = m0_v6 ]; then ens "$p" full 200 40000; else ens "$p" full 800 40000; fi
    ens "$p" drains_own 200 40000
    ens "$p" drains_mf 200 40000 "${mfjob[$p]}"
    ens "$p" enzymes_mf 200 40000 "${mfjob[$p]}"
    all+=("$(sub --time=04:00:00 --mem=16G $job popmean "$p" "$outdir")")
done
for p in "${refs[@]}"; do ens "$p" full 800 40000; done
for p in "${targets[@]}"; do ens "$p" full 200 50000; done
deps=$(IFS=:; echo "${all[*]}")
s=$(sbatch --parsable --export=ALL,SMOKE="$smoke" --output="$outdir/logs/%x_score_%j.log" \
    --time=00:30:00 --mem=8G --dependency=afterany:$deps $job score "$outdir")
echo "stages ${all[*]}, score $s"
