#!/bin/bash
# Issue #89, after B3: 2A again with independent seeds (§12 2026-10-23), into the
# archive the first run wrote (its f10_ref/m0_low/m0_high 800-cell ensembles and
# the targets' mean fields are reused).
#   dev/scripts/jensen_rethink_89_transfer.sh /fred/oz022/tkimpson/infercell-89-rethink-archive-<date>
# Each target: 200 full cells at its own seeds. Each of twenty posterior draws:
# a mean field and 40 full cells. Then transfer, calibrate and score (afterany).
set -euo pipefail
outdir=$1
case "$outdir" in */.worktrees/*|*/.claude/worktrees/*) echo "outdir is inside a worktree" >&2; exit 1;; esac
mkdir -p "$outdir/logs"
job=dev/scripts/jensen_rethink_89.slurm
sub() { sbatch --parsable --export=ALL,SMOKE=0 --output="$outdir/logs/%x_%A_%a.log" "$@"; }
all=()
while read -r pt seed0; do
    case "$pt" in
        f10_p*) m=$(sub --time=02:00:00 $job meanfield "$pt" "$outdir"); all+=("$m")
                all+=("$(sub --array=0-0 --time=01:00:00 $job ensemble "$pt" full 1 40 "$seed0" "$outdir")") ;;
        *)      all+=("$(sub --array=0-4 --time=01:30:00 $job ensemble "$pt" full 5 200 "$seed0" "$outdir")") ;;
    esac
done < <(julia --project dev/scripts/jensen_rethink_89.jl seeds)
deps=$(IFS=:; echo "${all[*]}")
for st in transfer calibrate score; do
    all+=("$(sbatch --parsable --export=ALL,SMOKE=0 --output="$outdir/logs/%x_${st}_%j.log" \
        --time=01:00:00 --mem=16G --dependency=afterany:$deps $job $st "$outdir")")
done
echo "jobs ${all[*]}"
