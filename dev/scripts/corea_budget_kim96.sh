#!/bin/bash
# KIM-96: run dev/scripts/corea_budget_kim96.jl as one short Slurm job, from a clean worktree
# at a pushed commit, so that its numbers carry a job id.
#   dev/scripts/corea_budget_kim96.sh <outdir>
# The job writes budget.md, budget_rows.tsv and inputs.tsv to <outdir>, then SHA256SUMS over
# them and over a copy of the script that ran. Its id goes to <outdir>/jobs.tsv, its log to
# <outdir>/logs. It reads only archived TSVs and simulates nothing.
set -euo pipefail
outdir=$(realpath -m "${1:?give the output directory}")
[ -e "$outdir/budget.md" ] && { echo "$outdir/budget.md already exists" >&2; exit 1; }
source "$(dirname "$0")/launcher_guard.sh"
script=dev/scripts/corea_budget_kim96.jl
commit=$(git rev-parse HEAD)
cp "$script" "$outdir/corea_budget_kim96.jl"
printf '%s\n' "$commit" > "$outdir/COMMIT"
j=$(sub --cpus-per-task=1 --mem=8G --time=00:30:00 --job-name=kim96_budget \
        --wrap "$run $script $outdir && cd $outdir && sha256sum budget.md budget_rows.tsv inputs.tsv corea_budget_kim96.jl COMMIT > SHA256SUMS")
record budget "$j"
echo "submitted at $(git rev-parse --short HEAD):"
cat "$outdir/jobs.tsv"
