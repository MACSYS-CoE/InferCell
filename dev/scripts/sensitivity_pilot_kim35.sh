#!/bin/bash
# KIM-35: lactate export's sensitivity to θ_I and θ_H
# (dev/scripts/sensitivity_pilot_kim35.jl), from a worktree at a pushed commit.
#   dev/scripts/sensitivity_pilot_kim35.sh /fred/oz022/tkimpson/infercell-m4-kim35-<date>
# One job precompiles on milan; then the four arms (16 threads each), the two
# local-sensitivity stages, and the merge once all six finish. Job ids go to
# <outdir>/jobs.tsv, logs to <outdir>/logs.
set -euo pipefail
outdir=$(realpath -m "$1")
case "$outdir" in "$(git rev-parse --show-toplevel)"/*|*/.worktrees/*|*/.claude/worktrees/*)
    echo "outdir $outdir is inside a worktree" >&2; exit 1;; esac
if [ -n "$(git status --porcelain -- src test dev/scripts Project.toml Manifest.toml)" ]; then
    echo "the tree is dirty; run from a clean worktree at a pushed commit" >&2; exit 1
fi
if [ -z "$(git branch -r --contains HEAD)" ]; then
    echo "HEAD $(git rev-parse --short HEAD) is on no remote branch; push it first" >&2; exit 1
fi
mkdir -p "$outdir/logs"
sub() { sbatch --parsable --partition=milan --exclude='dave[148-149]' --ntasks=1 \
               --output="$outdir/logs/%x_%j.log" "$@"; }
run="module load julia/1.10.5; julia --project"
: > "$outdir/jobs.tsv"
pre=$(sub --cpus-per-task=1 --mem=16G --time=01:00:00 --job-name=kim35_pre \
          --wrap "$run -e 'using InferCell'")
printf 'precompile\t%s\n' "$pre" >> "$outdir/jobs.tsv"
deps=""
for stage in arm_I_6300 arm_H_6300 arm_I_600 arm_H_600; do
    j=$(sub --dependency=afterok:$pre --cpus-per-task=16 --mem=48G --time=03:00:00 \
            --job-name=kim35_$stage \
            --wrap "$run --threads=16 dev/scripts/sensitivity_pilot_kim35.jl $stage $outdir")
    printf '%s\t%s\n' "$stage" "$j" >> "$outdir/jobs.tsv"; deps="$deps:$j"
done
for stage in local_6300 local_600; do
    j=$(sub --dependency=afterok:$pre --cpus-per-task=1 --mem=16G --time=06:00:00 \
            --job-name=kim35_$stage \
            --wrap "$run dev/scripts/sensitivity_pilot_kim35.jl $stage $outdir")
    printf '%s\t%s\n' "$stage" "$j" >> "$outdir/jobs.tsv"; deps="$deps:$j"
done
j=$(sub --dependency=afterok$deps --cpus-per-task=1 --mem=16G --time=00:30:00 \
        --job-name=kim35_merge \
        --wrap "$run dev/scripts/sensitivity_pilot_kim35.jl merge $outdir")
printf 'merge\t%s\n' "$j" >> "$outdir/jobs.tsv"
echo "submitted at $(git rev-parse --short HEAD):"
cat "$outdir/jobs.tsv"
