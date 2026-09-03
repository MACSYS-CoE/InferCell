#!/bin/bash
#SBATCH --job-name=infercell_ssa_ref
#SBATCH --partition=milan
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=01:00:00
#SBATCH --output=dev/scripts/ssa_ref_%j.log

# Regenerates test/fixtures/ssa_reference.jl on the checked-out tree and shows
# the result is deterministic: the script runs twice in two processes and the
# two files must be byte-identical. The committed header names a pre-phase-2
# commit because the generator was run against a worktree at that commit; see
# the comment at the top of make_ssa_reference.jl.
module load julia/1.10.5

echo "Generating seeded SSA reference at $(date) on $(hostname)"
echo "HEAD: $(git rev-parse --short HEAD)"

OUT=test/fixtures/ssa_reference.jl
TMP=$(mktemp -p dev/scripts ssa_ref_check_XXXX.jl)
trap 'rm -f "$TMP"' EXIT INT TERM

julia --project dev/scripts/make_ssa_reference.jl "$OUT" || exit 1
julia --project dev/scripts/make_ssa_reference.jl "$TMP" || exit 1

# The header carries a timestamp; compare everything after it.
if cmp <(grep -v '^# Julia' "$OUT") <(grep -v '^# Julia' "$TMP"); then
    echo "DETERMINISTIC: two independent runs produced identical bit patterns"
else
    echo "NOT DETERMINISTIC: the two runs differ" ; diff "$OUT" "$TMP"; exit 2
fi
echo "Finished at $(date)"
