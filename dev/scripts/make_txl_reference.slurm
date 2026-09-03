#!/bin/bash
#SBATCH --job-name=infercell_txl_ref
#SBATCH --partition=milan
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=01:00:00
#SBATCH --output=dev/scripts/txl_ref_%j.log

# Regenerates test/fixtures/txl_reference.jl on the unchanged tree and shows
# the result is deterministic: the script runs twice in two processes and the
# two files must be byte-identical.
module load julia/1.10.5

echo "Generating TXL byte-identity reference at $(date) on $(hostname)"
echo "HEAD: $(git rev-parse --short HEAD)"

OUT=test/fixtures/txl_reference.jl
TMP=$(mktemp -p dev/scripts txl_ref_check_XXXX.jl)

julia --project dev/scripts/make_txl_reference.jl "$OUT" || exit 1
julia --project dev/scripts/make_txl_reference.jl "$TMP" || exit 1

# The header carries a timestamp; compare everything after it.
if cmp <(grep -v '^# Julia' "$OUT") <(grep -v '^# Julia' "$TMP"); then
    echo "DETERMINISTIC: two independent runs produced identical bit patterns"
else
    echo "NOT DETERMINISTIC: the two runs differ" ; diff "$OUT" "$TMP"; rm -f "$TMP"; exit 2
fi
rm -f "$TMP"
echo "Finished at $(date)"
