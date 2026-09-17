#!/bin/bash
#SBATCH --job-name=extract_tx_genes
#SBATCH --partition=milan
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=00:20:00
#SBATCH --output=dev/scripts/extract_tx_genes_%j.log

# Regenerate the per-gene transcription extract and check it round-trips.
#
# The extract must reproduce byte for byte from the same checkout (spec §5), so
# this generates twice into two paths and compares. Unlike the SSA and TX/TL
# reference fixtures there is no timestamp to filter out — the header carries
# the upstream commit and nothing else that varies — so the comparison is over
# the whole file.

module load julia/1.10.5

CHECKOUT=${1:-/fred/oz022/tkimpson/Minimal_Cell}
OUT=src/organisms/coreA/data/transcription_genes.tsv
# Not under dev/scripts: a hard kill (Slurm's SIGKILL after the grace
# period, an OOM, a node failure) skips the trap, and a leftover file
# there is byte-identical to the real extract sitting next to it.
TMP=$(mktemp -t tx_genes_check_XXXX.tsv)
trap 'rm -f "$TMP"' EXIT INT TERM

echo "Regenerating $OUT from $CHECKOUT at $(date)"
echo "Upstream HEAD: $(git -C "$CHECKOUT" rev-parse HEAD)"
echo ""

julia --project=dev/scripts dev/scripts/extract_transcription_genes.jl "$CHECKOUT" "$OUT" || exit 1
julia --project=dev/scripts dev/scripts/extract_transcription_genes.jl "$CHECKOUT" "$TMP" || exit 1

if cmp "$OUT" "$TMP"; then
    echo "DETERMINISTIC: two independent runs produced identical bit patterns"
else
    echo "NOT DETERMINISTIC: the two runs differ"; diff "$OUT" "$TMP"; exit 2
fi

echo ""
echo "Finished at $(date)"
