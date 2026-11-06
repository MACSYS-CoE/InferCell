#!/bin/bash
# KIM-43: build the evidence archive for M5's runs (KIM-40, KIM-41, KIM-42).
#
# It copies and never moves, because the result notes cite the original paths. It checks every
# source before copying and every copy afterwards, and refuses a target that already exists.
#
# From a clean checkout at a pushed commit:
#   sbatch --chdir="$PWD" dev/scripts/archive_predict_kim43.sh
#SBATCH --job-name=kim43_archive
#SBATCH --partition=milan
#SBATCH --exclude=dave148,dave149,gina18
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=02:00:00
#SBATCH --output=/dev/null
set -euo pipefail

base=/fred/oz022/tkimpson
archive=$base/infercell-predict-archive-2026-11-06
m0=$base/infercell-m5-kim40-2026-11-01
budget=$base/infercell-m5-kim42-2026-11-05
secheck=$base/infercell-m5-kim42-secheck-2026-11-05
checks=$base/infercell-m5-kim42-checks-2026-11-06
kim35=$base/infercell-m4-kim35-2026-10-30
b7=$base/infercell-89-b7-archive-2026-10-27
v7=$base/infercell-16a9b-archive-2026-10-17/v7/sbc_v7_16a9b
here=$(pwd)

[ ! -e "$archive" ] || { echo "$archive exists; refusing" >&2; exit 1; }
[ -z "$(git -C "$here" status --porcelain)" ] || { echo "dirty checkout; refusing" >&2; exit 1; }
commit=$(git -C "$here" rev-parse HEAD)

mkdir "$archive"
exec > "$archive/build.log" 2>&1
echo "KIM-43 archive, job $SLURM_JOB_ID on $(hostname), at $commit, $(date -Is)"

echo "== the sources, before copying"
# The result notes quote these two hashes.
echo "29ca43c53b5eba874c9467c23107e121a0a480bf962115caedeaea49f1f8279a  $m0/SHA256SUMS" | sha256sum -c
echo "99d912fb68179769963a225b0c5122f5e20ce3e2573f4901f1f27c9d8e4576cd  $budget/SHA256SUMS" | sha256sum -c
(cd "$m0" && sha256sum -c --quiet SHA256SUMS && echo "m0: $(wc -l < SHA256SUMS) files pass")
(cd "$budget" && sha256sum -c --quiet SHA256SUMS && echo "budget: $(wc -l < SHA256SUMS) files pass")
(cd "$checks" && sha256sum -c --quiet SHA256SUMS && echo "checks: $(wc -l < SHA256SUMS) files pass")
(cd "$kim35" && sha256sum -c --quiet SHA256SUMS && echo "kim35: $(wc -l < SHA256SUMS) files pass")
sha256sum -c --quiet "$m0/inputs.sha256"
echo "M0's inputs: $(wc -l < "$m0/inputs.sha256") files pass"
(cd "$b7" && grep -E '  run/v6/(chain_[12]|prep)\.jls$' SHA256SUMS | sha256sum -c --quiet && echo "V6's chains pass B7's SHA256SUMS")

echo "== the runs, copied"
mkdir "$archive/runs"
cp -a "$m0" "$archive/runs/m0-kim40-kim41"
cp -a "$budget" "$archive/runs/budget-kim42"
cp -a "$secheck" "$archive/runs/secheck-kim42"
cp -a "$checks" "$archive/runs/checks-kim42"
cp -a "$kim35" "$archive/runs/pilot-kim35"
# The SE check's directory never had a SHA256SUMS. The new one is in the copy only.
(cd "$archive/runs/secheck-kim42" && find . -type f ! -name SHA256SUMS | LC_ALL=C sort | xargs -d '\n' sha256sum > SHA256SUMS)

echo "== the inputs, copied"
mkdir -p "$archive/inputs/v7-fits" "$archive/inputs/b7-chains" "$archive/inputs/b7-v6"
cp -a "$v7"/replicate*.jls "$archive/inputs/v7-fits/"
for r in $(seq 1 50); do
    mkdir "$archive/inputs/b7-chains/r$r"
    cp -a "$b7/run/r$r/chain_1.jls" "$b7/run/r$r/chain_2.jls" "$b7/run/r$r/prep.jls" "$archive/inputs/b7-chains/r$r/"
done
cp -a "$b7/run/v6/chain_1.jls" "$b7/run/v6/chain_2.jls" "$b7/run/v6/prep.jls" "$archive/inputs/b7-v6/"

# The pins are the ones written before the runs, with their paths rewritten to the copies.
{
    sed -e "s#  $v7/#  v7-fits/#" -e "s#  $b7/run/#  b7-chains/#" "$m0/inputs.sha256"
    grep -E '  run/v6/(chain_[12]|prep)\.jls$' "$b7/SHA256SUMS" | sed -e 's#  run/v6/#  b7-v6/#'
} > "$archive/inputs/inputs.sha256"

# Earlier directories this archive does not copy, pinned by their own SHA256SUMS.
sha256sum "$b7/SHA256SUMS" \
    "$base/infercell-m4-kim37-2026-10-30/SHA256SUMS" \
    "$base/infercell-m4-kim38-2026-10-30/SHA256SUMS" \
    "$base/infercell-m4-kim66-2026-10-30/SHA256SUMS" \
    "$base/infercell-m4-checks-2026-10-30/SHA256SUMS" \
    "$base/infercell-m4-smoke-2026-10-30/SHA256SUMS" \
    "$base/infercell-m5-kim68-2026-10-30/SHA256SUMS" > "$archive/inputs/earlier_archives.sha256"

echo "== the copies, against the sources' checksums"
(cd "$archive/runs/m0-kim40-kim41" && sha256sum -c --quiet SHA256SUMS && echo "m0: $(wc -l < SHA256SUMS) files pass")
(cd "$archive/runs/budget-kim42" && sha256sum -c --quiet SHA256SUMS && echo "budget: $(wc -l < SHA256SUMS) files pass")
(cd "$archive/runs/checks-kim42" && sha256sum -c --quiet SHA256SUMS && echo "checks: $(wc -l < SHA256SUMS) files pass")
(cd "$archive/runs/pilot-kim35" && sha256sum -c --quiet SHA256SUMS && echo "kim35: $(wc -l < SHA256SUMS) files pass")
# The SE check has no earlier checksums, so its copy is compared byte for byte.
diff -r --exclude=SHA256SUMS "$secheck" "$archive/runs/secheck-kim42"
echo "secheck: identical to the source"
(cd "$archive/inputs" && sha256sum -c --quiet inputs.sha256 && echo "inputs: $(wc -l < inputs.sha256) files pass")

echo "== the record"
cp -a "$here/dev/scripts/archive_predict_kim43_README.md" "$archive/README.md"
cp -a "$here/dev/scripts/archive_predict_kim43.tsv" "$archive/PROVENANCE.tsv"
printf 'commit\t%s\nbuild\t%s\n' "$commit" "$SLURM_JOB_ID" > "$archive/jobs.tsv"

echo "== the archive's own SHA256SUMS (all but itself and this log)"
cd "$archive"
find . -type f ! -path ./SHA256SUMS ! -path ./build.log | LC_ALL=C sort | xargs -d '\n' sha256sum > SHA256SUMS
sha256sum -c --quiet SHA256SUMS
echo "SHA256SUMS: $(wc -l < SHA256SUMS) files pass; its own hash: $(sha256sum SHA256SUMS | cut -d' ' -f1)"
echo "done, $(date -Is)"
