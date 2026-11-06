#!/bin/bash
# M0's after-run audit: ran_on records, inputs.sha256, then SHA256SUMS over the outputs.
# The job's cwd is the run tree, for `julia --project`; it writes only to the outdir.
#SBATCH --job-name=kim40_audit
#SBATCH --partition=milan
#SBATCH --exclude=dave148,dave149,gina18
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=03:00:00
set -u
out=$1
here=$2
module load julia/1.10.5
{
    echo "M0 after-run audit, job $SLURM_JOB_ID on $(hostname), $(date -Is)"
    echo "== ran_on, commit, SMOKE"
    OPENBLAS_NUM_THREADS=1 julia --project "$here/m0_audit.jl" "$out"; echo "audit exit: $?"
    echo "== sha256sum -c inputs.sha256"
    sha256sum -c --quiet "$out/inputs.sha256"; echo "inputs exit: $?, $(wc -l < "$out/inputs.sha256") files"
} > "$out/audit.txt" 2>&1

cd "$out" || exit 1
find . -path ./tmp -prune -o -type f ! -name SHA256SUMS -print | LC_ALL=C sort |
    xargs -d '\n' sha256sum > tmp/SHA256SUMS
mv tmp/SHA256SUMS SHA256SUMS
sha256sum -c --quiet SHA256SUMS; echo "SHA256SUMS exit: $?, $(wc -l < SHA256SUMS) files"
