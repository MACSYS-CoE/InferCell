#!/bin/bash
# Continue bulk V4 from a grid stage s whose box the previous merge wrote
# (task 16a.7c): grid s and merge s, then the FBA-fine and ENO-fine grids, the
# chains and the target check, then the final test. Usage: <s>. Run from a
# worktree at a pushed commit holding the earlier stages' files.
set -euo pipefail
s=$1
D=dev/scripts/block1_16a5_bulk
n=$(tr ',' '\n' < $D/rows$s.txt | wc -l)
S=dev/scripts/block1_16a5.slurm
E="--export=ALL,V4_BULK=1"
g=$(sbatch --parsable $E --time=01:00:00 --array=0-$((n - 1)) -J v4b_grid$s $S grid $s)
m=$(sbatch --parsable $E --time=00:20:00 -J v4b_merge$s --dependency=afterok:$g $S merge $s)
fb=$(sbatch --parsable $E --time=01:30:00 --array=0-40 -J v4b_fba --dependency=afterok:$m $S fba)
fm=$(sbatch --parsable $E --time=00:20:00 -J v4b_fbamerge --dependency=afterok:$fb $S fbamerge)
en=$(sbatch --parsable $E --time=01:30:00 --array=0-40 -J v4b_eno --dependency=afterok:$m $S eno)
em=$(sbatch --parsable $E --time=00:20:00 -J v4b_enomerge --dependency=afterok:$en $S enomerge)
ch=$(sbatch --parsable $E --time=06:00:00 --array=0-19 -J v4b_chain --dependency=afterok:$m $S chain)
ck=$(sbatch --parsable $E --time=00:30:00 -J v4b_check --dependency=afterok:$m $S check)
fi=$(sbatch --parsable $E --time=00:20:00 -J v4b_final --dependency=afterok:$ch:$fm:$em $S final)
echo "grid$s $g, merge$s $m, fba $fb/$fm, eno $en/$em, chains $ch, check $ck, final $fi"
