# InferCell's evidence archive (KIM-43)

The outputs of M5's runs, and the inputs they read. These are M0's crossed comparison and ptsG
panel (KIM-40, KIM-41), scored against a truth, and the reduced model's budget (KIM-42), which
is **unscored**. The budget is not a calibration result.

Every file here is a copy. The result notes cite the original directories, which stay in place.

## Checking it

    cd <this directory> && sha256sum -c --quiet SHA256SUMS

`SHA256SUMS` covers every file except itself and `build.log`, the build job's own record.
`jobs.tsv` gives that job's id and the commit it ran at. The script is
`dev/scripts/archive_predict_kim43.sh` at that commit. Each directory in `runs/` also keeps the
`SHA256SUMS` written when it ran, except `secheck-kim42`, which had none: its `SHA256SUMS` was
written for this archive, after a byte-for-byte comparison with the source.

## What is here

`PROVENANCE.tsv` gives each directory's issue, commit, job ids, original path and result note.

| directory | what | files |
|---|---|---|
| `runs/m0-kim40-kim41/` | M0's run, its score and its audit | 3,663 |
| `runs/budget-kim42/` | the budget at widths ½, 1 and 2 | 19 |
| `runs/secheck-kim42/` | KIM-42's SE check | 5 |
| `runs/checks-kim42/` | login-node checks behind KIM-42's derived figures; no job ids | 28 |
| `runs/pilot-kim35/` | the sensitivity pilot KIM-42's note compares with | 17 |
| `inputs/v7-fits/` | V7's 50 exact fits; M0 reads all of them, and KIM-42 reads `replicate1.jls` | 50 |
| `inputs/b7-chains/` | B7's mixed fits, `r1` to `r50`, 3 files each; M0 reads them | 150 |
| `inputs/b7-v6/` | the `v6` chains that KIM-35 and KIM-42's checks read | 3 |
| `inputs/inputs.sha256` | the pins for all 203 inputs, with paths relative to `inputs/` | 1 |
| `inputs/earlier_archives.sha256` | the `SHA256SUMS` hash of each earlier directory not copied here | 1 |

`inputs/inputs.sha256` carries M0's pins, written before its run, and B7's own checksums for the
`v6` chains. `runs/m0-kim40-kim41/inputs.sha256` is unchanged, so it still names the original
paths.

## Where the headline numbers are

| number | file |
|---|---|
| Q's coverage, width and z, per arm | `runs/m0-kim40-kim41/score.md`, `score_datasets.tsv` |
| ΔQ's coverage, width and z, per arm | the same files |
| the 88.2% reference for coverage at nominal 90% | `dev/scripts/m0_evidence_kim40_result.md`, "How to read the coverage"; the note derives it, and `score.md` does not state it |
| held share 0.631, 0.908 and 0.980 at widths ½, 1 and 2 | `runs/budget-kim42/budget_rows.tsv`, `budget.md` |
| the budget's derived figures | `runs/budget-kim42/w*.jls`, `plain_w*.jls`, `prediction_w*.jls`, read by `runs/checks-kim42/` |
| the SE check's readings | `runs/secheck-kim42/secheck.tsv` |
| cost, nodes and memory | each run's `sacct.tsv` and `logs/` |

The KIM-42 note's "Where the numbers come from" and the M0 note's run section go further. The
`.jls` files open only at the commit their run used; see `PROVENANCE.tsv`.

## What this archive does not hold

- **The paper's §5.5.** It had not been written when this was built. KIM-46 writes it, and is
  where its numbers are traced to the files above.
- **The earlier results** the paper cites, such as V6, V7's ranks, F10 and the Jensen reports,
  and B7's merged tables. Those stay in their own archives. Only the fits M5 read are copied.
- **The M4 pilots** other than KIM-35. `inputs/earlier_archives.sha256` pins them.
- **Any published copy.** This sits on the cluster's filesystem. Public deposit is KIM-54's.
