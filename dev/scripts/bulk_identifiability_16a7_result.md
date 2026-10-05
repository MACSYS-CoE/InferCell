# Bulk identifiability: is F6 possible with population-average metabolites?

Job 17934962 at `64ecf5f`, `dev/scripts/bulk_identifiability_16a7.jl`, from 15.8's
sensitivity run (job 17704779 at `7835c20`: 1470 seeds, the sampler model, D11's
six freed, δ = 0.2). No new simulation.
- **Each row** is an observable's 200-cell ensemble mean at one save point, in units
  of its resolution.
- **Metabolite rows** are floored at a bulk assay CV: 1% (15.8's own), 5%, 10% or
  20%.
- **Transcript rows** keep 15.8's resolution, the standard error of a 200-cell mean
  floored at 1%.
- **Each cell** of the tables below is F6's quantity, posterior over prior SD,
  linearised: Fisher information JᵀJ plus the lognormal prior. In brackets is the
  same ratio with the split-half noise matrix in place of J. A cell no lower than
  its bracket has not shrunk.
- **Save points are treated as independent.** "Sparse" keeps every 15th, 7 save
  points, to show how much rests on that.

## F6's two load-bearing cells

| Cell | CV 1% | CV 5% | CV 10% | CV 20% |
|---|---|---|---|---|
| ptsG promoter, metabolites only | 0.02 (0.05) | 0.05 (0.10) | 0.10 (0.19) | 0.19 (0.36) |
| The same, sparse | 0.08 (0.19) | 0.18 (0.35) | 0.34 (0.59) | 0.57 (0.81) |
| ENO, transcripts only | 0.17 (0.18) | the same | the same | the same |
| The same, sparse | 0.59 (0.56) | the same | the same | the same |

- **The ptsG cell survives population averaging.** At a 10% or 20% assay CV it
  shrinks to about half its noise bracket, sparse or not. Metabolites-only rank is
  6 of 6 at CV 5 to 20% (5 of 6 at 1% and when sparse).
- **The ENO cell is not resolved, at any CV.** Transcript rows do not depend on the
  metabolite model, so this is unchanged from 15.8's own run.
  - Transcripts-only rank is 4 of 6, and ENO's ratio sits at its noise bracket.
  - 15.8 tested only the joint rank, so this went unseen.
  - The transcript rows' Monte Carlo noise (floor 54.6) is as large as ENO's
    signal, so 15.8's run cannot tell "no information" from "too few seeds".
  - Ensemble means also discard what exact per-cell counts carry. So the ENO cell
    is **unresolved, not refuted**.

## The full tables

From the log, per CV, rows: metabolites only, transcripts only, joint, metabolites
only (sparse), transcripts only (sparse). The column order is S_0607, S_0445,
S_0779, krnadeg, ENO, FBA.

| CV | Design | Rank | S_0607 | S_0445 | S_0779 | krnadeg | ENO | FBA |
|---|---|---|---|---|---|---|---|---|
| 1% | metabolites | 5 | 0.02 (0.04) | 0.03 (0.04) | 0.02 (0.05) | 0.01 (0.04) | 0.03 (0.26) | 0.00 (0.14) |
| 1% | transcripts | 4 | 0.01 (0.03) | 0.01 (0.03) | 0.01 (0.03) | 0.00 (0.03) | 0.17 (0.18) | 0.10 (0.11) |
| 1% | joint | 6 | 0.01 (0.02) | 0.01 (0.02) | 0.01 (0.03) | 0.00 (0.02) | 0.03 (0.15) | 0.00 (0.08) |
| 10% | metabolites | 6 | 0.06 (0.18) | 0.14 (0.16) | 0.10 (0.19) | 0.03 (0.16) | 0.13 (0.75) | 0.04 (0.51) |
| 10% | joint | 5 | 0.01 (0.03) | 0.01 (0.03) | 0.01 (0.03) | 0.00 (0.03) | 0.10 (0.18) | 0.04 (0.10) |
| 20% | metabolites | 6 | 0.13 (0.34) | 0.28 (0.30) | 0.19 (0.36) | 0.05 (0.30) | 0.26 (0.91) | 0.08 (0.76) |
| 20% | joint | 4 | 0.01 (0.03) | 0.01 (0.03) | 0.01 (0.03) | 0.00 (0.03) | 0.14 (0.18) | 0.06 (0.11) |

The 5% rows and the sparse rows are in the log,
`dev/scripts/bulk_identifiability_16a7_17934962.log` in the run's worktree (logs are ignored; the job id is the record).

**Under bulk metabolites, ENO is identified from the metabolites.** It shrinks to
0.13 at 10% CV, against a bracket of 0.75. So the seam claim's second direction,
transcripts informing ENO, is a separate question from the observation model.
