# KIM-95 archive: scored predictions on the reduced model

The archive is `/fred/oz022/tkimpson/infercell-kim95-archive-2026-11-19/`. It is a copy of
the run outdir `/fred/oz022/tkimpson/infercell-kim95-2026-11-18/`, plus this README, and
`SHA256SUMS` covers every file. The result note is
`dev/scripts/corea_predict_kim95_result.md`.

## Code and jobs

| stage | job | commit | what |
|---|---|---|---|
| precompile | 21122608 | 473ecf8 | `using InferCell` |
| prep | 21122609 | 473ecf8 | `prep.jls`, `prep.md` |
| sizing | 21122610 | 473ecf8 | 40 draws × 8 cells per arm: `sizing.*` |
| resolve | 21126961 | 3b02ea4 | 40 draws × 100 cells, posterior only: `sizing_resolve.*` |
| truth | 21154497 (array 1–4) | baf5331 | 10,810 cells at θ* (`rows_truth_*.jls`) |
| posterior | 21154498 (array 1–20) | baf5331 | 100 draws × 973 cells (`rows_posterior_*.jls`) |
| prior | 21154499 (array 1–4) | baf5331 | 100 draws × 200 cells (`rows_prior_*.jls`) |
| collect + score | 21154500 | baf5331 | the TSVs, `collect.md`, a first `score.md` |
| score | login node | 2333ae2 | `score.md` and `score_targets.tsv` as archived, with each arm's MC share added |

All jobs ended COMPLETED 0:0, as `sacct.tsv` records. Julia came from the cluster's
`julia/1.10.5`. Logs are in `logs/`.

## Files

**Inputs.**
- **`prep.jls`** holds three sets of draws:
  - the posterior: KIM-92's 1600 draws, from `infercell-kim92-archive-2026-11-17`, checked to the bit against its `pilot_draws.tsv`;
  - the prior: 1600 draws on seed 17950001;
  - θ*: the truth, from dataset seed 17200001.
- **`prep.md`** records the same.

**Sizing.**
- `sizing.{jls,tsv,md}` and `sizing_resolve.{jls,tsv,md}` read only the predictive, never the truth.

**Per-chunk rows.**
- **`rows_<arm>_<k>.jls`** holds, for each draw and target, μ̂ (the mean over finished cells), σ̂² and n.
- **Seeds:** the cell seeds, which `collect` checks are unique and disjoint from the dataset's and the reference ensemble's.

**Plain TSVs** (the scorer reads only these):
- **`targets.tsv`:** id, family, species and t for each of the 1,805 targets:
  - Q, Q_kd and dQ;
  - 17 pools and 17 transcripts × 53 saves, at 3180–6300 s.
- **`predictive_<arm>_{mu,s2,n}.tsv`:** one row per draw (`draw`, `chain`, `iteration`), then one column per target, at %.17g.
- **`truth.tsv`:** the truth's mean, SD, n and SE for each target, pooled over its 4 chunks.
- **`status.tsv`:** cell statuses, by arm and by intervention.

**Scores.**
- `score.md` and `score_targets.tsv`, written by `dev/scripts/corea_predict_kim95_score.jl`.

## Reproduce the scores

```bash
cp -r /fred/oz022/tkimpson/infercell-kim95-archive-2026-11-19 /tmp/k95 && cd <repo at 2333ae2 or later>
julia --project dev/scripts/corea_predict_kim95_score.jl /tmp/k95
diff /tmp/k95/score.md /fred/oz022/tkimpson/infercell-kim95-archive-2026-11-19/score.md
diff /tmp/k95/score_targets.tsv /fred/oz022/tkimpson/infercell-kim95-archive-2026-11-19/score_targets.tsv
```

The scorer needs only `Statistics`, `Printf` and `Random`. It runs its CRPS self-checks first.
