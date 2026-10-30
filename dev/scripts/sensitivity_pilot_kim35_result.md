# KIM-35: lactate export's sensitivity to θ_I and θ_H (M4 pilot)

**Run:** `dev/scripts/sensitivity_pilot_kim35.sh` at `a742224` (pushed), into
`/fred/oz022/tkimpson/infercell-m4-kim35-2026-10-30/` (`SHA256SUMS`, `jobs.tsv`).
- Jobs: precompile 19560137; arms 19560138 (I, 6,300 s), 19560139 (H, 6,300 s),
  19560140 (I, 600 s), 19560141 (H, 600 s); local sensitivities 19560142
  (6,300 s) and 19560143 (600 s); merge 19560144.
- The merge's table is `sensitivity.md` in that directory.

**Held out:** θ_I comes from B7's `v6` fit (m0_v6, not one of the 50 evidence
datasets), and θ_H from its prior. Every rng seed is from M4's block (the
script's header lists them).
- **The v6 chains** are
  `/fred/oz022/tkimpson/infercell-89-b7-archive-2026-10-27/run/v6/chain_{1,2}.jls`
  (`dev/notes/2026-10-29-archived-posteriors.md`). Their SHA256s are
  `8e1a2d2e…bbe82cd` and `7103a478…4ce4e5`. A copy elsewhere means editing the
  script's `V6`.
- **The archive's `SHA256SUMS`** has SHA256 `8bb49b92ad4ba969…`.
- **The `.jls` files open only at `a742224`.** The arms store whole
  `Prediction`s, and KIM-36 added fields to that type afterwards. So the merge
  stage reruns from a worktree at `a742224`. `sensitivity.md` carries every
  number below.

**Design and rule** (declared in the script's header before the run):
- **Arm I:** 50 v6 draws × 4 cells, with θ_H nominal. The rows are independent,
  so var − MC is Var_I.
- **Arm H:** θ_I at v6's posterior mean of log θ, with 50 prior draws of the 15
  held promoters × 4 cells. var − MC is Var_H.
- **Pass:** each group's net variance is above 2 bootstrap SEs, and at least 5%
  of Var_I + Var_H.
- **How a failure is read:** a failure on θ_I alone means "change Q" only if
  θ_I's local elasticities are under a quarter of θ_H's. Otherwise it reflects
  the posterior's width.
  - **When it was written:** after the smoke pass, before this run.
  - **What the smoke had already shown:** job 19551114 at `657081e` (tag
    `m4-smoke-657081e`), 3 draws × 2 cells at 630 s and 60 s; its
    `kim35/sensitivity.md` has SHA256 `9ff2684a…`. It printed:
    - θ_I's share of `lactate_out` as 0 at both horizons. That was noise:
      Var_I was negative net of the MC term;
    - deterministically, mean-field elasticities of the same size for the two
      inferred and two held parameters it sampled: 0.037 and 0.020, against
      0.013 and 0.021;
    - its pass column at the long horizon, "no (I false, H true)": a failure on
      θ_I alone;
    - the log SDs, 0.11–0.19 for θ_I against 0.69.

    The file labels its sections with the nominal H (6300 s and 600 s), a label
    fixed at `a742224`.
  - **So the reading and its quarter threshold were written knowing roughly how
    they would come out.** The approved rule (`dev/notes/2026-10-30-m4-plan.md`)
    is the blind one. Under it alone, this note would recommend changing Q.

## Results

| H | Q | Var_I ± SE | Var_H ± SE | share I | share H | rule |
|---|---|---|---|---|---|---|
| 6,300 s | `lactate_out` | 3.26e10 ± 1.3e10 (2.5 SE) | 9.70e11 ± 2.3e11 | 0.033 | 0.967 | fails on θ_I's share alone |
| 6,300 s | `exported` | −7.6e−8 ± 2.4e−8 | −1.5e−8 ± 4.2e−8 | — | — | fails on both |
| 600 s | `lactate_out` | 1.07e9 ± 9.2e8 (1.2 SE) | 6.9e9 ± 4.0e9 (1.7 SE) | 0.134 | 0.866 | fails on both (neither above 2 SE) |
| 600 s | `exported` | 2.6e−7 ± 8.3e−7 | 4.4e−6 ± 1.9e−6 | 0.056 | 0.944 | fails on θ_I |

- **θ_I's share is 3.3% ± 1.5%** (the delta method on the two bootstrap SEs).
  Its 2-SE band reaches past the rule's 5%, so the failure is itself not
  resolved.
- **The MC term is as large as Var_I.** At 6,300 s it is 3.29e10 against
  Var_I's 3.26e10, with 4 cells per pair.
- **Mean-field local elasticities** (d ln Q / d ln θ at the reference, 6,300 s):
  - inferred: `krnadeg` −0.62, the largest of all 20; S_0607 0.085; S_0779 0.056;
    ENO and FBA about 0. The median |ε| is 0.056;
  - held: 0.011 to 0.128, with S_0213 0.128 and S_0475 0.112. The median |ε| is
    0.047.

  θ_I's elasticities are therefore not small next to θ_H's. Weighted by each
  parameter's log SD, the local variance shares are about 7% for θ_I and 93% for
  θ_H. That agrees in order with the simulated 3.3%.
- **`exported` is nearly constant:** 0.9985 at 6,300 s and 0.992 at 600 s, with
  every |ε| ≤ 5e−3.
- **No cell starved.** At 6,300 s, 7 and 11 of 200 cells clipped, in arms I and H.
  Under the candidate rules (clipping allowed, starved at 0.5, no floors), every
  pair is viable, whether min_viable_fraction is 0.5 or 1.0.
  - The width-1 evidence is arm H alone: one θ_I × 50 θ_H, and zero of 50 has a
    95% upper bound of about 6%.
  - Arm I's θ_H are nominal.

## Recommendation (the freeze decides)

- **By the reading, keep Q = `lactate_out`, at the reduced-model horizon of
  6,300 s.** By the approved rule alone, change it. The user adopted the reading
  on 2026-10-30. At that gate they were told it had been added after the plan's
  approval, but not what the smoke had previewed. That came out in the PR #123
  review, after which the user re-confirmed keeping `lactate_out`, disclosed.
  - Both groups move it, at comparable elasticities.
  - θ_I's small share comes from its posterior's width (log SDs 0.05–0.19,
    against the prior's 0.69).
  - θ_I's variance is resolved, at 2.5 SE. That is the split the paper reports.
- **Drop `exported`.** It barely moves with any parameter.
- **Do not use 600 s on the reduced model.** At 4 cells, neither group's
  variance is resolved at 2 bootstrap SE.
