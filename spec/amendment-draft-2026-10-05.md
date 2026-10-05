# Amendment, 2026-10-05: bulk metabolites, and a one-way seam (the reviewed text)

**Status: approved and applied 2026-10-05**, with every ⚠️ default taken. This file
is kept as the reviewed text. The record is each spec's §12 entry of 2026-10-05.
(It was drafted as: "draft for review. Nothing here is applied.") On approval, each change
below is made in place in `spec/spec.md` or `spec/phases/16-recovery.md`. Each
file gets a dated §12 entry, built from the "Trigger" section, and the "Last
amended" date moves. Items marked ⚠️ were choices, all taken at their defaults. After the #81 review,
the applied text also scopes the ENO null (one truth, linearised, 40 cells scaled
to 200), narrows "one way" to ENO with FBA open, states K2's relation, and moves N
from 10 to 20. See each spec's §12. Part B's "15.6 and 15.8 are re-measured" was
applied to 15.8 only. 15.6's panel choice stands, since the bulk check finds all
six targets identifiable from bulk metabolites.

## Trigger (the evidence, in order)

All in `dev/scripts/` on `phase-16a-sampler`:
1. **KP fired twice.** No path-update variant mixes at the data's σ: PG, PGAS and
   truncated PGAS at N = 50, and annealed windows at K ≤ 100.
   - The cause is the per-cell panel's near-empty pools (Pi, 3PG, 2PG) observed
     at σ ≈ 0.1. They pin a cell's history to a sliver.
   - Annealing that does reach it needs K ≥ 200, about 42k CPU-h per M0
     posterior. (`csmc_variants_16a7_result.md`)
2. **The per-cell panel is not measurable.** At the best cited single-cell
   detection limit, about 1.2 × 10⁵ molecules (Lin et al., *Anal. Chem.* 2011), 16
   of the 17 pools in a syn3A cell are below detection. At that floor the
   panel's information about a history falls to 0.00 ± 0.02 nats. (Same file,
   "Option 1".)
3. **A cut posterior at the history was rejected.** It forbids metabolites →
   gene expression by construction, for the reason Amendment 4 (2026-09-03)
   dropped the cut. The metabolic constants also sit in both modules.
4. **Bulk metabolites keep the ptsG cell.** From metabolites alone, ptsG's promoter
   shrinks to posterior over prior SD 0.10 at a 10% assay CV, against a noise
   bracket of 0.19, and 0.19 (0.36) at 20%. It holds on 7 save points.
   (`bulk_identifiability_16a7_result.md`)
5. **The ENO cell is refuted.** Even a cell's full path leaves ENO's posterior SD at
   99.9% of its prior over 200 cells (bootstrap 0.998 to 1.000). A 20% change in
   ENO moves a path's log-density by about 0.003 nats. Transcripts are a function
   of the path, so they carry no more. (`eno_path_information_16a7_result.md`)

## Part A — the claim: the seam carries information one way

**Parent §1, the "crosses the boundary" paragraph.** Its last sentence becomes:

> §6 F6 is where it is cashed out. In the metabolites → gene-expression direction,
> it is the ptsG promoter strength shrinking under metabolite-only data. In the
> other direction, the enolase constant under transcript-only data is shown
> **not** to shrink: a complete-data information bound (§12 2026-10-05) leaves its
> posterior SD at 99.9% of its prior over 200 cells. The pools reach the
> stochastic block only through rate constants rebuilt once a minute. That
> direction is reported as a measured null, not omitted.

**Parent §1, the bar table.** The last row's right column becomes "Uncertainty
crosses the boundary correctly where information crosses it, and the information
in each direction is measured."

**Parent §0, "The result the work exists for."** The two load-bearing cells become
one load-bearing cell and one measured null, worded as above.

**Parent §6 F6.**
- **The table is unchanged:** six targets × three data configurations.
- **The "number it delivers" column becomes:** "One load-bearing cell: the ptsG
  promoter under metabolites only (shrinks ⟹ information crossed from the ODE
  block's data into the stochastic block's parameter). One measured null: ENO under
  transcripts only, which the information bound says cannot shrink. A shrinkage
  there would mean a sampler error, so it is a check, not a hope."

**Parent §6, the abstract scalars.** "The tight-prior control's shrinkage under
transcript-only data" is replaced by "the information bound on ENO from
transcripts".

**Parent K4** (the tight-prior control, scored in 16b.1 and 16b.3) is **unchanged.**
It tests recovery under a tight prior, not transcript-only shrinkage.

## Part B — the observation model: bulk metabolites, joint inference

**What is observed.**
- **Transcripts: unchanged.** Exact per-cell counts every 60 s.
- **Metabolites: one bulk measurement per pool per save time,** of the population
  mean over the dataset's C cells, lognormal at a bulk assay CV σ_b:

  ```
  log z_jk  ~  N( log max( (1/C) Σ_c x_cjk , floor ), σ_b )
  ```

  `x_cjk` is the carry-inclusive particle count, as now. The floor stays one
  particle, which a 200-cell mean never approaches.
- ⚠️ **σ_b:** default truth 0.10, prior `LogNormal(log 0.2, 1)`, the same prior as
  σ now. The bulk check covers 5 to 20%.
- ⚠️ **Which cells are averaged:** default, the same C cells at every save, which is
  synthetic and idealised. The alternative is a fresh cohort per save, which is
  realistic but makes the data unlinkable to a cell's history and changes the
  target. The default is stated as an idealisation in T2.
- ⚠️ **Particle count or concentration:** default, the mean particle count. Volume
  grows about 1.07× over the cycle, so the two differ by under 7%. The choice is
  stated in T2.

**The target** (sub-spec §3). The per-cell metabolite product is replaced by

```
Π_k Π_j  N( log z_jk ; log max( x̄_jk(θ, X_1:C), 1 ), σ_b )   with   x̄_jk = (1/C) Σ_c x_cjk
```

**Cells are no longer conditionally independent given θ.** They couple through x̄.

**The three blocks** (sub-spec §3 table, and D16.3).
- **B1:** unchanged in cost. Replay every cell, form x̄, score it. σ_b's conditional
  needs no replay.
- **B2:** unchanged.
- **B3:** per cell, conditional on the other cells' current histories.
  - Cell c's window weight uses the bulk likelihood at x̄ with cell c's particle
    substituted for its current contribution: `x̄ = (S_{−c,k} + x_ck)/C`.
  - **Cells are updated one at a time** in a sweep, not in parallel. This is exact,
    and it costs wall-clock, not CPU-h.
  - One cell moves x̄ by about 1/C of its deviation, so the weight is gentle. The
    expectation is that PG mixes, as the metabolite-off runs did (97 to 100% of
    windows changed). That expectation is checked, not assumed (Part D, 16a.7c).
  - The annealed window (§12 2026-10-02) stays in the code, off by default.

**Phase 15** (parent §11).
- **15.6 and 15.8 are re-measured under bulk data,** by the bulk check (item 4). At
  5 to 20% CV, metabolites-only rank is 6 of 6, and joint rank is 6 at 5%, 5 at
  10% and 4 at 20%.
  - ⚠️ The joint rank falling to 4 or 5 at higher CV is the 1%-floored transcript
    rows' Monte Carlo noise dominating the floor (54.6), not lost information.
  - The amendment records it, and annotates 15.8 rather than re-opening it.
- **15.7's dataset gains bulk observations.**
  - The latent record is regenerable (about 1.8 CPU-h for 200 cells), and the bulk
    z comes from it with a recorded noise seed.
  - The per-cell panel stays in the record as the old observation model, for
    reference.
  - M0's `m0_dataset` gains the same bulk modality.

## Part C — the sub-spec's checks

- **V5's metabolite-on half,** against brute force on a small population. On the toy
  and on M0 cut to ptsG, C = 3 cells are simulated jointly. The rejection reference
  keeps runs matching all three cells' transcripts, weighted by the bulk
  likelihood. The kernel is a sweep updating the three cells in turn. The tolerance
  is unchanged: within 3 SE, with σ_b inflated only if the ESS needs it.
- **V7 (SBC):** unchanged in form. σ_b replaces σ among the drawn quantities.
- **F10 and V6:** unchanged.
- **D16.6's cost projection:** PG at about N cycles per cell per sweep, as first
  projected. Sequential cells multiply wall-clock by C on one core; CPU-h is
  unchanged.

## Part D — tasks

- **16a.7:** annotated. Variant chosen: **PG** (lag 0), under bulk data, at an N set
  by 16a.7c's rates.
- **16a.7a:** annotated. It is superseded for production by the bulk model; its
  code and V5 pass stand.
- **New 16a.7c:**
  - Implement the bulk modality (dataset, B1, B3 with sequential cells).
  - Pass V5 in both cases under it.
  - Measure per-window update rates and cost per sweep on M0 at N ∈ {5, 10, 20, 50}.
  - Choose N, or record KP.
- **Phase 15:** 15.7 gains the bulk observations (a new task, 15.7b), and 15.8 is
  annotated with the bulk check.
- **16b.4 (F6):** the ENO cell is reported against the bound, as Part A says.

## What this amendment does not do

- It does not make the transcript data realistic. Exact counts every 60 s from one
  living cell are also an idealisation. T2 states it, and nothing here depends on
  changing it.
- It does not test the ENO bound at another truth. One truth is a stated limit.
  It is cheap to repeat (about 5 CPU-h a truth) if you want a second.
