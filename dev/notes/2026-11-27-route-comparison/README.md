# KIM-195: the candidate routes for the full cell's block, compared on one frame

Linear: [KIM-195](https://linear.app/kimpson/issue/KIM-195/compare-the-candidate-routes-for-the-full-cells-block-what-each).
Feeds the route decision, KIM-187, and T1. The frame (the cases each route is costed
on) is in [`frame.md`](frame.md), committed at d465a35 before any timing in this issue ran.

Status: in progress. Sections are added one checklist item at a time.

## 1. The candidates: what each approximates, and where its error is measured

Five routes are compared. The first three are KIM-187's original options. The last two
were added on 2026-11-27 (Decision log).

| # | route | in one line |
|---|---|---|
| R1 | **mean field + NUTS** (the mixed route carried over) | the bulk factor's mean from a deterministic mean field of the whole cell; the transcript factor reads its rates from the mean-field pools; NUTS on the block |
| R2 | **data plug-in** | stage 1: the transcript factor with its rates taken from the observed panel's pools, so no ODE is solved; stage 2: a block-local bulk factor for metabolism over stage-1 draws (spec D7, variant 2) |
| R3 | **mean field + learned surrogate** for the part that does not close | R1, with the non-closing reactions' contribution to the means replaced by a learned map trained on SpeedyCell runs |
| R4 | **emulated synthetic likelihood** | R1's likelihood, with the bulk factor's offset δ and covariance S learned as functions of θ from SpeedyCell runs, in place of fixing them at θ_ref |
| R5 | **neural SBI on SpeedyCell** | a neural posterior or likelihood estimator trained on simulations of the full cell; no closure at all |

### 1.1 How each route handles each row of KIM-177's inventory

The row ids are those of `dev/notes/2026-11-26-kim177-nonclosure/README.md` §4. Key:
- **PM**: product of means, the naive closure. It is not exact.
- **CA**: conservation-aware treatment of a 0/1 site (O + F₁ = 1), as spec D14 describes.
- **exact**: the route makes no approximation for that row.
- **sim**: the row is handled inside SpeedyCell simulations, which are exact draws from
  the model.
- **δ(θ)**: the row's bias in the bulk mean is absorbed by the learned offset.
- **surr**: the learned surrogate stands in for the row.
- **obs**: the row's inputs are read from observed data.
- **—**: the row does not enter that route's likelihood.

| row | what it is (count) | R1 mean field | R2 data plug-in | R3 MF + surrogate | R4 emulated SL | R5 neural SBI |
|---|---|---|---|---|---|---|
| 1a | DnaA × a 0/1 site (99) | CA, or PM as the second arm (D14) | — (replication is not in stage 1) | surr | δ(θ) | sim |
| 1b | DnaA × a multi-copy site (2) | PM | — | surr | δ(θ) | sim |
| 1c | site × site within one cascade (177) | exact (these close) | — | exact | exact | sim |
| 2 | initiator × chromosome (2) | PM | — | surr | δ(θ) | sim |
| 3 | tRNA charging (62) | PM | obs, if charged tRNA is on the panel; otherwise fixed at θ_ref | surr | δ(θ) | sim |
| 4 | cost payments (20) | PM; reads nothing downstream, so it cannot change an observable | — | PM | PM | sim |
| 5 | gene copy number × transcription (493 genes, 486 reactions) | mean dosage from the mean replication state | **latent and not handled**: fixed at one copy, as KIM-185 timed it | surr | δ(θ) for bulk; the transcript factor still takes the mean dosage | sim |
| 6 | fork steps with count-dependent constants (491) | rates rebuilt from the mean state each minute | — | as R1 | as R1 | sim |
| 7 | `max(1, pool)` floor (29 nodes, 1,434 rules) | applied to the mean | applied to the observed pool | as R1 | as R1 | sim |
| 8 | own-protein promoter, `min(c·P/180, 180)` (29 rules) | applied to the mean protein count | obs, if the protein is observed; otherwise from a stage-1 pass | as R1 | as R1 | sim |
| 9 | rounded `CellV` and the volume cap (1 + 2 nodes) | continuous volume from the lipid means, rebuilt each minute | volume from the known cell number and the panel | as R1 | as R1 | sim |
| 10 | complex AND rule, `min` over subunits (9 nodes, 18 ODE reactions) | `min` over the subunit means | stage 2 only, as R1 | as R1 | as R1 | sim |
| 11 | hook rules: counter clamp, rounding, `CellSA` | dropped or made continuous (item 3 lists which) | — | as R1 | δ(θ) | sim |
| 12 | doubling time as a threshold crossing of `CellSA` | not in the likelihood. It is the phenotype that picks the block (D9), read from SpeedyCell. | — | — | — | — |

What the routes approximate beyond the inventory:
- **R1 and R3.** The bulk factor's δ and S are fixed at θ_ref (spec, "Assumptions").
  One draw of paths per θ.
- **R2.** It ignores that metabolism feeds back on the pools between saves. The stage-1
  rates use only the observed pools. On the reduced model this shifted promoter 0607 by
  +0.65 joint SD and failed KIM-116's "close" test.
- **R4.** The emulator's error in δ(θ) and S(θ). That error grows with distance from
  the training design.
- **R5.** The amortisation error of the network, and the summary statistics it is given.
  The posterior is exact only in the limit of infinite simulations and a well-specified
  estimator.

### 1.2 Where each route's error can be measured

| route | M0 (exact sampler, SBC) | reduced model (joint posterior) | M13 (exact check with copy number and second-order reactions) |
|---|---|---|---|
| R1 | yes. Done: passes SBC, about 80× cheaper. | yes. Done: it is the joint (KIM-92). | **yes, its main test**. Row 1a under CA and PM, row 3 PM, row 5 the mean dosage (D14). |
| R2 | yes, but M0 has no pool feedback to cut | **yes. Done: fails "close" (+0.65 SD, KIM-116)** | row 5 untested: it holds copy number fixed by construction |
| R3 | no: nothing on M0 needs a surrogate | partly: the surrogate's error on the reduced model's means | yes. The surrogate is trained where M13 is exact and is scored there. M17 is its own check. |
| R4 | **yes**: δ(θ) against the exact sampler. KIM-90's postmean arm already tests δ's θ-dependence on M0. | yes: against the joint | yes |
| R5 | **yes, and needed**: SBC and the exact sampler on M0. This would be a new pilot, proposed in item 6 if R5 stays in. | yes: against the joint | yes |

### 1.3 What each route changes in the Claim, the non-goals and the positioning

| route | Claim and contribution sentence | non-goals | against PreCell (amortised SBI on the full cell, 28 group multipliers) |
|---|---|---|---|
| R1 | unchanged. "Inferred across the ODE–stochastic boundary" holds. Gate 3 (M13) decides whether the closure is good enough. | unchanged | distinct: a likelihood-based posterior on a targeted block, verified at each tier. PreCell has no exact check. |
| R2 | **changes**: it does not cross the boundary. §3.4 and the contribution sentence have to be reworded (spec risk). | unchanged | distinct, but weaker: a plug-in is not compositional inference |
| R3 | unchanged. **M17 moves ahead of M15**, because the surrogate becomes an enabler rather than §4.2's demonstration. | unchanged. §4.2's swap becomes part of the route. | distinct. The learned part is local and checked on M13. |
| R4 | unchanged in substance: the likelihood still spans the boundary, and its bulk mean is learned. "Synthetic likelihood" enters the methods. | unchanged | between the two: learned like PreCell, but local to δ and S, with a mechanistic transcript factor |
| R5 | **changes**: the mechanistic likelihood goes, and the contribution becomes amortised inference plus verification | **breaks** "no neural SBI comparison beyond the one surrogate swap" | **same ground** as PreCell. The paper becomes a comparison with it. |
