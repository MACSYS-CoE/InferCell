# A True Three-Module Minimal Cell for the ISAB Talk

**Date:** 2026-08-06
**Talk:** MACSYS ISAB Meeting, 12 August 2026 (`dev/talks/ISAB/talk.tex`)
**Branch:** `minimal-cell-example`

## Problem

The deck motivates itself with a cell-level number — slide 3 pushes a parameter
range through a forward model, slide 8 asks *"Is doubling time 110 ± 5 min, or
60–200?"* — and then reports two parameter corner plots for a rig that is not a
cell. Two of its three modules are the same four reactions in two formalisms
(`TranscriptionTranslation` and `StochasticGeneExpression` share all four
parameters and register under the same names). The deck's own source concedes
this at `talk.tex:1071`, where the slide was retitled from "test cell" to "test
harness" to head off the obvious question.

Two consequences:

1. **The story is muddled for a general audience.** The setup slide promises a
   mini whole-cell model; the results slide delivers an abstract methods check.
   The transition between them is the part that reads worst when practised.
2. **The headline figure overstates what was measured.** The
   conditioned-vs-unconditioned contrast in `step3_ssa_posteriors.png` is
   dramatic *because* the simulation block shares all four of its parameters
   with the differentiable block. Presented as a boundary result between
   distinct modules, it claims a transfer effect that the experiment did not
   test. "How many parameters do the two blocks share?" undoes the slide.

## Decision

Switch the deck to the **bursty-regulator cell**, which already exists in the
repository and is exercised by `test/run_integ_iterative.jl`. It is a genuine
three-module cell: three distinct subsystems, not one subsystem described twice.

| Module | Formalism | Inference | Biology |
|---|---|---|---|
| `LightMetabolism` | ODE | NUTS | ATP / NTP / AA; production gated by enzyme via Michaelis–Menten |
| `TranscriptionTranslation` | ODE | NUTS | expression of the enzyme gene; its protein feeds metabolism |
| `BurstyGeneExpression` | SSA (jump) | ABC-SMC | telegraph promoter, discrete counts, low copy number |

Modules 1 and 2 are **state-coupled in both directions** — the enzyme raises ATP
production, and mRNA drives the TX/TL sinks — which closes the autocatalytic
loop. Module 3 is **parameter-coupled** to module 2 through `k_tl`,
`gamma_mRNA` and `gamma_protein`: the regulator is translated by the same
ribosomes and cleared by the same degradation machinery, so those are literally
the same three numbers.

**11 unique parameters. 3 shared.** (8 in the differentiable block after
deduplication, plus `k_on`, `k_off`, `k_tx_burst`.)

Rejected alternatives: building a growth/division module (best story — its
observable would be the doubling-time distribution the deck already asks about —
but new modelling, new summary statistics and fresh synthetic data six days out,
with no fallback); and relabelling the existing simulation module as a distinct
subsystem without rerunning anything (cheapest, but see consequence 2 above).
Growth/division remains the right move after the talk.

## Inference flow

Unchanged from `sequential_infer` — no new machinery:

1. Joint NUTS over the differentiable block (both ODE modules, all 8 parameters,
   conditioned on the bulk metabolic/expression time series).
2. `boundary_condition(chain, [bge]; method=:kde)` fits KDE priors for the three
   shared parameters and leaves `k_on`, `k_off`, `k_tx_burst` at their own
   priors.
3. ABC-SMC over the simulation block's 6 parameters under those priors.

Plus an **unconditioned baseline**: ABC-SMC on the same simulation block, same
data, with all six parameters at their original priors. The contrast between the
two is the result.

## The claim the figure supports

Three of the six simulation-block parameters should sharpen onto truth; three
should stay as broad as their priors — because nothing in the bulk metabolic
data knows about promoter switching. The figure shows information crossing a
module boundary *and stopping exactly where the coupling stops*.

This is a stronger claim than the current figure's, and a safer one: the
mechanism is visible in the picture rather than asserted in the caption.

## Deliverables

1. `test/run_talk_bursty_boundary.jl` + `.slurm` — runs the conditioned and
   unconditioned inferences, exports particles, weights, the conditioned mask
   and ground truth as CSV. Follows the existing `test/run_integ_*.jl` pattern.
2. `dev/talks/ISAB/figures/plot_bursty_boundary.py` — corner plot, blue
   unconditioned / green conditioned / orange truth, matching the existing
   figure style (`corner` 2.2.3).
3. Slide edits in `dev/talks/ISAB/talk.tex`:
   - **Setup slide** (`A Three-Module Test Harness` → a cell, not a harness):
     three distinct modules; drop the "same reactions, two formalisms"
     annotation and the θ-shared caveat. Replace the closing "Goal: the joint
     posterior" bullet with the question the results answer — *"The bulk data
     never sees the regulator gene. Can it still tell us anything about
     single-molecule bursting?"* This is the fix for the awkward transition:
     the results slide currently answers a question the setup slide never asked.
   - **Results slide**: one figure, the 6 simulation-block parameters,
     unconditioned vs conditioned. The 8-parameter NUTS corner moves to backup.
   - **Backup**: the NUTS corner (identifiability), retitled so it still pays
     off the identifiability bullet on the "What Becomes Possible" slide.
   - Update the stale cautions at `talk.tex:1071` and `talk.tex:1231`, which
     currently warn against exactly the swap being made.

## Risks

- **The effect may be visually weaker.** With 3 of 6 shared rather than 4 of 4,
  the collapse will be less dramatic than `step3_ssa_posteriors.png`. Mitigation:
  the run happens first, before any slide is rewritten. If the three shared
  parameters do not visibly sharpen, we fall back to the current figure with the
  honest "same reactions, two formalisms" framing and lose only the slide work.
  `test/run_integ_iterative.jl` asserts tightening on these three parameters, so
  there is prior evidence the effect exists.
- **ABC in 6 dimensions with 200 particles is thin.** Running at 400 particles
  for figure quality; walltime set accordingly.
- **Wall-clock.** Six days to the talk. The compute is the long pole and is
  queued before the slide work starts.

## Run log

**2026-08-06, smoke run (job 15152689)** — 30 particles, 2 populations, 10
replicates, 100 NUTS samples. Script ran end to end and exported cleanly.
Boundary behaved as designed: `Conditioned: Bool[0,0,0,1,1,1]`, 3 of 6.

| param | shared | truth | uncond mean±sd | cond mean±sd | sd ratio |
|---|---|---|---|---|---|
| k_on | no | 0.2 | 0.854 ± 1.605 | 1.577 ± 1.469 | 1.09 |
| k_off | no | 0.5 | 0.5 ± 0.495 | 1.324 ± 2.243 | 0.22 |
| k_tx_burst | no | 10.0 | 2.546 ± 0.468 | 3.217 ± 1.173 | 0.40 |
| k_tl | **yes** | 2.0 | 0.649 ± 0.808 | 2.151 ± 0.004 | 210 |
| gamma_mRNA | **yes** | 0.5 | 0.907 ± 0.377 | 0.531 ± 0.003 | 127 |
| gamma_protein | **yes** | 0.1 | 1.16 ± 0.484 | 0.098 ± 0.001 | 875 |

Two findings, both from inspecting the smoke figures rather than the table:

1. **Axis ranges must use weighted quantiles.** ABC-SMC leaves near-zero-weight
   survivors far out in the tails; ranging on unweighted quantiles let those set
   the limits and squashed every posterior into the corner of its panel. Fixed
   in `plot_bursty_boundary.py` (`weighted_quantile`).
2. **The smoke result was flattered by a non-converged chain.** At 100 NUTS
   samples the differentiable posterior spans 0.960–0.976 for `k_tx` with truth
   at 1.0 — far too narrow, and it *excludes* truth. The boundary faithfully
   passed that overconfidence on, which is why the conditioned sd is 0.001–0.004
   and the conditioned mean for `k_tl` sits at 2.151 rather than 2.0. Do not
   quote smoke numbers. Production uses 1000 samples.

**2026-08-06, production run (job 15152954)** — 400 particles, 6 populations,
50 replicates, 1000 NUTS. 525 s conditioned + 529 s unconditioned; 18 min wall.
Boundary again `Bool[0,0,0,1,1,1]`.

| param | shared | truth | unconditioned | conditioned | sd ratio |
|---|---|---|---|---|---|
| k_on | no | 0.2 | 2.106 ± 1.979 | 2.494 ± 2.031 | 0.97 |
| k_off | no | 0.5 | 0.998 ± 1.080 | 0.893 ± 0.883 | 1.2 |
| k_tx_burst | no | 10.0 | 5.109 ± 4.445 | 4.259 ± 1.447 | 3.1 |
| k_tl | **yes** | 2.0 | 4.502 ± 3.430 | **2.074 ± 0.086** | 40 |
| gamma_mRNA | **yes** | 0.5 | 0.437 ± 0.403 | **0.514 ± 0.031** | 13 |
| gamma_protein | **yes** | 0.1 | 0.384 ± 0.328 | **0.098 ± 0.002** | 187 |

The headline claim holds and is now measured: the three shared parameters
collapse onto truth (13–187× tighter), the other three barely move. The
differentiable chain converged this time — truth sits inside every marginal and
the corner shows clean tilted ellipses, so the backup identifiability slide
earns its place.

**The promoter rates are not recovered, and that is not a bug.** `k_on` (truth
0.2, posterior ≈ 2.1) and `k_tx_burst` (truth 10, posterior ≈ 4.3) miss badly in
*both* runs. The reason is identifiability, not sampler failure: the summary
statistics are replicate-**mean** trajectories, and a telegraph promoter enters
the mean only through the effective burst rate
`k_tx_burst · k_on/(k_on + k_off)`. That combination is identifiable and is
recovered — truth 2.857, unconditioned 3.307 ± 2.986, conditioned
2.793 ± 0.244, so the boundary sharpens it 12× as well. This also explains why
`k_tx_burst`'s own marginal tightened 3.1× despite not being conditioned.

Consequences taken: the results slide says the promoter rates "barely move ---
the bulk data never sees them, and mean trajectories cannot separate them
anyway" rather than implying they are well constrained, and the numbers plus the
full explanation sit in a comment above that slide so the answer is to hand in
Q&A. Recovering individual promoter rates would need distributional summaries
(Fano factor, autocorrelation) instead of means — a real next step for the
summary-statistic layer, and arguably a better demo of why the SSA module exists
at all.

## Out of scope

Growth/division module; iterative (Level-2) boundary protocol, which stays in
backup; any change to `src/`. This is a talk deliverable built from existing
framework code.
