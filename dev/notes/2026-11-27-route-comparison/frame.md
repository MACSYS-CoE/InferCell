# KIM-195 frame: the cases every candidate route is costed on

Linear: [KIM-195](https://linear.app/kimpson/issue/KIM-195/compare-the-candidate-routes-for-the-full-cells-block-what-each),
item 1. This file is committed before any timing job of KIM-195 is submitted. It fixes
what is costed. It does not fix any threshold: KIM-195 sets no pass line, and the route
is KIM-187's call.

## 1. The axes

A **case** is one value on each of the four axes below. Every route in the comparison
gets a cost for every case, or a sentence saying why the case does not apply to it.

### 1.1 Block size B

B = 10, 20 and 40 free parameters. The blocks are nested: the B = 10 block is part of
the B = 20 block, and that is part of the B = 40 block. B = 20 is the headline. It is
close to the reduced model's 21 (spec D1).

### 1.2 Where the block sits

| id | placement | the parameters θ_b | how the block is ranked |
|---|---|---|---|
| **GE** | inside gene expression | log multipliers on CME reaction rates, applied **after** each minute rebuild. The rebuild cannot then overwrite them. | integrated flux over the fit window at θ = 0, as KIM-184 ranked it |
| **MET** | inside metabolism | log multipliers on ODE `kcatF` constants (158 named in `data/model/params/ode_constants.csv`), through the parameterised rate laws | integrated ODE flux over the fit window at θ = 0 |
| **BOTH** | spanning both | the top B/2 of GE, then the top B/2 of MET | as above, for each half |

The scientific block is not this ranking. KIM-180 and KIM-186 pick it by doubling time
(spec D9). The ranking here only gives each placement a block that moves states
throughout the window, so the cost is not flattered by inert parameters.

Why a multiplier on the rebuilt rate rather than on a named CME constant: the 14 named
CME constants, such as `riboKcat` and `rnaPolKcat`, are global. Each one moves every
gene at once. A per-reaction multiplier is the closer analogue of the reduced model's
per-gene promoter strengths. The cost per gradient is the same either way, because one
pass of B partials is one dual solve.

### 1.3 θ_H, the parameters outside the block

| id | θ_H | what it changes |
|---|---|---|
| **fixed** | held at the published defaults | the estimand is the block's posterior given θ_H = defaults. This is the reduced model's setting. |
| **drawn** | drawn from KIM-57's priors | the estimand is the block's posterior with θ_H marginalised, which is an honest posterior if θ_H is uncertain |

How each route handles "drawn":
- **Mean field and data plug-in.** Each likelihood evaluation needs one value of θ_H.
  "Drawn" means either a pseudo-marginal over θ_H (several solves per evaluation,
  roughly M× the cost for M draws, with a noisy gradient that NUTS cannot use) or
  putting θ_H into the sampled vector (B grows to thousands, which is the joint this
  paper rules out). So these routes are costed at **fixed** only. For "drawn" they get
  a sentence, not a number.
- **Emulator and neural SBI.** Drawing θ_H changes only the training design: the
  simulations draw θ_H from its prior, and the learned likelihood or posterior
  integrates over it. Drawn costs the same number of simulations. It may need more of
  them for the same accuracy, which is argued, not measured.

### 1.4 The data design

KIM-92's generous design (spec D2, D6 (a)), scaled to the full cell:

| | reduced model (KIM-92) | full cell (this frame) |
|---|---|---|
| cells | 200 | 200 |
| saves | every 60 s, 0 < t ≤ 3,150 s: 52 saves | the same 52 |
| fit window | 0 to `COREA_FIT_S` = 3,150 s | 0 to 3,150 s, in the full cell's own time from procid 1's initial state |
| transcript stream | all 17 genes' transcripts, exact | all **455** genes' transcripts, exact |
| bulk stream | `M0_PANEL`'s 17 pools: 1,785 rows, 884 in the fit window | the same 17 pools, where the full cell has them. A wider panel is a later rung (M10). |
| cell number | known | known |

The transcript factor's cost on this design was measured in KIM-185: 455 genes × 200
cells × 52 windows.

## 2. What a cost per posterior is built from

**CPU-h per posterior** = gradients per posterior × seconds per gradient + one-off costs
(dataset, MAP, reference ensemble, metric). **Wall-clock** is the longest chain. The
chains run in parallel, and a chain's gradients run in series.

### 2.1 Gradients per posterior

These are taken from KIM-92's run, not assumed. KIM-92 used 4 chains of 200 + 200
iterations, plus one extension of 200 per chain, from a warm diagonal metric. Its
measurements:
- leapfrog steps per 25-iteration block: 295–588 in adaptation (iterations 101–200) and
  319–375 when kept;
- 42.8–48.4 s per step;
- 451.7 CPU-h of chains and extensions over 4 chains.

That gives about 113 CPU-h per chain over 600 iterations, or about **9,200 gradient
evaluations per chain** (113 h × 3,600 / 44 s). The figure includes Turing's
per-iteration density evaluation. It sits close to the issue's 8,600.

**The projection uses 9,000 gradients per chain and 4 chains.** It is an argued figure,
scaled from KIM-92's measured counts (`dev/scripts/corea_pilot_kim92_result.md`, the
"Gradient cost" and "Cost" sections). How many steps the full cell's posterior needs is
not known. It depends on conditioning, which KIM-194's item 3 would address.

### 2.2 Seconds per gradient, by route

| term | applies to | source in this issue |
|---|---|---|
| mean-field solve, B partials | the mean-field route, and its surrogate variant | **measured**: item 3, per placement and B |
| transcript factor, **full width** (every gene carries all B partials, because its rates are read from the mean-field pools; spec D7) | the mean-field route, every placement | KIM-185's full-coordinate pass: 194.11 s per pass at Dual width 12 (job 21665221), scaled to ⌈B/chunk⌉ passes. This is **argued**, unless item 3's job times it directly at B = 10, 20 and 40. |
| transcript factor, **per gene** (pools fixed from the panel, so each gene carries only its own 2 partials) | the data plug-in's stage 1 | KIM-185: 53.61 s serial (job 21665221), 7.28 s at 16 threads (job 21665222). **Measured.** |
| block-local bulk factor for metabolism | the data plug-in's stage 2 | item 3's MET mean-field solve, or a smaller ODE-only solve. Stated in the note. |
| simulations of the full cell | emulator, neural SBI | KIM-183: 49.5 s per cycle on one core (job 21664418). **Argued** at scale. |

For the mixed route, the transcript term does not depend on the placement. A MET block
still moves the pools that every gene's transcription rate reads, through the NTPs, so
every gene's factor carries all B partials.

### 2.3 Serial and threaded

Each measured term is reported serial and threaded, at the thread count its own timing
used:
- the mean-field passes: one pass per thread, up to 4 (KIM-184);
- the transcript factor: genes over threads, up to 16 (KIM-185).

Wall-clock per chain uses the threaded figure. CPU-h uses threads × wall.

## 3. What is not in the frame

- **Accuracy.** Cost and accuracy are separate. This frame costs the routes. Where each
  route's error can be measured is item 2's table. That error is measured on M0, on the
  reduced model or in M13, not here.
- **Real data, the degradation ladder's rungs, and a second organism.** All are out of
  scope (project non-goals; M10 for the ladder).
- **The forward and predictive runs after the posterior** (M16). They cost the same for
  every route that ends in posterior draws, so they do not separate the routes.
