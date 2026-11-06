# KIM-42: the reduced model's uncertainty budget (M5). UNSCORED

**This is not a calibration result.** No truth is compared with anything here. It shows
that the method runs at this scale, and what the budget looks like on the reduced model.

**Run:** `dev/scripts/reduced_budget_kim42.sh` at `f23fd06` (main), into
`/fred/oz022/tkimpson/infercell-m5-kim42-2026-11-05/` (`budget.md`, `budget_rows.tsv`,
`jobs.tsv`, `sacct.tsv`, `SHA256SUMS`). The design is the M4 freeze's (`93abbce`), with
KIM-70's K_d diagnostic and the §12 2026-11-05 pre-registration.

**Where the numbers come from.** The budget's numbers come from `budget.md` or
`budget_rows.tsv`, unless they are marked as derived; derived ones are computed from
`w*.jls`, `plain_w*.jls` or `prediction_w*.jls`. The SE check's readings come from
`/fred/oz022/tkimpson/infercell-m5-kim42-secheck-2026-11-05/secheck.tsv` (job 20008627, at
`d3a544e`). `predict`'s hours come from `logs/`, and core-h, `TotalCPU` and MaxRSS from
`sacct.tsv`. A regeneration should be compared against `budget_rows.tsv`. Julia came from
the cluster's `julia/1.10.5` module, while the Manifest was resolved on 1.10.11.

- **Jobs:** precompile 20041038; widths 20041039 (array), whose tasks 0, 1 and 2 have the
  raw ids 20041245, 20041246 and 20041039; report 20041040. All completed with exit 0.
- **Time and cost:** `predict` took 9.39 h, 9.42 h and 12.03 h at widths ½, 1 and 2,
  against 8–9 h estimated and 11–13 h allowed. The jobs held 1,978 core-h, within the
  2,500 allowed and above the 1,000–1,600 estimated. Their `TotalCPU` is about 1,295 h,
  against the freeze's 736 CPU-h of work. The largest MaxRSS was 43.4 GiB, against 192G.
- **After the run:** all three widths ran on `znver3` (dave118, dave116, dave308), and
  `SHA256SUMS` covers 18 files. Its own SHA-256 starts `99d912fb68179769`. It replaces the
  first version (`d2245461a99ead0d`) because `sacct.tsv` was re-exported on 2026-11-06 with
  steps, `TotalCPU`, MaxRSS and `SubmitLine`; the other 17 files' hashes did not change, and
  all 18 pass `sha256sum -c`. V7 replicate 1's hash was checked against M0's pin at the start of each
  width. `report` checked that the widths share their cells' seeds. No separate
  `ran_on` audit was run; none was declared for this run.

## The headline

**At the default prior width, the held parameters account for 90.8% ± 4.1% of the
variance of μ** (I first, `per_theta_I`), and 89.8% ± 4.2% with H first. The share is
of Var(μ), the parameter uncertainty in `lactate_out` at 6,300 s. It is not a share of a
new cell's variance, which also holds the cells' own spread: of that, the held
parameters' part is 72% at width 1 (derived, from the two tables below).

**The number depends strongly on the prior width,** so it is not a property of the
model alone:

| width | mean Q ± SE | √Var(μ) ± SE | between θ_I ± SE | held share, I first | held share, H first | gap ± SE |
|---|---|---|---|---|---|---|
| ½ | 8.7986e+06 ± 1.7e+05 | 7.518e+05 ± 8.1e+04 | 2.08e+11 ± 4.8e+10 | 0.631 ± 0.102 | 0.640 ± 0.098 | -5.05e+09 ± 4.9e+09 |
| 1 | 9.5765e+06 ± 3.3e+05 | 1.473e+06 ± 2e+05 | 1.99e+11 ± 5.6e+10 | **0.908 ± 0.041** | 0.898 ± 0.042 | 2.23e+10 ± 1.3e+10 |
| 2 | 1.2835e+07 ± 8.5e+05 | 3.835e+06 ± 5.3e+05 | 2.9e+11 ± 1.4e+11 | 0.980 ± 0.009 | 0.978 ± 0.008 | 3.63e+10 ± 1e+11 |

- **Why the share moves:** θ_I's between-draw variance stays within its SEs across the
  widths, and θ_I's draws are the same at each. The held variance grows, from 3.57e+11
  to 1.97e+12 to 1.44e+13.
- **Order:** the two orders differ by about a point or less, and no gap is more than
  about 1.7 of its SEs from zero. The interaction is small on this scale.
- **The mean moves too:** Q's mean rises about 46% from width ½ to width 2 (derived).
- **The widths are not independent estimates.** They share their cells' seeds and their
  θ_H draws, scaled. Twenty columns set the held variance at all three.

**How to read the SEs.** These are the declared two-way jackknife, on 100 rows that are
about 39 effective.
- **Width 1:** 0.908 sits near `secheck`'s held-dominant shape (shares 0.92–0.93, ratios
  of mean SE to empirical SD 1.01–1.07), so its SE is about right on those made-up
  matrices.
- **Width 2:** 0.980 is near the ceiling, between `secheck`'s cases at 0.97 (ratios 1.04)
  and those at 0.994, where it failed narrowly (1.51 and 1.56). So its held-share SE may
  be too wide, by up to about 50% (§12 2026-11-05). Its between-θ_I SE is not to be read
  as a normal interval there: in those failing cases, ±2 SE held that estimate only about
  80% of the time.

**The sentence the plan aimed for is not claimed as written.** "x% of this prediction's
uncertainty comes from parameters no measurement supports" overstates this run in two
ways. θ_H are 15 promoter strengths marked `:asserted`
(`src/organisms/coreA/transcription.jl`): their medians rest on proteomics through a proxy,
protein copy number over 180, but no measurement sets their width. And θ_I's posterior
comes from V7 replicate 1, a fit to M0's synthetic data, so the 9% is not what real
measurements would leave. What the run supports:

*At the default prior width, about 91% ± 4% of Var(μ), the parameter uncertainty in
`lactate_out` at 6,300 s, comes from the 15 promoter strengths that M0's fit cannot reach
because their genes are cut from M0. Their medians are a proteomics proxy, their spread
(geometric SD 2) is this project's choice, and the share runs from 63% at half that spread
to 98% at double.*

Every other constant is held fixed or derived, so this is a split among these 20 drawn
parameters only. The widths scale the log-scale σ at a fixed median, which gives geometric
SDs of 1.41, 2 and 4.

**What the 9% rests on.** θ_I's share comes from one replicate of V7's 50, fitted at 600 s
and carried to 6,300 s. The range from 63% to 98% varies only θ_H's width, not θ_I's
posterior. And the transfer assumes θ_I's posterior would not move if θ_H's counterparts in
M0 did. The predictive note (§3.2) states that this is untested and that the paper says so.

## The filter

| width | acceptance | ESS, `per_theta_I` | ESS, `joint` | largest abs corr(log θ_I, log θ_H) | clean | clipped | starved |
|---|---|---|---|---|---|---|---|
| ½ | 1.0000 | 2000 | 2000 | 0.000 | 2402 | 5598 | 0 |
| 1 | 1.0000 | 2000 | 2000 | 0.000 | 2101 | 5871 | 28 |
| 2 | 0.9355 | 1869 | 1871 | 0.046 | 788 | 6508 | 704 |

The cell columns count μ̂'s 8,000 baseline cells at each width, before V. No row failed,
and no cell failed or fell below a floor.
- **V acted only at width 2,** and removed 6.45% of pairs. At ½ and 1 it kept everything,
  so the two weight sets are identical there and the correlations are zero by
  construction. Even at width 2 it induces almost none.
- **Clipping is the common state here.** 70%, 74% and 90% of those cells clipped at least
  once (derived, clipped and starved cells together, since a starved cell has clipped),
  which V's frozen rule allows. Starved cells are 0%, 0.35% and 8.8% (derived).
- **Why so much more than KIM-35, probably: θ_I, not the drawing of θ_H.** KIM-35's arm H
  drew the same 15 promoters from the same prior and clipped 5.5% of 200 cells at 6,300 s,
  with θ_I at V6's fit. The evidence, all derived and none of it in the declared report:
  - Width ½ still clips 70%, so narrowing θ_H does not remove it.
  - krnadeg is 9.66 at V6's point and 3.59–5.70 on these 100 rows.
  - Across those rows' krnadeg deciles at width 1, the clipped fraction falls at every
    step, from 84% to 61%.
  - In mean-field solves (a review check on the login node at `f23fd06`, with no job id),
    minimum GTP is 7,280 particles at V6's point and 4,404 at replicate 1's. Swapping
    krnadeg alone to V6's value gives 7,292. Swaps of S_0607, S_0779 and kcatF_R_ENO
    leave it within about 30 of 4,404; kcatF_R_FBA's swap did not finish.

  That swap was not run on stochastic cells, so the link is unproven at the cell level.
  Replicate 1 is not an outlier: its krnadeg ranks 32nd of V7's 50, and only 6 of the 50
  reach V6's value. **So the clipping, K_d and Var(μ) figures here hold for replicate 1's
  θ_I, and another replicate could give different ones.**

## The cells' own spread

| width | Var(μ) | E[σ²] | Var(μ) / (Var(μ) + E[σ²]) |
|---|---|---|---|
| ½ | 5.652e+11 | 4.446e+11 | 0.5597 |
| 1 | 2.169e+12 | 5.671e+11 | 0.7928 |
| 2 | 1.47e+13 | 1.741e+12 | 0.8941 |

σ_b is not used: it enters only for a future noisy observation, not for Q.

## ΔQ under the knockdown (secondary): its split is not resolved

**At 4 cells per pair, this run cannot resolve the parameter variance of ΔQ, so it has no
variance split for ΔQ.** Cell noise swamps it. The estimate of Var(μ) came out negative
and is reported as 0 at widths 1 and 2, and at width ½ it is about 0.25 of its SE from
zero. The components are clipped one by one and so no longer sum, which is why the report
prints `NaN` and `Inf` for the shares at widths 1 and 2 and `NaN` for the SEs of all
three. That is arithmetic, not a failed computation, and those entries are not numbers
to quote.

ΔQ's mean is defined, on the kept pairs at dose 0.3:

| width | mean ΔQ ± SE |
|---|---|
| ½ | -1.1831e+05 ± 9.4e+03 |
| 1 | -98421 ± 1.4e+04 |
| 2 | -39525 ± 3.2e+04 |

The interval of ±2 SE reaches zero at width 2. No intervened cell failed at any width. M0's
ΔQ has the same kind of limit, stated in prose beside its numbers (§12 2026-11-05).

## K_d diagnostic, for KIM-47

Not scored, and it gates nothing. The fraction of cells whose `minima[:M_gtp_c]` fell
below 2,018 particles (a lower bound) and 2,160 (an upper bound):

| width | baseline, all 8 of every pair | baseline, kept pairs | knockdown, kept pairs |
|---|---|---|---|
| ½ | 0.9694–0.9776 (n = 16,000) | 0.9694–0.9776 (16,000) | 0.9624–0.9706 (8,000) |
| 1 | 0.9467–0.9576 (16,000) | 0.9467–0.9576 (16,000) | 0.9391–0.9500 (8,000) |
| 2 | 0.9325–0.9441 (16,000) | 0.9278–0.9403 (14,968) | 0.9176–0.9340 (7,484) |

**So nearly every cell reaches the threshold at least once; the run does not say for how
long.** A cell counts here if its GTP was *ever* below the threshold. The freeze ties the
one-way reduction's basis failure to fractions of rebuild instants: the dose screen, with
θ_I drawn from its prior rather than this posterior, found GTP below K_d for 21–55% of
them in 12 of 20 draws, and for none in 6 (§12 2026-10-30). This run stores no such
fraction, so it does not show that the failure reaches these cells.

**Many go well below the threshold, usually briefly** (derived from μ̂'s 8,000 baseline
cells at each width):
- 68%, 69% and 72% reach exactly 0 GTP, a real state of the model. It is the same event as
  a GTP clip, in all but 3 of the 24,000 cells.
- Among the cells that clip GTP, the median time clipped is 10, 27 and 274 of the run's
  6,300 drains.
- At width 2, 12.2% of them, the 704 starved cells, are clipped for over half the run.

KIM-47's table should carry these numbers.

## Attribution (descriptive)

First-order correlations with the line means, at width 1. Twenty columns cannot separate
the held parameters.
- **θ_I:** krnadeg r = -0.95, S_JCVISYN3A_0607 r = -0.37, S_JCVISYN3A_0779 r = 0.23,
  kcatF_R_ENO r = -0.10, kcatF_R_FBA r = 0.03.
- **θ_H, largest 5:** S_JCVISYN3A_0221 r = 0.67, S_JCVISYN3A_0220 r = 0.55,
  S_JCVISYN3A_0651 r = 0.50, S_JCVISYN3A_0131 r = 0.49, S_JCVISYN3A_0606 r = 0.48.

## What this does not show

- It is not calibration, on the reduced model or anywhere else.
- It says nothing about how the reduced model compares with real measurements.
- Its code had no adversarial review before the run: the user waived `/check-PR` for
  #151. What tested it was the SE check and the two smoke passes.
