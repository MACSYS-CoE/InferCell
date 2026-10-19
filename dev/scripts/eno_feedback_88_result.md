# Issue #88 A3: the cost of reducing the reverse feedback

**Run:** arrays 18825531–18825546 (400 tasks) and merge 18825547, all at `8813dca`,
submitted from a detached worktree at that pushed commit. The `src` hook they use,
`override_rebuild_pools!`, passed the full suite at `6da2398` (job 18821310, 28361
of 28361), and `src` and `test` are unchanged from there to `8813dca`.
**Script:** `dev/scripts/eno_feedback_88.jl`.
**Re-merged** after the #92 review, by job 18836376 at `2b039a0`. The merge now also checks
that every input file comes from one commit (`8813dca`, confirmed) and one kind of run.
Every number is identical; only the header and a caption line changed. The original
merge's output is kept as `merged_8813dca.md`.
**Cost:** 86.97 CPU-h: 17.4 for the ensembles, 69.6 for the path replays, 0.015
for the merge.
**Raw outputs:** `/fred/oz022/tkimpson/infercell-88-a3-archive-2026-10-19/`: one
directory per stage, arm and truth (400 files, `MD5SUMS`), plus `merged.md`,
`merged.jls` and the Slurm logs.
**Regenerate:** `dev/scripts/eno_feedback_88.sh <outdir outside the worktrees>`.

## Setup

Reduced syn3A model (`d11_models`, published drains). Four arms:

| Arm | Transcription and translation rate constants are rebuilt from |
|---|---|
| `live60` | the live pools every 60 s (the published model) |
| `nominal` | the initial ATP, GTP and charged tRNA, every 60 s |
| `popmean` | the 200-cell mean of those pools at the truth, from `live60`'s ensemble, at each 60 s rebuild. A fixed trajectory: mean field's reverse edge in the ensembles, which run at the truth, but frozen pools in the path replays at ln θ ± δ, where a mean-field loop would move the mean with θ |
| `live6` | the live pools every 6 s (spec R11). Checked separately on the login node at `a934e98`'s `src`: 20 refreshes in 120 s |

Two truths:
- **Truth 1:** 15.7's, `Xoshiro(1507)`, with ptsG 4.36, ENO 57.5, FBA 32.1.
- **Truth 2:** `Xoshiro(8802)`, with ptsG 0.80, ENO 85.5, FBA 73.1.
- The merge checked that every arm and stage used the same truth.

**`live60` at truth 1 reproduces 16a7 exactly** at δ = 0.2, so the script is a faithful
extension of it:
- observed-information ratio, ENO 0.996 (0.978 at mean + 2 SE) and FBA 0.842 (0.722);
- expected-information ratio, ENO 0.999 and FBA 0.995.

Two scores:
1. **Path information.** As `eno_path_information_16a7.jl`: 40 cells (seeds
   30001–30040), each path replayed at ln θ ± δ.
   - Parameters: ENO, FBA and D11's three promoters, at δ ∈ {0.05, 0.2}.
   - The ratio is the complete-data Fisher information of 200 cells, as posterior
     over prior SD.
   - It bounds what the transcripts can carry.
2. **The bulk metabolite likelihood.**
   - 200-cell forward ensembles per arm and truth (seeds 40001–40200).
   - The score is Σ_t (log x̄_arm − log x̄_live60)²/(2σ_b²) at σ_b = 10%, the KL
     between the two arms' bulk likelihoods, in nats.
   - The noise column beside it is the same sum's expected value with no shift,
     from `live60`'s split halves. It is conservative, because the arms share cell
     seeds with `live60`.

## Verdict 1: the 60 s cadence does not explain the ENO null

Expected-information ratio, ENO and FBA from their own 2 × 2 block (the 16a7
quantity), δ = 0.05:

| Truth | `live60` ENO, FBA | `live6` ENO, FBA |
|---|---|---|
| 1 | 0.999, 0.995 | 0.999, 0.994 |
| 2 | 0.896, 0.754 | 0.923, 0.904 |

- **Rebuilding ten times as often adds no path information** at either truth. The
  bootstrap intervals overlap; δ = 0.2 agrees, and the full table is below.
- **It leaves the bulk metabolites unchanged too.** `live6` against `live60` totals
  81.9 nats against a noise of 115 at truth 1, and 42.8 against 55.7 at truth 2.
- **So spec §1's "likely reason", the 60 s rebuild, is refuted at these two truths.**
  What remains is the weak reverse gain. The discussion note puts that down to
  saturated rate laws. Saturation itself was not varied here.
- **The ENO null is truth-dependent.**
  - Truth 1 reproduces 16a7's null: ENO 0.999 (16a7: 0.999 at the same truth).
  - **At truth 2, ENO may carry a little path information.**
    - The outer-product (expected-information) estimator gives 0.90–0.93,
      consistent across δ.
    - The curvature (observed-information) estimator does not confirm it: the
      mean curvature over 40 cells is negative, which reports as NaN.
  - **FBA at truth 2 is undetermined.** It reads 0.75 at δ = 0.05 and 0.92 at
    δ = 0.2, with a bootstrap of 0.68–0.99.
  - **Desaturation is a hypothesis, not a finding.** Truth 2 has a low ptsG (0.80
    against 4.36), so glucose uptake and the pools may be lower and the rate laws
    less saturated. `rate_constant_elasticity` at truth 2 would test it.

## Verdict 2: the one-way reduction is not uniformly acceptable

Bulk metabolite likelihood, total over the 17-pool panel, in nats at σ_b = 10%. The
statistic scales as 1/σ_b², so divide by 4 for σ_b = 20%:

| Truth | noise | `nominal` | `popmean` | `live6` |
|---|---|---|---|---|
| 1 | 115 | 152 | 74.8 | 81.9 |
| 2 | 55.7 | **707** | **198** | 42.8 |

- **`popmean` (population-mean pools, note §7's proposal):**
  - **Truth 1:** indistinguishable from the full model. It is below the noise in
    the total and in every pool but external lactate (0.69 against 0.60, marginal).
  - **Truth 2:** about 140 nats above the noise at σ_b = 10%, or about 35 at 20%.
  - The shift is concentrated in G6P and F6P (31 and 33 against noise 1.6 each),
    charged tRNA (22.5 against 3.1), lactate (24.5 against 5.0) and Pi (26.9
    against 9.2).
  - So dropping the cell-to-cell part of the feedback is not second order
    everywhere. It moves the bulk likelihood by many nats where the reverse edge
    is live.
- **`nominal` (fixed reference pools):**
  - About 37 nats above the noise at truth 1: G6P, F6P, external lactate and FDP
    shift.
  - About 650 at truth 2.
  - Cutting the loop at fixed pools is not acceptable for the bulk metabolites.
- **Path information for the promoters barely moves under frozen pools**
  (expected-information ratio, δ = 0.05). For `popmean` these are frozen-pool
  numbers, not mean field's: a mean-field path would add a promoter → mean pool →
  rate constant channel, which neither reduced arm has.

  | Truth | Arm | GAPD | PGI | ptsG |
  |---|---|---|---|---|
  | 1 | `live60` | 0.027 | 0.059 | 0.030 |
  | 1 | `nominal` | 0.025 | 0.056 | 0.032 |
  | 1 | `popmean` | 0.028 | 0.062 | 0.024 |
  | 2 | `live60` | 0.018 | 0.037 | 0.075 (0.063–0.101) |
  | 2 | `nominal` | 0.017 | 0.035 | 0.049 (0.042–0.066) |
  | 2 | `popmean` | 0.017 | 0.034 | 0.061 (0.050–0.085) |

  - **Truth 1:** GAPD and PGI agree within their bootstrap intervals.
    - ptsG under `popmean` (0.024) falls just below `live60`'s interval
      (0.026–0.040), though `live60`'s 0.030 lies inside `popmean`'s (0.019–0.034).
  - **Truth 2:** the same holds for ptsG.
    - `popmean`'s 0.061 is just below `live60`'s interval (0.063–0.101), and
      `live60`'s 0.075 is inside `popmean`'s (0.050–0.085).
    - `nominal`'s 0.049 (0.042–0.066) moves further.
- **What frozen pools lose outright:** under `nominal` or `popmean`, ENO and FBA no
  longer reach the path, so their ratio is exactly 1.000.
  - This holds by construction and is the expected check.
  - At truth 2 it discards whatever ENO information the path had (0.90 by the
    outer-product estimator, not confirmed by the curvature).
  - It is not what mean field would lose. A mean-field loop keeps ENO reaching the
    rebuild through the θ-dependent mean pools, and this run does not measure that.
  - The metabolites identify ENO anyway (A2: 0.13 at 10%).

**Implication for track B (#89):**
- **Verdict 2 rests on the bulk score**, not on the uncertain truth-2 path
  information.
- **In the ensembles, run at the truth, `popmean` is the mean-field description of
  the reverse edge.** A mean-field gene-expression module can only feed back
  population-mean pools.
- **It costs nothing at truth 1 and about 140 nats at truth 2.** Note §14 step 4 asks
  this question, together with R and the Jensen gap (B3): is mean field the wrong
  second description? This says it is, at least at a low-ptsG truth.
- **The options:**
  - keep per-cell feedback, so path sampling, where the reverse gain is strong;
  - correct the population-mean feedback to second order;
  - limit the approximate method's claim to the part of the prior where the
    reduction is measured to be cheap.
- V7's datasets span ptsG's prior, which makes them the place to find where
  `popmean` breaks. Two truths cannot say what fraction of the prior that is.

## Caveats

- Two truths, 40 path cells and 200 ensemble cells per arm.
- The bulk statistic is distributional: arms are not paired path by path, because
  rebuilds change how the jump block consumes random numbers.
- `popmean`'s mean pools come from a 200-cell ensemble, so they carry that
  ensemble's sampling noise.
- **The noise column's split halves and the bootstrap draws follow the input files'
  sorted order.** Regenerating with other than 10 ensemble tasks per arm gives
  other halves, so the noise column changes within its own noise.
- **An observed-information cell is NaN where the mean curvature over 40 cells is
  negative,** which is noise.
  - Its bracket, for example "NaN (0.271)" for ENO at truth 2, assumes the
    curvature is 2 SE more positive than measured. It is not evidence of
    information.
  - It is reported but not used.

## The full merge output

### Truth 1 (`Xoshiro(1507)`)

S_JCVISYN3A_0607 2.973, S_JCVISYN3A_0445 0.5672, S_JCVISYN3A_0779 4.364, krnadeg 8.775, kcatF_R_ENO 57.46, kcatF_R_FBA 32.11, kcatR_R_FBA 0.3012, kcatR_R_ENO 21.67 — the same in every arm and stage (200 files).

#### Path information (200 cells)

| Arm | δ | cells | Expected-information ratio: ENO | FBA | S_0607 (GAPD) | S_0445 (PGI) | S_0779 (ptsG) | ENO, FBA alone | Observed-information ratio: ENO | FBA | S_0607 (GAPD) | S_0445 (PGI) | S_0779 (ptsG) |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| live60 | 0.05 | 40 | 0.999 (0.999–1.000) | 0.996 (0.993–0.998) | 0.027 (0.023–0.034) | 0.059 (0.052–0.075) | 0.030 (0.026–0.040) | 0.999, 0.995 | 0.987 (0.839) | 0.608 (0.373) | 0.023 (0.023) | 0.059 (0.059) | 0.028 (0.028) |
| live60 | 0.2 | 40 | 0.999 (0.999–1.000) | 0.996 (0.994–0.998) | 0.027 (0.023–0.034) | 0.059 (0.052–0.076) | 0.030 (0.025–0.040) | 0.999, 0.995 | 0.996 (0.978) | 0.842 (0.722) | 0.023 (0.023) | 0.059 (0.059) | 0.028 (0.028) |
| nominal | 0.05 | 40 | 1.000 (1.000–1.000) | 1.000 (1.000–1.000) | 0.025 (0.022–0.032) | 0.056 (0.048–0.071) | 0.032 (0.025–0.045) | 1.000, 1.000 | 1.000 (1.000) | 1.000 (1.000) | 0.022 (0.022) | 0.058 (0.058) | 0.028 (0.028) |
| nominal | 0.2 | 40 | 1.000 (1.000–1.000) | 1.000 (1.000–1.000) | 0.025 (0.022–0.032) | 0.055 (0.048–0.071) | 0.032 (0.026–0.044) | 1.000, 1.000 | 1.000 (1.000) | 1.000 (1.000) | 0.022 (0.022) | 0.057 (0.057) | 0.027 (0.027) |
| popmean | 0.05 | 40 | 1.000 (1.000–1.000) | 1.000 (1.000–1.000) | 0.028 (0.024–0.039) | 0.062 (0.055–0.078) | 0.024 (0.019–0.034) | 1.000, 1.000 | 1.000 (1.000) | 1.000 (1.000) | 0.023 (0.023) | 0.059 (0.059) | 0.028 (0.028) |
| popmean | 0.2 | 40 | 1.000 (1.000–1.000) | 1.000 (1.000–1.000) | 0.028 (0.023–0.039) | 0.062 (0.054–0.078) | 0.024 (0.019–0.033) | 1.000, 1.000 | 1.000 (1.000) | 1.000 (1.000) | 0.023 (0.023) | 0.059 (0.059) | 0.028 (0.028) |
| live6 | 0.05 | 40 | 0.999 (0.999–0.999) | 0.995 (0.994–0.997) | 0.029 (0.025–0.035) | 0.065 (0.054–0.085) | 0.029 (0.024–0.039) | 0.999, 0.994 | NaN (0.876) | 0.960 (0.406) | 0.023 (0.023) | 0.059 (0.059) | 0.028 (0.028) |
| live6 | 0.2 | 40 | 0.999 (0.999–1.000) | 0.995 (0.994–0.997) | 0.029 (0.025–0.036) | 0.065 (0.054–0.084) | 0.029 (0.024–0.039) | 0.999, 0.994 | NaN (0.982) | 0.859 (0.720) | 0.023 (0.023) | 0.059 (0.059) | 0.028 (0.028) |

Expected-information cells are joint over the five parameters, with a 95% bootstrap over cells. "ENO, FBA alone" is the same from their 2 × 2 block only, the 16a7 script's quantity. Observed-information cells are each parameter's own curvature alone, at the mean and, in brackets, at mean + 2 SE (the most information the cells allow).

#### Bulk metabolite likelihood (200-cell means, 60 s saves)

Per pool: Σ_t (log x̄_arm − log x̄_live60)² / (2σ_b²) at σ_b = 10% (it scales as 1/σ_b²), with `live60`'s split-half noise on the same scale in the first column. The arms share cell seeds with `live60`, so until their paths part the difference is below independent noise, and the noise column is conservative. Above the noise is a real shift; above ~1 nat the reduction moves the bulk likelihood materially.

| Pool | noise | nominal | popmean | live6 |
|---|---|---|---|---|
| `M_g6p_c` | 5.03 | 12.9 | 1.98 | 2.87 |
| `M_f6p_c` | 4.63 | 12.8 | 1.76 | 2.77 |
| `M_lac__L_c` | 14.5 | 16.5 | 8.84 | 6.93 |
| `M_g3p_c` | 2.45 | 1.42 | 0.964 | 0.793 |
| `M_dhap_c` | 0.788 | 0.374 | 0.306 | 0.243 |
| `M_trna_chg_c` | 0.62 | 0.741 | 0.477 | 0.469 |
| `M_pi_c` | 25.4 | 29.2 | 18.4 | 18.5 |
| `M_3pg_c` | 30.4 | 37.9 | 20.8 | 23.1 |
| `M_2pg_c` | 26.8 | 31.1 | 18.4 | 21.2 |
| `M_gtp_c` | 2.32 | 2.31 | 1.17 | 2.88 |
| `M_fdp_c` | 0.0073 | 0.258 | 0.00557 | 0.00497 |
| `M_atp_c` | 0.233 | 0.559 | 0.135 | 0.346 |
| `M_lac__L_e` | 0.604 | 4.71 | 0.693 | 0.523 |
| `M_amp_c` | 0.833 | 0.855 | 0.499 | 0.656 |
| `M_gmp_c` | 0.355 | 0.388 | 0.174 | 0.378 |
| `M_gdp_c` | 0.225 | 0.189 | 0.14 | 0.152 |
| `M_adp_c` | 0.0556 | 0.0311 | 0.0274 | 0.0292 |
| **total** | **115** | **152** | **74.8** | **81.9** |

### Truth 2 (`Xoshiro(8802)`)

S_JCVISYN3A_0607 5.455, S_JCVISYN3A_0445 1.931, S_JCVISYN3A_0779 0.7992, krnadeg 2.698, kcatF_R_ENO 85.53, kcatF_R_FBA 73.14, kcatR_R_FBA 0.6861, kcatR_R_ENO 32.26 — the same in every arm and stage (200 files).

#### Path information (200 cells)

| Arm | δ | cells | Expected-information ratio: ENO | FBA | S_0607 (GAPD) | S_0445 (PGI) | S_0779 (ptsG) | ENO, FBA alone | Observed-information ratio: ENO | FBA | S_0607 (GAPD) | S_0445 (PGI) | S_0779 (ptsG) |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| live60 | 0.05 | 40 | 0.913 (0.879–0.965) | 0.781 (0.677–0.989) | 0.018 (0.015–0.023) | 0.037 (0.032–0.048) | 0.075 (0.063–0.101) | 0.896, 0.754 | NaN (0.271) | NaN (0.085) | 0.018 (0.017) | 0.033 (0.033) | 0.068 (0.067) |
| live60 | 0.2 | 40 | 0.932 (0.893–0.981) | 0.929 (0.885–0.990) | 0.017 (0.015–0.022) | 0.036 (0.031–0.047) | 0.073 (0.061–0.099) | 0.920, 0.924 | NaN (0.697) | NaN (0.578) | 0.018 (0.017) | 0.033 (0.033) | 0.068 (0.067) |
| nominal | 0.05 | 40 | 1.000 (1.000–1.000) | 1.000 (1.000–1.000) | 0.017 (0.014–0.024) | 0.035 (0.031–0.045) | 0.049 (0.042–0.066) | 1.000, 1.000 | 1.000 (1.000) | 1.000 (1.000) | 0.017 (0.017) | 0.031 (0.031) | 0.064 (0.064) |
| nominal | 0.2 | 40 | 1.000 (1.000–1.000) | 1.000 (1.000–1.000) | 0.017 (0.014–0.023) | 0.035 (0.031–0.045) | 0.049 (0.041–0.067) | 1.000, 1.000 | 1.000 (1.000) | 1.000 (1.000) | 0.017 (0.017) | 0.031 (0.031) | 0.064 (0.064) |
| popmean | 0.05 | 40 | 1.000 (1.000–1.000) | 1.000 (1.000–1.000) | 0.017 (0.014–0.023) | 0.034 (0.029–0.045) | 0.061 (0.050–0.085) | 1.000, 1.000 | 1.000 (1.000) | 1.000 (1.000) | 0.017 (0.017) | 0.032 (0.032) | 0.067 (0.067) |
| popmean | 0.2 | 40 | 1.000 (1.000–1.000) | 1.000 (1.000–1.000) | 0.017 (0.014–0.023) | 0.034 (0.029–0.045) | 0.061 (0.051–0.085) | 1.000, 1.000 | 1.000 (1.000) | 1.000 (1.000) | 0.017 (0.017) | 0.032 (0.032) | 0.067 (0.067) |
| live6 | 0.05 | 40 | 0.925 (0.890–0.976) | 0.906 (0.839–0.989) | 0.019 (0.017–0.025) | 0.040 (0.034–0.051) | 0.078 (0.068–0.102) | 0.923, 0.904 | 0.511 (0.350) | NaN (0.172) | 0.018 (0.017) | 0.033 (0.033) | 0.068 (0.067) |
| live6 | 0.2 | 40 | 0.904 (0.831–0.978) | 0.986 (0.977–0.995) | 0.019 (0.017–0.025) | 0.040 (0.034–0.052) | 0.078 (0.068–0.101) | 0.899, 0.986 | 0.790 (0.681) | NaN (0.398) | 0.018 (0.017) | 0.033 (0.033) | 0.068 (0.067) |

Expected-information cells are joint over the five parameters, with a 95% bootstrap over cells. "ENO, FBA alone" is the same from their 2 × 2 block only, the 16a7 script's quantity. Observed-information cells are each parameter's own curvature alone, at the mean and, in brackets, at mean + 2 SE (the most information the cells allow).

#### Bulk metabolite likelihood (200-cell means, 60 s saves)

Per pool: Σ_t (log x̄_arm − log x̄_live60)² / (2σ_b²) at σ_b = 10% (it scales as 1/σ_b²), with `live60`'s split-half noise on the same scale in the first column. The arms share cell seeds with `live60`, so until their paths part the difference is below independent noise, and the noise column is conservative. Above the noise is a real shift; above ~1 nat the reduction moves the bulk likelihood materially.

| Pool | noise | nominal | popmean | live6 |
|---|---|---|---|---|
| `M_g6p_c` | 1.56 | 83.7 | 31.3 | 2.77 |
| `M_f6p_c` | 1.59 | 83.3 | 32.6 | 2.81 |
| `M_lac__L_c` | 4.95 | 124 | 24.5 | 3.05 |
| `M_g3p_c` | 1.22 | 8.85 | 3.14 | 0.692 |
| `M_dhap_c` | 0.119 | 0.562 | 0.41 | 0.105 |
| `M_trna_chg_c` | 3.1 | 96.1 | 22.5 | 2.64 |
| `M_pi_c` | 9.2 | 101 | 26.9 | 8.83 |
| `M_3pg_c` | 11.9 | 71.3 | 20.7 | 8.52 |
| `M_2pg_c` | 10.2 | 73.7 | 15.4 | 7.49 |
| `M_gtp_c` | 7.88 | 18.2 | 13.2 | 4.86 |
| `M_fdp_c` | 0.00309 | 0.886 | 0.0948 | 0.000816 |
| `M_atp_c` | 0.649 | 0.862 | 1.06 | 0.131 |
| `M_lac__L_e` | 0.117 | 41.3 | 4.93 | 0.203 |
| `M_amp_c` | 1.94 | 1.05 | 0.548 | 0.305 |
| `M_gmp_c` | 0.927 | 0.631 | 0.265 | 0.167 |
| `M_gdp_c` | 0.19 | 0.222 | 0.288 | 0.183 |
| `M_adp_c` | 0.0907 | 0.664 | 0.115 | 0.0546 |
| **total** | **55.7** | **707** | **198** | **42.8** |

