# KIM-95: scored predictions on the reduced model (later horizon, ptsG knockdown)

**What was run.** KIM-92's joint posterior (1600 draws) was pushed forward by `predict` on
the reduced syn3A model over the full 6300 s cycle. Each run had both arms of do(ptsG):
`Knockdown` with dose 0.3 on `S_JCVISYN3A_0779` (spec/prereg/kim-92.md §4). The outputs
were scored against the true simulator at θ*, with the prior predictive beside them.

**Targets** (prereg §4):
- Q = `lactate_out` over 6300 s, Q under the knockdown (Q_kd), and their difference ΔQ;
- the 17 pools and 17 transcripts at the forecast saves 3180–6300 s. These are 53 saves each, after the fit window, which is the same 6300 s run cut at `COREA_FIT_S`.

**Target form.** Each target is a population mean: the mean over a draw's cells. The truth
is the mean over the truth's cells at θ*. No V is applied on either arm; cell statuses are
counted below.

**Provenance.**
- **Archive:** `/fred/oz022/tkimpson/infercell-kim95-archive-2026-11-19/` (README:
  `corea_predict_kim95_README.md`). Its `SHA256SUMS` lists 85 files, all of which pass
  `sha256sum -c`. The file's own SHA-256 starts `e0717b53605265b2`.
- **Jobs:**
  - sizing: 21122608–21122610 (473ecf8) and 21126961 (3b02ea4);
  - runs: truth 21154497, posterior 21154498 and prior 21154499 (arrays);
  - collect and score: 21154500.

  The runs were at baf5331. Every job ended COMPLETED 0:0 (`sacct.tsv`).
- **The tables below** are the archive's `score.md`, written by
  `corea_predict_kim95_score.jl` at 2333ae2. Rerunning the scorer on a copy of the
  archive's TSVs gives a byte-identical `score.md` and `score_targets.tsv`.
- **Cost:** 5,447 CPU-h of cell work: posterior 4,260, prior 732, truth 456 (`collect.md`).

## Sizes, frozen before any truth was read

M0's declared rules (KIM-69) are:
- the MC term is at most 10% of the predictive spread, which gives n_cells = ⌈9σ²/Var_θ⌉;
- the truth's SE is at most 10% of the predictive SD, which gives n_truth = ⌈100σ²/Var_θ⌉.

Var_θ is the spread of μ̂ across draws that comes from θ, and σ² is one cell's variance.

- **The first sizing pass did not resolve the posterior's Var_θ** (job 21122610, 40 draws × 8 cells).
  - The ratio of var(μ̂) across draws to its MC term was 0.81 on Q, and 0.96 at the median of the 1,802 per-save targets. A ratio of 1 means the draws' spread is all cell noise.
  - The prior's Var_θ was resolved, and its rule needed at most 202 cells per draw.
- **The resolving pass** (job 21126961, 40 draws × 100 cells) gave:
  - Q_kd: Var_θ/σ² = 0.0093, a ratio of 1.93 ± 0.23;
  - Q: still unresolved, with a ratio of 0.99;
  - ΔQ: still unresolved, with a ratio of 0.81.
- **Tom's call (2026-11-18): size on Q_kd.** The frozen sizes are:
  - 100 draws per arm;
  - 973 cells per posterior draw and 200 per prior draw;
  - n_truth = 10,810, for each knockdown arm.

  They are frozen at baf5331.

The full run has more power than either sizing pass. With it, Q's posterior Var_θ is resolved: its MC share is 0.21, which puts Var_θ/σ² at about 0.004. ΔQ's is not resolved: its MC share is 1.09.

## Q and the knockdown

| target | truth ± SE | posterior 5–50–95% | prior 5–50–95% | post 50/90 | prior 50/90 | CRPS post | CRPS prior | skill | SE/SD | MC share post | MC share prior |
|---|---|---|---|---|---|---|---|---|---|---|---|
| Q | 1.498e+07 ± 1.828e+04 | 1.481e+07 – 1.501e+07 – 1.524e+07 | 6.264e+06 – 1.138e+07 – 1.738e+07 | in/in | out/in | 3.293e+04 | 2.207e+06 | 0.9851 | 0.143 | 0.2108 | 0.0008168 |
| Q_kd | 1.504e+07 ± 1.755e+04 | 1.49e+07 – 1.51e+07 – 1.532e+07 | 6.191e+06 – 1.106e+07 – 1.708e+07 | in/in | out/in | 4.155e+04 | 2.386e+06 | 0.9826 | 0.1389 | 0.1917 | 0.0007377 |
| dQ | 6.423e+04 ± 2.32e+04 | -2.714e+04 – 8.801e+04 – 2.183e+05 | -5.645e+05 – -1.707e+05 – 2.331e+05 | in/in | out/in | 1.857e+04 | 1.417e+05 | 0.869 | 0.3276 | 1.085 | 0.1916 |

- **Q and Q_kd.**
  - The posterior's 90% interval is ±1.4% of Q. The truth sits 0.26 posterior SD below the posterior median on Q, and 0.48 posterior SD below it on Q_kd. Both are inside the 50% intervals.
  - The prior's median is 24–27% low, and its 90% interval spans a factor of 2.8.
  - Skill is 0.98 on both. It is near 1 mostly because the prior is so wide, so the main finding is that the narrow posterior interval holds the truth.
- **ΔQ.**
  - The truth is +6.4e4 ± 2.3e4, about 0.4% of Q.
  - The posterior's interval runs from −2.7e4 to +2.2e5, and all of that width is cell noise at 973 cells (MC share 1.09). So the posterior shows that the knockdown effect is small next to one draw's cell noise, and it cannot resolve the effect beyond that.
  - The prior's median, −1.7e5, has the wrong sign, and its interval only just reaches the truth.
- **The truth's SE is larger than the rule allows.** The rule asks for SE/SD ≤ 0.1. It is 0.14 on Q, 0.14 on Q_kd and 0.33 on ΔQ, because the posterior is narrower than the sizing could see. So the hits on these three are read with a truth that is not sharp against the predictive.

## The later horizon (3180–6300 s)

| family | species | saves | post hit 50 | post hit 90 | prior hit 50 | prior hit 90 | median species skill | SE/SD median | SE/SD max | MC share post (median) | MC share prior (median) |
|---|---|---|---|---|---|---|---|---|---|---|---|
| pool | 17 | 901 | 34% | 92% | 29% | 73% | 0.9452 | 0.1856 | 0.3413 | 0.3717 | 0.004651 |
| transcript | 17 | 901 | 80% | 100% | 64% | 85% | 0.9514 | 0.2474 | 0.3398 | 0.6753 | 0.007484 |

The per-species table is in `score.md`. Saves within a species are strongly correlated, so a family's hit fraction rests on about 17 species, not 901 independent checks.

- **Transcripts.**
  - Every save is inside the posterior's 90% interval, and 80% are inside its 50% interval. At each species' median save the truth is within 0.63 posterior SD.
  - Their posterior predictive is about two-thirds cell noise (median MC share 0.68). That widens the intervals, so the high hit fractions overstate how sharp the θ-predictive is.
  - Skill is 0.89–0.99 for every transcript.
- **Pools.**
  - 92% of saves are inside the posterior's 90% interval, but only 34% are inside its 50% interval. That is a shift, not scatter. At each species' median save the truth sits about 1–1.4 posterior SD from the posterior median, as follows.
  - **Triose pools read high, with the truth above the posterior.** G3P is +1.39 SD (MC share 0.11) and DHAP +1.23 SD (MC share 0.02). Neither has a single save inside the 50% interval. Their skill is the lowest of all species: 0.54 and 0.64.
  - **Nucleotides split by charge.**
    - ATP, ADP and GDP read low (truth below): −1.33, −1.00 and −1.40 SD.
    - AMP and GMP read high: +1.14 and +0.97 SD.
  - **Upper glycolysis reads high:** G6P, F6P and FDP at +0.9 to +1.0 SD.
  - **The triose and FDP errors are θ-driven:** their MC shares are small (0.02–0.11), so these are biases of the posterior, not cell noise. The nucleotides' MC shares are larger (0.33–0.62).
  - Skill is still 0.90–0.99 everywhere except the trioses.
- **Reading against the pilot.** The pilot's promoters sat low (mean −0.81 SD; Decision log 2026-11-17), yet the transcripts here show no comparable shift. The pool shifts, especially the trioses and the energy charge, are where the posterior's error shows in the forecast.

## Cell statuses (`status.tsv`)

| arm | cells per knockdown arm | clipped | starved | failed | clean |
|---|---|---|---|---|---|
| posterior | 97,300 | 99.3% / 99.5% | 656 / 464 | 0 | 0 |
| prior | 20,000 | 74.9% / 74.7% | 1,815 / 1,512 | 0 | 3,197 / 3,544 |
| truth | 10,810 | 99.3% / 99.3% | 77 / 78 | 0 | 0 |

Each pair is baseline / knockdown. `clipped` is the GTP clip that KIM-42 reported. Starved cells stay in μ̂, because there is no V.

## For F5

F5's claim is that "the posterior's predictions cover / miss the truth by this much, against the prior predictive". The evidence for it:
- **They cover:** Q and Q_kd within ±1.4% intervals, every transcript save at 90%, and 92% of pool saves at 90%.
- **They miss:** the triose pools by about 1.3 SD and the energy charge by about 1–1.4 SD, which leaves only 34% of pool saves inside the 50% interval.
- **Against the prior:** the CRPS skill is 0.98 on Q, 0.87 on ΔQ, and a median of 0.95 across species.
- **Not resolved:** ΔQ under the knockdown is not resolved by the posterior beyond "small next to cell noise".
