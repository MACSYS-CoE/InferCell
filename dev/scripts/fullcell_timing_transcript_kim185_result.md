# KIM-185: the transcript factor on 455 genes

Full-cell timing 3 (Linear M11). This times **stage 1 only**, the gene-expression block, of the
full syn3A model by the data plug-in. It uses the reduced model's design and src's exact
transcript factor.

- **The block is a stand-in:** the reduced model's gene law applied to the full model's 455
  genes. It has one gene copy, the reduced death law, and pools fixed from the panel.
- **Stage 2 is outside the threshold.** That is the bulk block; it was 198.7 of the plug-in's
  365.2 CPU-h in KIM-116.

**PASS against the 2026-11-26 threshold, provided a per-gene gradient replaces ForwardDiff.gradient over all coordinates in both the sampler and the metric stage (hessian_fd).** The block projects to **539 CPU-h**
per posterior on the pre-named headline row, against at most 5,000 (Decision log 2026-11-26,
closing KIM-176). No Linear issue owns the per-gene-gradient work yet.

- **One serial dual pass takes 53.61 s.** That is ForwardDiff with `Chunk{2}` per gene, over 456
  coordinates, 200 cells and 52 windows. The density takes 28.41 s.
- **Threading over genes** brings the dual pass to 7.28 s at 16 workers (6.8×), and no lower at
  32.
- **The condition matters.** Naive ForwardDiff over all 456 coordinates projects to 74,158 CPU-h,
  a FAIL. In the metric stage alone it costs about 1,869 CPU-h (912 gradients × 38 passes ×
  194.11 s).
- **The sensitivity range runs from 787 CPU-h to 1,192–1,800 CPU-h.**
  - 787 CPU-h is the ratio to src's 17-gene gradient, a margin of about 6×.
  - 1,192 CPU-h adds NUTS's (D/18)^¼ step-size scaling. About 1,800 CPU-h adds the chains'
    1.50× per-step gap on top (not pre-named). That is a margin of about 3×.
  - The 1.50× gap is mostly accounted for; see the projection.
- **The step count at 456 coordinates was not measured.** The PASS holds up to about 9.3×
  KIM-116's leapfrog steps per iteration, or 6.2× with the 1.50× gap.
- **Affordable is not right.** A pass says the route can be paid for, not that it is accurate.
  KIM-116 found the plug-in biased on the reduced model (gate 3).

## Run

| stage | Slurm job | node | code | CPU-h |
|---|---|---|---|---|
| data (item 1) | 21664771 | dave5 | 7b06c43 | 0.14 |
| serial (item 2) | 21665221 | dave23 | 7b06c43 (output headers read 49426cf; see item 2) | 0.51 |
| threaded (item 3) | 21665222 | dave104, 32 CPUs | 7b06c43 (output headers read 49426cf) | 9.40 |
| projection | 21666966 | dave28 | 817a648 | 0.01 |

- **Hardware and software.** Every node is an AMD EPYC 7543 32-Core Processor on milan, running
  Julia 1.10.5 and ForwardDiff 0.10.39. BLAS ran on 1 thread.
- **Outputs:** `/fred/oz022/tkimpson/infercell-kim185-2026-11-26/run/`. It holds `data`, `serial`,
  `threaded` and `projection`, each as `.jls` and `.md`, plus `jobs.tsv`, `COMMIT` and `logs/`.
  `SHA256SUMS` in that directory covers its 14 files; its own sha256 is
  `6cf258f39f21d4b6e21e446525f8e6c454666ea713c4c651f56af8e7591be74f`. `data.jls` serialises src
  types, so it may not deserialise at a later commit; archive the directory before KIM-187 relies
  on it.
- **The committed driver differs from the one the jobs ran,** but only by the review's guards
  (f48dd2e).
  - The jobs ran 7b06c43; the projection job ran 817a648, which has the same code.
  - The guards are: the data stage's checks 1 and 2 now raise an error when they fail; the
    projection now refuses inputs that disagree, and prints SMOKE; and a header comment is
    corrected.
- **Smokes** (`smoke1–3/`) are plumbing only, never a result.
  - Two failed for technical reasons, each fixed in its own commit: the rebuild lookup (21663573,
    fixed in a6771d4) and the pool slice (21663913, fixed in 8198913).
  - So did a parse of an empty `WORKERS` (21664483, fixed in 7b06c43).
  - Their dependants 21663574–76 and 21663914–16 are left pending. 21663576 and 21663916 show
    "Dependency" and the other four show DependencyNeverSatisfied. They need a `scancel`.
- **Total compute:** 10.7 CPU-h, smokes included.
- **Rerun** from the repo root of a clean worktree at a pushed commit:
  `dev/scripts/fullcell_timing_transcript_kim185.sh {data|serial|threaded|project|all} <outdir>`.
  The `data` stage reads the KIM-92 archive (`infercell-kim92-archive-2026-11-17`). The `project`
  stage reads the KIM-116 archive (`infercell-kim116-archive-2026-11-26`).

## What was built

`src/` is unchanged. The driver is `dev/scripts/fullcell_timing_transcript_kim185.jl` and its
launcher is `fullcell_timing_transcript_kim185.sh`.

- **The factor.** src's `TranscriptFactor` and `transcript_loglik` are called as they are. The
  plug-in in src, `plugin_loglik`, is tied to the reduced model's 17-gene `MixedSetup`. So the
  455-gene case is assembled in the script from 455 independent genes:
  - each gene gets one single-gene `CoreATranscription`, whose `rate_constants` is the reduced
    model's rate law;
  - each gene and cell gets one single-gene `TranscriptFactor`;
  - the log-likelihood is the sum over genes of `transcript_loglik` over that gene's 200 factors.
- **The rates** are the plug-in's. Gene g's constant at each switch is the rate law at its
  promoter strength S_g, with the pools fixed from the observed panel. Those pools are the reduced
  model's plug-in states on KIM-92's dataset (`panel_plugin_states` at θ_ref): ATP and GTP from
  the panel, and CTP and UTP from the chemostats.
- **The block per gene** is the reduced model's: log S_g, plus one log krnadeg shared by all
  genes, with μ_g = krnadeg / len_g. **455 genes give 456 free coordinates.**
- **The priors** are the reduced model's. log S_g ~ Normal(log S_ref,g, log 2), and krnadeg has its
  reduced-model prior.
- **The support** is sized as `mixed_setup` sizes it: kmax = TURNOVER_CEILING / (len − 1), and
  μmin is the 1e-6 quantile of the krnadeg prior divided by len.
- **Checks against src**, at KIM-92's truth:
  - On the reduced model's own 17 genes and 200 cells, the per-gene construction reproduces src's
    plug-in stage-1 density (`stage_logdensity` with `plugin_loglik`, 18 coordinates) and its
    gradient. The numbers are in item 2.
  - A single-gene module's `rate_constants` reproduces the 17-gene rebuild
    (`transcription_rates_from_states`) at the truth, with a max relative difference of 0.

## Item 1: synthetic transcripts for 455 genes

Data job **21664771** at commit 7b06c43, on dave5 (AMD EPYC 7543 32-Core Processor), Julia
1.10.5. The output is `data.jls` and `data.md` in `/fred/oz022/tkimpson/infercell-kim185-2026-11-26/run/`.

- **Genes.** The source is syn3a-ensemble `data/model/cme/`, at commit
  `2075809cbb825140cb31a7bec8bc3529770c8e13`.
  - `reactions.csv` has sha256 `63fbeca180fe8bdebec646df0ddbffb0149958fc2ff9fc3e6ec8803cd24c6aaf`.
  - `species.csv` has sha256 `c1acb80878e8a500a86e0e1d9e195a491b46a01c337788c0a6554df19cbe6ec9`.
  - Exactly 455 transcription reactions make an mRNA. Each gives the gene's length (its
    `ATP_trsc` increment), its four base counts, and its first-minute rate constant. The mRNA's
    initial count comes from procid 1.
  - The first two bases are taken as AU for every gene. That holds for all 17 reduced genes, but
    not for 15 of the 455: syn3A.gb starts 12 genes with GTG (GU) and 3 with TTG (UU). The review
    of PR #206 put the largest effect on k at about 2.4e-5 relative, so the timing does not see it.
  - S_ref,g is the first-minute rate divided by the same law at S = 1. That law uses the published
    base mapping and the first-minute concentrations, which is how `MinCell_CMEODE.py` builds
    minute 1.
- **Check 1.** The 17 reduced genes, read from the full model's tables, have the extract's
  lengths and base counts. Their S_ref matches ptn_count / 180 to a max relative difference of
  1.15e-11.
- **The design** is the reduced model's (KIM-92's dataset, `corea_pilot_kim92.jl`: `ncells = 200`,
  `fit_to = COREA_FIT_S`).
  - There are 200 cells and 52 windows, with saves every 60 s from 60 s to 3120 s.
  - Every cell starts at the same initial counts, as the reduced model's cells share theirs.
- **The truth.**
  - log S_g is drawn from its prior with seed 18500001. The median S / S_ref is 0.986.
  - krnadeg = 4.26926, which is KIM-92's truth.
- **The counts are sampled exactly.** Each window runs 59 s at k[:, m], then 1 s at k[:, m + 1].
  A linear birth-death step of length Δ from count n is Binomial(n, e^{−μΔ}) + Poisson(k (1 − e^{−μΔ}) / μ).

**Array sizes:**

| array | size |
|---|---|
| counts `y` | 200 cells × 455 genes × 52 windows (`200×455×52 Array{Int64, 3}`) |
| initial counts `y0` | 455 |
| birth rates `k` | 455 × 53 (genes × windows + 1) |
| death rates `μ` | 455 |
| free coordinates | 456 |

What the genes and counts look like:

| quantity | value |
|---|---|
| gene length | median 852 nt, range 114–5388 |
| S_ref | median 0.617, range 0.0556–7.53 |
| birth rate k at the truth | median 6.58e-4 /s, range 1.02e-5–0.0602 |
| counts | mean 0.309, max 17; y0 mean 0.295, max 5 |
| cap over genes × cells | min 20, median 50, 90% 87, 99% 140, max 172 |

The full model's own first-order mRNA degradation constant is 0.821× the reduced model's law at
KIM-92's krnadeg, for every gene. This is reported only; the block uses the reduced model's law.

## Item 2: one density and one dual pass, serial

Serial job **21665221**, on dave23 (AMD EPYC 7543 32-Core Processor), Julia 1.10.5,
ForwardDiff 0.10.39, 1 thread, BLAS 1 thread. The output is `serial.jls` and `serial.md`.

- **The commit.** The job was submitted at 7b06c43, and the outdir's `COMMIT` file records that.
  The `.md` and `.jls` headers read 49426cf instead, because `pilot_common.jl` runs
  `git rev-parse` when Julia loads it. 49426cf was committed 26 s after both timing jobs started
  and adds only this note, so the code is 7b06c43's.
- **A hazard this ran into.** The jobs read the live worktree while commits landed, which the
  project's notes warn against. Here the only change was the note.
- **How each call ran.** The 455-gene calls ran 5 times after compiling on a one-cell copy. The
  17-gene harness compiled at full size. A full GC ran before each call.

**The harness, checked on the reduced model** (17 genes, 200 cells, 52 windows, 18 coordinates,
at KIM-92's truth, same node). The per-gene construction reproduces src's plug-in stage density
to a relative 2.36e-16 and its gradient to 3.26e-15 (max |Δ| over max |g|).

| call | AD | median, s |
|---|---|---|
| src plug-in density | — | 0.96 |
| src plug-in gradient | `ForwardDiff.gradient`, default chunk 9 (2 passes) | **11.26** |
| per-gene density | — | 1.09 |
| per-gene dual pass | ForwardDiff, `Chunk{2}` per gene | 2.48 |

KIM-116's profile (job 21196497, dave14) measured 11.26 s for the src gradient, and this node
reproduces the gradient. It does not reproduce the density: 0.96 s here against KIM-116's 1.18 s,
19% lower. The per-gene dual pass is 4.5× faster on the reduced model's own block.

**455 genes** (200 cells, 52 windows, 456 coordinates, at the truth):

- **Check against the whole factor.** The per-gene sum equals the whole factor (one
  `TranscriptFactor` per cell over all 455 genes) to a relative 2.24e-15. The per-gene gradient
  equals ForwardDiff's on the whole factor, on genes 1–6 plus krnadeg, to 8.09e-17.
- **At the truth,** the log density is −1557198.88 and the gradient norm is 742.45.

| call | AD | each, s | median, s | GiB allocated (median) | GC s (median) |
|---|---|---|---|---|---|
| density | — | 28.15, 28.41, 28.42, 28.39, 28.41 | **28.41** | 37.19 | 4.64 |
| dual pass (value + gradient, 456 coordinates) | ForwardDiff, `Chunk{2}` per gene | 52.87, 53.69, 53.33, 53.75, 53.61 | **53.61** | 104.86 | 8.31 |
| naive: one chunk pass | ForwardDiff, Dual width 12 over all 456 coordinates, whole factor | 191.92, 194.22, 194.11 (3 repeats) | 194.11 | 431.01 | 10.18 |

- **Seconds per dual pass: 53.61 serial.** That is one ForwardDiff pass per gene with two
  partials, which gives the value and the full 456-coordinate gradient.
- **The naive gradient** is 38 chunk passes × 194.11 s = 7,376 s, about 2 h. It is a
  multiplication, not a timed gradient. It is 138× the per-gene dual pass, because every chunk
  carries 12 partials through every gene, though each gene depends on only 2 of them.
- **Cost is uneven over genes.** The per-gene dual pass has a median of 0.077 s and a max of 1.54 s.
  The 10% most expensive genes take 37% of the total.
- **Against the reduced model.** The dual pass is 21.6× the 17-gene per-gene pass (2.48 s), for
  26.8× the genes. It is 4.8× src's 17-gene gradient (11.26 s).

## Item 3: the same, threaded over genes

Threaded job **21665222**, on dave104 (AMD EPYC 7543 32-Core Processor), Julia 1.10.5,
ForwardDiff 0.10.39. It had 32 threads, `--gcthreads=16,1` (17 GC threads, 1 sweep) and BLAS
on 1 thread. It was submitted at 7b06c43. Its output headers read 49426cf for the reason given in
item 2.
The output is `threaded.jls` and `threaded.md`.

- **How it is threaded.** W tasks pull genes off a shared atomic counter, heaviest genes first
  (by summed cap). Each gene's dual pass is ForwardDiff with `Chunk{2}`.
- **Results do not depend on W.** The results land in gene order and are summed in that order.
  At every W the value and gradient equal the one-worker result bit for bit; the script checks
  this.
- **Timing.** Each call ran 5 times, with a full GC before each. The speed-up is against W = 1 in
  the same job, on the same node.

| workers | density median, s | speed-up | dual pass median, s | speed-up | efficiency | GC s in the dual pass (median) |
|---|---|---|---|---|---|---|
| 1 | 25.37 | 1.00 | 49.64 | 1.00 | 1.00 | 3.97 |
| 4 | 6.66 | 3.81 | 14.92 | 3.33 | 0.83 | 3.10 |
| 16 | 2.91 | 8.70 | **7.28** | **6.82** | 0.43 | 2.48 |
| 32 | 3.89 | 6.53 | 7.65 | 6.48 | 0.20 | 2.78 |

- **The speed-up saturates** at about 6.8× by 16 workers, and 32 workers is no faster.
- **GC is a large share at high W, but it does not explain all the loss.** At W = 16:
  - GC is 2.48 s of the 7.28 s;
  - ideal scaling from W = 1 would give 3.10 s;
  - about 1.7 s is beyond both, and its cause is unexplained.
- **At W = 32, the gene is the floor.** The ideal time per worker (1.55 s) equals the heaviest
  gene's dual pass (1.54 s), so per-gene granularity limits the speed-up there.
- **Nothing was tuned.**
- **The serial job's 8% excess is GC marking, not the node.** Serial (dave23) took 53.61 s and
  W = 1 here took 49.64 s.
  - Their non-GC time agrees within 1%: 45.30 s against 45.67 s.
  - Their GC time does not: 8.31 s with 1 GC thread, against 3.97 s with 17.
  - So the threaded CPU-h rows at W = 1 and 4 count W CPUs but not the GC mark threads on the
    job's other cores. They are understated.

## The projection: the gene-expression block's cost per posterior

The projection job is **21666966** at 817a648. It reads only `serial.jls`, `threaded.jls` and
KIM-116's archive, and its output is `projection.md`.

**Gradients per posterior** are counted from the leapfrog steps (`n_steps`) in KIM-116's plug-in
stage-1 chains and extensions (`infercell-kim116-archive-2026-11-26/run/plugin/`). This count is
slightly low; see the 1.50× gap below.

| chain | iterations (chain + extension) | leapfrog steps (chain + extension) | max tree depth | hours |
|---|---|---|---|---|
| 1 | 400 + 200 | 5564 + 2760 | 8 | 42.9 |
| 2 | 400 + 200 | 6007 + 2856 | 7 | 42.3 |
| 3 | 400 + 200 | 6039 + 2816 | 7 | 38.8 |
| 4 | 400 + 200 | 6280 + 2960 | 8 | 41.0 |

- **Leapfrog steps:** 35,282 in total, over 165.1 h of sampling.
- **The warm metric:** a central-difference Hessian, which took 2·18 = 36 gradients at 18
  coordinates. At 456 coordinates it takes 912.
- **So one posterior is 36,194 gradients.**
- **The 1.50× gap per step.** The chains spent 16.85 s per leapfrog step, against the profile's
  11.26 s per gradient. The archive accounts for most of it. The review of PR #206 traced this; it
  was not re-measured here.
  - **`n_steps` undercounts by +6.6% to +10%.** AdvancedHMC's refresh recomputes one gradient per
    iteration, which adds 2,400 gradients: +6.6% (verified in AdvancedHMC 0.8.3's source). Turing
    may also re-evaluate the density once per iteration, which would make it about +10% (not
    verified). At 455 genes the headline would then be about 575–594 CPU-h.
  - **Suspended time was counted.** Each of the four chain jobs (21208681–84) was Suspended for
    52 min, which `@elapsed` counted: about 3.5 h of the 165.1 h.
  - **The rest is a slower gradient: 1.37×.** The sampler-free metric job 21208680 ran its 36
    gradients at about 15.4 s each, 1.37× the 11.26 s profile. So the slowdown is not sampler
    overhead. Its cause (the memory request, or node load) was not tested.
  - **Together these close the gap.** (165.1 − 3.5) h over 37,682 gradients is 15.4 s a
    gradient, the same 1.37×.

| projection | CPU-h |
|---|---|
| **headline (fixed before the timing was read):** 36,194 gradients × 53.61 s, the serial per-gene dual pass | **539** |
| the same, with the leapfrog steps × (456 / 18)^¼ = 2.24 for NUTS's step size at higher dimension | 1,192 |
| the ratio to src's 17-gene gradient: KIM-116's stage-1 165.4 CPU-h × (53.61 s ÷ 11.26 s, on the same node). It carries KIM-116's per-gradient slowdown and undercount, and is 1.46× the headline | 787 |
| both adjustments: the dimension row × 1.50 (not pre-named) | about 1,790 |
| naive AD everywhere: 36,194 gradients × 38 chunk passes × 194.11 s | 74,158 |
| naive AD in the metric stage only: 912 gradients × 38 × 194.11 s, added to any row above | 1,869 |
| threaded, 4 workers: 36,194 × 14.92 s × 4 CPUs; wall 37.5 h per chain of 9,048 gradients | 600 |
| threaded, 16 workers: 36,194 × 7.28 s × 16 CPUs; wall 18.3 h per chain | 1,171 |
| threaded, 32 workers: 36,194 × 7.65 s × 32 CPUs; wall 19.2 h per chain | 2,463 |

The "both adjustments" and "metric stage only" rows are this note's arithmetic, not rows of
`projection.md`.

**Verdict: PASS against the 2026-11-26 threshold, provided a per-gene gradient replaces ForwardDiff.gradient over all coordinates in both the sampler and the metric stage (hessian_fd).**

- **The threshold** is at most 5,000 CPU-h for the gene-expression block over 455 genes by the data
  plug-in (Decision log 2026-11-26, closing KIM-176).
- **The verdict** is on the pre-named headline, 539 CPU-h.
- **The sensitivity range** runs from 787 CPU-h (the ratio row, a margin of about 6×) to
  1,192–1,800 CPU-h with dimension scaling (a margin of about 3×).
- **Naive AD.** Every row except naive AD everywhere is under 5,000 CPU-h. Leaving naive AD in the
  metric stage alone, on top of both adjustments, gives about 3,700 CPU-h, still under.
- **No Linear issue owns the per-gene-gradient work yet.**

**Break-even on the unmeasured step count.**

- The PASS holds up to about 9.3× KIM-116's leapfrog steps per iteration (14.7 → about 136).
  With the 1.50× gap it holds up to about 6.2× (about 91).
- That many steps is unlikely here: only krnadeg couples the genes, and the metric is warm and
  diagonal. But max depth 10 (up to 1,023 steps) allows it.
- **Cost per gradient also depends on krnadeg.** The review estimated ×1.42 at krnadeg × 2 (one
  prior SD) and ×2.17 at × 4. These are not Slurm measurements. On top of about 1,790 the × 4 case
  is about 3,900, still under.

**Wall-clock per chain** (9,048 gradients):

| case | serial | 16 workers |
|---|---|---|
| headline | 135 h | 18 h |
| with the 1.50× gap | about 202 h | about 27 h |
| with the dimension factor | about 300 h | about 41 h |
| with both | about 450 h | about 61 h |

The serial column is dave23's time and the 16-worker column is dave104's, so their ratio (7.4×)
is not the one-node speed-up (6.82×).

- **milan's limit is 7 days (168 h).** KIM-116 ran each chain and its extension as separate jobs,
  so the limit applies per job.
- **The headline and the 1.50× case fit, serial.** Under the headline, the longest chain job is
  about 94 h and its extension about 44 h. With the 1.50× gap they are about 140 h and 66 h.
- **The dimension factor does not fit.** It makes the longest chain job about 210 h serial. With
  it, a serial chain needs checkpointing across jobs, or threading.

### What the projection assumes

1. **The chain design is KIM-116's plug-in stage 1:**
   - 4 chains;
   - NUTS(200 adaptation, 0.8), max depth 10;
   - 200 kept draws per chain, plus the one extension of 200;
   - the warm diagonal metric.

   It also assumes that one extension is enough. A worst-R̂ / min-ESS gate over 456 coordinates is
   stricter than over 18. Each extra extension adds about a third to the step count, which is
   inside the break-even.
2. **The same leapfrog steps per iteration at 456 coordinates as at 18.** This was not measured.
   - The (D/18)^¼ row is the usual scaling of the step size with dimension.
   - A posterior with stronger correlations or worse conditioning than the reduced model's would
     need more steps. Nothing here measures that.
3. **The cost per gradient is the serial dual pass at the truth,** on one EPYC 7543 core.
   - The caps are fixed by the data and the support. Elsewhere in the posterior only the term
     counts move (through Λ = k + μ·cap).
   - The headline leaves out KIM-116's 1.37× per-gradient slowdown and the `n_steps` undercount.
     The ratio row carries both.
   - Cost per gradient also depends on krnadeg; see the break-even.
4. **The AD is per gene.** Each gene's term depends only on its own promoter and on krnadeg, so
   a dual pass with two partials per gene is exact.
   - Neither KIM-116's sampler nor src's `plugin_loglik` does this today: they run
     `ForwardDiff.gradient` over all coordinates.
   - At 456 coordinates that is 38 chunk passes and 74,158 CPU-h, which fails the threshold.
   - The verdict therefore depends on wiring a per-gene gradient into the sampler and into the
     metric stage. KIM-116's `hessian_fd` calls `ForwardDiff.gradient` over all coordinates, and
     naive AD there costs about 1,869 CPU-h.
   - No Linear issue owns that work yet.
5. **Stage 2 is excluded.** That is the bulk/metabolic block, at the frozen threshold's scope.
6. **The model is the reduced model's gene block applied to 455 genes:**
   - one gene copy, no replication over the window;
   - death rate krnadeg / len (the full model's own constant is 0.821× that at the truth);
   - first two bases taken as AU (12 genes really start GU and 3 UU; the effect on k is at most
     about 2.4e-5);
   - pools fixed from the panel, as the plug-in fixes them.
7. **The data design is the reduced model's:** 200 cells, 52 windows, every transcript observed
   exactly.
   - Cost grows about linearly with windows.
   - It grows strongly sublinearly with cells. `transcript_loglik` shares a window distribution
     across cells with the same (gene, window, starting count, cap).
   - So the cost depends on the count distribution: the mean count is 0.31 here. The full model,
     with replication and growth, would give higher counts and more distinct starts.

### What this does not say

- **Affordable is not right.** A pass here does not make the route right; that is gate 3
  (project description). KIM-116 found the plug-in biased on the reduced model: promoter 0607 sat
  +0.65 joint SD off, and 16 of 17 promoters shifted up. This timing does not touch that.
- **Threading.** The threaded figures are as measured. GC takes about a third of the time at
  16 workers, and nothing was tuned.
- **The per-gene code is a timing harness, not a sampler's density.**
  - `gene_ll` calls `transcript_loglik` directly.
  - Outside the support it throws an error.
  - Outside the ±6 SD box it stays finite, where src's plug-in returns −Inf.
  - A per-gene gradient for the sampler must keep src's −Inf handling.
