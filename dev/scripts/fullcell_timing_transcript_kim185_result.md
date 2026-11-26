# KIM-185: the transcript factor on 455 genes

Full-cell timing 3 (Linear M11). This times the gene-expression block of the full syn3A model by
the data plug-in: 455 genes, the reduced model's design, and src's exact transcript factor.

**PASS against the 2026-11-26 threshold.** The block projects to **539 CPU-h** per posterior,
against at most 5,000 (Decision log 2026-11-26, closing KIM-176).

- **One serial dual pass takes 53.61 s.** That is ForwardDiff with `Chunk{2}` per gene, over 456
  coordinates, 200 cells and 52 windows. The density takes 28.41 s.
- **Threading over genes** brings the dual pass to 7.28 s at 16 workers (6.8×), and no lower at
  32.
- **The verdict rests on the per-gene gradient.** The sampler would need one that it does not have
  today. Naive ForwardDiff over all 456 coordinates projects to 74,158 CPU-h.
- **The adjustments the headline leaves out keep it under the threshold.** NUTS's step-size scaling
  with dimension and KIM-116's sampler overhead, together, give about 1,790 CPU-h.

## Run

| stage | Slurm job | node | code | CPU-h |
|---|---|---|---|---|
| data (item 1) | 21664771 | dave5 | 7b06c43 | 0.14 |
| serial (item 2) | 21665221 | dave23 | 7b06c43 (`COMMIT` 49426cf; see item 2) | 0.51 |
| threaded (item 3) | 21665222 | dave104, 32 CPUs | 7b06c43 (`COMMIT` 49426cf) | 9.40 |
| projection | 21666966 | dave28 | 817a648 | 0.01 |

- **Hardware and software.** Every node is an AMD EPYC 7543 32-Core Processor on milan, running
  Julia 1.10.5 and ForwardDiff 0.10.39. BLAS ran on 1 thread.
- **Outputs:** `/fred/oz022/tkimpson/infercell-kim185-2026-11-26/run/`. It holds `data`, `serial`,
  `threaded` and `projection`, each as `.jls` and `.md`, plus `jobs.tsv`, `COMMIT` and `logs/`.
- **Smokes** (`smoke1–3/`) are plumbing only, never a result.
  - Two failed for technical reasons, each fixed in its own commit: the rebuild lookup (21663573,
    fixed in a6771d4) and the pool slice (21663913, fixed in 8198913).
  - So did a parse of an empty `WORKERS` (21664483, fixed in 7b06c43).
  - Their dependants 21663574–76 and 21663914–16 are left pending with
    DependencyNeverSatisfied.
- **Total compute:** 10.7 CPU-h, smokes included.
- **Rerun** from a clean worktree at a pushed commit:
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
  - The first two bases are taken as AU for every gene (the ATG start codon). All 17 reduced genes
    are AU.
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
ForwardDiff 0.10.39, 1 thread, BLAS 1 thread. It was submitted at 7b06c43. The job's `COMMIT`
reads 49426cf, the branch head when it started; that commit adds only this note, so the code is
7b06c43's. The output is `serial.jls` and `serial.md`.

Each call ran 5 times after compiling on a one-cell copy, with a full GC before each.

**The harness, checked on the reduced model** (17 genes, 200 cells, 52 windows, 18 coordinates,
at KIM-92's truth, same node). The per-gene construction reproduces src's plug-in stage density
to a relative 2.36e-16 and its gradient to 3.26e-15 (max |Δ| over max |g|).

| call | AD | median, s |
|---|---|---|
| src plug-in density | — | 0.96 |
| src plug-in gradient | `ForwardDiff.gradient`, default chunk 9 (2 passes) | **11.26** |
| per-gene density | — | 1.09 |
| per-gene dual pass | ForwardDiff, `Chunk{2}` per gene | 2.48 |

KIM-116's profile (job 21196497, dave14) measured 11.26 s for the src gradient. So this node and
this harness reproduce what KIM-116 measured. The per-gene dual pass is 4.5× faster on the
reduced model's own block.

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
on 1 thread. It was submitted at 7b06c43; its `COMMIT` reads 49426cf for the reason given in item 2.
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
- **GC takes a large share at high W.** At W = 16, GC is 2.5 s of the 7.3 s, and GC stops every
  thread. The dual pass allocates heavily, so the scaling is bound by allocation and GC, not by the
  number of genes. This is a measurement, not a tuned figure. Fewer allocations in
  `_propagate`/`_uniform_step`, or a process per group of genes, would be the levers (KIM-194).
- **Node-to-node spread.** One worker here took 49.64 s, against 53.61 s for the serial job on
  dave23 (8%).

## The projection: the gene-expression block's cost per posterior

The projection job is **21666966** at 817a648. It reads only `serial.jls`, `threaded.jls` and
KIM-116's archive, and its output is `projection.md`.

**Gradients per posterior** are counted, not inferred. They are the leapfrog steps (`n_steps`) in
KIM-116's plug-in stage-1 chains and extensions (`infercell-kim116-archive-2026-11-26/run/plugin/`):

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
- **Sampler overhead in the chains:** 16.85 s per leapfrog step, against the profile's 11.26 s per
  gradient. That is 1.50× more.

| projection | CPU-h |
|---|---|
| **headline (fixed before the timing was read):** 36,194 gradients × 53.61 s, the serial per-gene dual pass | **539** |
| the same, with the leapfrog steps × (456 / 18)^¼ = 2.24 for NUTS's step size at higher dimension | 1,192 |
| KIM-116's stage-1 165.4 CPU-h × (53.61 s ÷ 11.26 s, src's 17-gene gradient on the same node). This carries the 1.50× sampler overhead | 787 |
| both adjustments: the dimension row × 1.50 | about 1,790 |
| naive AD: 36,194 gradients × 38 chunk passes × 194.11 s | 74,158 |
| threaded, 4 workers: 36,194 × 14.92 s × 4 CPUs; wall 37.5 h per chain of 9,048 gradients | 600 |
| threaded, 16 workers: 36,194 × 7.28 s × 16 CPUs; wall 18.3 h per chain | 1,171 |
| threaded, 32 workers: 36,194 × 7.65 s × 32 CPUs; wall 19.2 h per chain | 2,463 |

The "both adjustments" row is this note's arithmetic, not a row of `projection.md`.

**Verdict: PASS against the 2026-11-26 threshold.** The threshold is at most 5,000 CPU-h for the
gene-expression block over 455 genes by the data plug-in (Decision log 2026-11-26, closing
KIM-176). The headline is 539 CPU-h, about 9× under it. Every row except naive AD is also under
5,000 CPU-h, including both adjustments together (about 1,790). Serial, one chain is 9,048
gradients × 53.61 s = 135 h, inside milan's 7-day limit. With 16 workers it is about 18 h.

### What the projection assumes

1. **The chain design is KIM-116's plug-in stage 1:**
   - 4 chains;
   - NUTS(200 adaptation, 0.8), max depth 10;
   - 200 kept draws per chain, plus the one extension of 200;
   - the warm diagonal metric.
2. **The same leapfrog steps per iteration at 456 coordinates as at 18.** This was not measured.
   - The (D/18)^¼ row is the usual scaling of the step size with dimension.
   - A posterior with stronger correlations or worse conditioning than the reduced model's would
     need more steps. Nothing here measures that.
3. **The cost per gradient is the serial dual pass at the truth,** on one EPYC 7543 core.
   - The caps are fixed by the data and the support. Elsewhere in the posterior only the term
     counts move (through Λ = k + μ·cap).
   - Sampler overhead is not in the headline. The ratio row carries it.
4. **The AD is per gene.** Each gene's term depends only on its own promoter and on krnadeg, so
   a dual pass with two partials per gene is exact.
   - Neither KIM-116's sampler nor src's `plugin_loglik` does this today: they run
     `ForwardDiff.gradient` over all coordinates.
   - At 456 coordinates that is 38 chunk passes and 74,158 CPU-h, which fails the threshold.
   - The verdict therefore depends on wiring a per-gene gradient into the sampler.
5. **Stage 2 is excluded.** That is the bulk/metabolic block, at the frozen threshold's scope.
6. **The model is the reduced model's gene block applied to 455 genes:**
   - one gene copy, no replication over the window;
   - death rate krnadeg / len (the full model's own constant is 0.821× that at the truth);
   - first two bases AU;
   - pools fixed from the panel, as the plug-in fixes them.
7. **The data design is the reduced model's:** 200 cells, 52 windows, every transcript observed
   exactly. A design with more cells or windows scales the cost about linearly.

### What this does not say

- **Affordable is not right.** A pass here does not make the route right; that is gate 3
  (project description). KIM-116 found the plug-in biased on the reduced model: promoter 0607 sat
  +0.65 joint SD off, and 16 of 17 promoters shifted up. This timing does not touch that.
- **Threading.** The threaded figures are as measured. GC takes about a third of the time at
  16 workers, and nothing was tuned.
