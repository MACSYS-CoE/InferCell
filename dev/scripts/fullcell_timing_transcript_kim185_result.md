# KIM-185: the transcript factor on 455 genes

Full-cell timing 3 (Linear M11). The gene-expression block of the full syn3A model, by the data
plug-in: 455 genes, the reduced model's design, src's exact transcript factor.

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
