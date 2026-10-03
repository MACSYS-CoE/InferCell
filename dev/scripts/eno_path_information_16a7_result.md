# F6's ENO cell: refuted by an upper bound

Jobs 17938700 to 17938702 at `2b9dba9`, `dev/scripts/eno_path_information_16a7.jl`.
40 Core A′ cells at 15.7's truth (seeds 30001 to 30040), each a full 6,300 s path
replayed at ln θ ± δ for ENO and for FBA. About 5.1 CPU-h.

**Why it is a bound.** The transcripts are a function of a cell's full jump path. So
the Fisher information they carry about a metabolic constant cannot exceed the
path's own (the data-processing inequality). The path density sees ENO and FBA only
through the rebuilt rate constants: transcription from ATP and GTP, translation
from charged tRNA.

**How little the path moves.** Changing a constant by 20% (δ = 0.2) changes one
cell's whole-path log-density by:

| Constant | Mean \|Δ log p\| per cell | Largest of 40 |
|---|---|---|
| ENO, δ = 0.2 | 0.0025 to 0.0030 nats | 0.013 |
| FBA, δ = 0.2 | 0.0015 to 0.0021 nats | 0.010 |

**Scaled to 200 cells:** posterior over prior SD from the complete-data Fisher
information, under each target's lognormal prior (ENO SD 0.157 and FBA SD 0.558 in
ln θ).

| δ | Score mean (ENO, FBA) | ENO ratio, expected information | FBA ratio | ENO, 95% bootstrap | ENO ratio, observed | FBA ratio, observed |
|---|---|---|---|---|---|---|
| 0.05 | −0.001, +0.001 | 0.999 | 0.995 | 0.998 to 1.000 | 0.987 | 0.608 |
| 0.2 | −0.002, +0.002 | 0.999 | 0.995 | 0.999 to 1.000 | 0.996 | 0.842 |

- **ENO cannot shrink from transcripts.** Even the full path, every event and its
  time, would leave ENO's posterior SD at 99.9% of its prior over 200 cells. The
  transcripts carry no more. F6's second load-bearing cell, ENO shrinking under
  transcripts only, is **refuted at 15.7's truth**, whatever the observation model
  or sampler.
- **The scores average to zero,** as they should at the truth.
- **FBA's observed-information ratio disagrees with its expected one** (0.61 and
  0.84 against 0.995). The curvature of differences of order 1e-3 nats is noisy, and
  the path density may not be quadratic across a clip (16a.6). The expected
  information is the robust figure. ENO's two agree.
- **Scope.** This is one truth, and pools in another regime could be more
  sensitive. It is linearised at the truth, from 40 cells scaled to 200, and on the
  published (clamped) model.

**What it means.** Information crosses the seam in one direction only. Metabolite
data inform the gene-expression parameters (the bulk check: ptsG shrinks to about
half its noise bracket). Gene-expression data say essentially nothing about the
metabolic constants. The rate constants that read the pools are rebuilt once a
minute, so they damp the effect.
