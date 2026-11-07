# Paper 1 draft: status

**Branch:** `paper-1-draft`, PR #100. **Build:** `latexmk -pdf main.tex` in `paper/`, on OzStar with
TinyTeX on `PATH` (`~/.TinyTeX/bin/x86_64-linux`, TeX Live 2026). The system TeX is too old for the
PLOS template. **Venue:** PLOS Computational Biology, Methods, in PLOS's template (v3.8).
**Title:** *InferCell: a compositional inference framework for whole-cell modelling, demonstrated on
reduced JCVI-syn3A*. **Tracking:** Linear, milestone "M6 Paper 1 complete".
**Review:** `PRE_SUBMISSION_REVIEW_2026-11-07.md` holds what was fixed and what is left.

**Numbers:** every quoted number is a `\pn{key}` in `numbers.tex`, with its source.
**Figures:** `figures/make_figures.py` (data, from the KIM-43 archive) and `figures/make_schematics.py`.
**Tables:** the model-form counts are generated: `tables/reductions_dump.jl` (Slurm job 20165612,
package source identical to main at `3c677a1`), then `tables/make_modelform.py`. Table 1's evidence
is `tables/frameworks_sources.md`.

## Open before submission

| Item | Owner |
|---|---|
| `\placeholder`: affiliation, email, thanks. Funding, ORCID and contributions go in Editorial Manager | the author; KIM-49 |
| `\placeholder`: the Zenodo DOI, and what goes in the deposit | KIM-54 |
| The scope risk and likely referee questions in the review | the author |
| Spot-check Table 1, read the supplements | the author |
| Co-authors, if any: sole author for now, by the author's decision of 2026-11-06 | the author |
| Flatten `\input`, upload figures separately as TIFF or EPS | at submission |

## Standing flags

- **AD and "one ODE":** only route (ii) and all-ODE compositions are one `ODEProblem`. The hybrid is not
  differentiable. Keep §2's scoping.
- **The framing (the author, 2026-11-01):** lead with the reduced model, then score on M0 as the
  largest checkable slice. The transfer is assumed, not tested, and a test is likely paper 2.
- **On the exact sampler the forward direction is unmeasured,** so do not claim ptsG shrinks under
  metabolites only on that route. Option (b), about 550 CPU-h, would be needed.
- **Word rule:** "detectable" or "declared" discrepancies, not "resolved". "The reduced syn3A model", not
  "Core A′", in user-facing text.
- **B8 is withdrawn** with §5.4 and F6, parked outside this paper (KIM-17). B7's matched KS null is
  sound (KIM-67).
- **M0 and the mixed route:** 17 of 50 mixed replicates and 23 of 50 exact replicates fall short of the
  convergence targets. The text says so.
- **Companion paper:** syn3a-ensemble is a `@misc` "in preparation"; update it when it is public.
