# Paper 1 draft: status

**Build:** `latexmk -pdf main.tex` in `paper/`, on OzStar with TinyTeX on `PATH`
(`~/.TinyTeX/bin/x86_64-linux`, TeX Live 2026). The system TeX is too old. **Format:** two-column
preprint with `arxiv.sty`, as in the companion paper; this replaced the PLOS template on 2026-11-07, and
the PLOS preamble and Author Summary are in git history. **Venue:** open.
**Story (2026-11-23, KIM-98):** `paper/story.md`, from the Linear project description; supersedes the 2026-11-07 story below. **Title:** *InferCell: a compositional inference framework for whole-cell modelling, demonstrated on
reduced JCVI-syn3A*. **Tracking:** Linear, milestone "M6 Paper 1 complete".
**Review:** the pre-submission review is kept outside the repository, in
`/fred/oz022/tkimpson/infercell-m6-review-2026-11-07/`. It predates the restructure and assumed PLOS.
**Dependencies:** Python with matplotlib 3.7 or later, scipy and numpy; TinyTeX with the packages in
`main.tex`'s preamble (`authblk` comes from TeX Live's `preprint` package).

**Numbers:** every quoted number is a `\pn{key}` in `numbers.tex`, with its source.
**Figures:** `figures/make_figures.py` (data, from the KIM-43 archive and `figures/data/`) and `figures/make_schematics.py`. `figures/budget_dump.jl` writes `figures/data/budget_mu.tsv` from the archive. `fig_trajectories` reads
`figures/data/traj_*.tsv`, copied from `/fred/oz022/tkimpson/infercell-trajectories-2026-11-08/tsv/`
(KIM-87, `dev/scripts/trajectories_predict_result.md` on main); that archive is not yet in the KIM-54 deposit.
**Tables:** the model-form counts are generated: `tables/reductions_dump.jl` (Slurm job 20165612,
package source identical to main at `3c677a1`), then `tables/make_modelform.py`. Table 1's evidence
is `tables/frameworks_sources.md`.

## Open from the KIM-98 draft (2026-11-23)

The paper was restructured around the reduced model (`paper/story.md`). Decisions for the author:
- **Contribution sentence without M10:** it reads as one sentence only with the ladder clause; if paper 1 ships before M10, reword it (release-vs-M10 question, parked).
- **tab:frameworks:** InferCell's "Inference on dynamics" cell changed from "On M0" to "Both blocks (h)". Confirm.
- **Name of Q:** defined in §2 as extracellular lactate = lactate a cell exports over the 6,300 s cycle, in particles, population mean (`:lactate_out`). "Lactate exported" would be the more literal name; changing it touches abstract, intro, §2, §6, §8.
- **Discussion vs sec:info:** the Discussion draws a working rule from `sec:info`, whose opening says no main-text claim rests on it. One should give way.
- **§2 code listing** illustrates releasing constants (an appendix use) and shows `COREA_CYCLE_S`; consider moving it to app_m0predict.
- **T1:** the M0 mixed-route wall-clock is `\DRAFT{?}` (no source found). The 492 CPU-h counts the whole run including the 24 CPU-h probe.
- **Precision fixes** recorded in story.md's amendment log (gap parameters S_0779/krnadeg/ENO; pools 0.9–1.4 SD; skill 0.98–0.99; "θ-driven" only for trioses and FDP; "about 80×" is 77–103×).
- **Citations:** `[CITE NEEDED: CRPS, Gneiting and Raftery 2007]` in app_methods; §6's CRPS sentence also uncited.

## Open before submission

| Item | Owner |
|---|---|
| `\placeholder`: thanks. Check the affiliation, ORCID and email, filled from the companion paper. Funding and contributions depend on the venue | the author; KIM-49 |
| `\placeholder`: the Zenodo DOI, and what goes in the deposit | KIM-54 |
| The scope risk and likely referee questions in the review | the author |
| Spot-check Table 1, read the supplements | the author |
| Co-authors, if any: sole author for now, by the author's decision of 2026-11-06 | the author |
| Convert to the chosen venue's format | at submission |

## Standing flags

- **AD and "one ODE":** only route (ii) and all-ODE compositions are one `ODEProblem`. The hybrid is not
  differentiable. Keep the scoping in the InferCell section.
- **The story (the author, 2026-11-07, replacing the 2026-11-01 "lead with the reduced model"):** a test
  system, then posteriors (M0 only), then predictions (scored on M0, then unscored on the reduced model).
  The transfer is assumed, not tested, and a test is likely paper 2.
- **On the exact sampler the forward direction is unmeasured,** so do not claim ptsG shrinks under
  metabolites only on that route. Option (b), about 550 CPU-h, would be needed.
- **Word rule:** "detectable" or "declared" discrepancies, not "resolved". "The reduced syn3A model", not
  "Core A′", in user-facing text.
- **B8 is withdrawn** with the old outline's §5.4 and F6, parked outside this paper (KIM-17). B7's matched KS null is
  sound (KIM-67).
- **M0 and the mixed route:** 17 of 50 mixed replicates and 23 of 50 exact replicates fall short of the
  convergence targets. The text says so.
- **Companion paper:** syn3a-ensemble is a `@misc` "in preparation"; update it when it is public.
