# Handoff — 2026-08-05

## Goal

Make the Step 3 corner-plot figures in the talks legible on a projector. The
posteriors were fine; the fonts were sized for on-screen viewing, so at the size
the figures are actually shown (6.2cm and 3.35cm tall) the tick numbers and axis
labels were unreadable.

## Status

Done, on branch `talk-figure-legible-fonts` (off `main`). Figures regenerated and
deployed to `docs/talks/figures/`. No inference was re-run — the plotting script
is deterministic and reads cached results from `examples/results/`.

## The main finding: the generator was untracked

The script that produces these figures is
**`examples/step3_plot_results.py`** — **Python + `corner` 2.2.3 + matplotlib**,
living only on the OzSTAR cluster at `/fred/oz022/tkimpson/InferCell`.
`.gitignore` had a bare `examples/`, so it appeared in no branch of the GitHub
remote and the talk figures were unreproducible from a fresh clone. A prior
session lost real time searching for it.

Two misconceptions worth recording, because both were plausible and both wrong:

- It is **not** PairPlots.jl / CairoMakie. `Project.toml` on
  `origin/issue-11-calibration-sweep` declares those, but nothing uses them for
  these figures.
- `examples/step3_plot_results.jl` exists and is a **stale sibling** — it is not
  what generated the committed PNGs. The `.py` file is.

Identity was verified rather than assumed: running the script untouched
reproduced both PNGs byte-identically (`md5 95fabd86…` for
`step3_ssa_posteriors.png`, `538c3e0f…` for `step3_corner.png`). It is
reproducible because the ABC-SMC resampling uses `default_rng(42)`.

**This branch un-ignores `examples/*.py`, `*.jl`, `*.slurm`** so it cannot recur.
Data, figures and logs (`examples/results/`, `examples/figures/`, `*.log`) stay
ignored.

## What changed this session

**`examples/step3_plot_results.py`**
- Font sizes hoisted to named constants sized for the *on-slide* size, not the
  on-screen size: `FS_LABEL=38`, `FS_TICK=28`, `FS_TITLE=26`, `FS_LEGEND=30`.
  Tick labels went 11pt → 28pt, axis labels 15pt → 38pt. The 4-panel corner
  exports ~9.4in wide and is shown 6.2cm tall, i.e. shrunk ~4x, so anything
  below ~25pt in the figure is illegible on a projector.
- `max_n_ticks=3` (from corner's default 5) — 5 ticks crowd at this scale.
- `labelpad=0.28` so the enlarged axis labels clear the enlarged tick numbers;
  without it `$\gamma_{mRNA}$` overlapped its own tick values.
- Histogram linewidths 1.5 → 3.0, thicker tick marks, legend moved to
  `(0.99, 0.99)` and kept at 30pt so it clears the diagonal panel titles.
- **`TITLE_MEDIAN_ONLY = True`** (new flag) replaces corner's stacked
  `0.97^{+0.04}_{-0.03}` titles with median-only `$k_{tx} = 0.97$` via the new
  `set_median_titles()` helper. Reason: corner renders the sub/superscripts at
  ~0.7x the title size — about 5pt on a 6.2cm slide, never readable — and the
  stacked string is ~2x wider than its panel, so it overlapped neighbouring
  panels at any legible size. Set the flag `False` to restore corner's
  behaviour. The dashed 16/84% quantile lines still convey spread.
- **New export `step3_corner_reduced.{png,pdf}`** — a 4x4 cut of
  `{k_tx, k_tl, gamma_mRNA, k_ntp}`.

**`.gitignore`** — track `examples/` scripts (see above).

**`docs/talks/figures/`** — regenerated `step3_ssa_posteriors.png`
(2820x2852 → 2918x2989) and `step3_corner.png` (5311x5378 → 5432x5474); added
`step3_corner_reduced.png` (2919x3004). PNG and PDF both refreshed in
`examples/figures/` too. Both talks use `\graphicspath{{../}}`, so
`docs/talks/figures/` is a single shared copy — not one per talk.

## Why the reduced corner contains those four parameters

Not an arbitrary subset. From the correlation matrix of `examples/results/ode_chain.csv`:

- `{k_tx, k_tl, gamma_mRNA, k_ntp}` form one tightly degenerate block, pairwise
  |r| = 0.68–0.88 (strongest: `k_tx`–`k_ntp` at +0.88, `k_tx`–`gamma_mRNA` at
  +0.81, `k_tx`–`k_tl` at −0.80).
- `{gamma_prot, k_aa}` are a weaker second pair (+0.68).
- **`sigma_obs` is uncorrelated with everything** (|r| ≤ 0.04) — it contributes
  nothing but panels.

The 4x4 cut carries the degeneracy story at a size where the ellipse tilts read
clearly even at 3.35cm.

## The full 8x8 is texture by design

`step3_corner.png` is shown at 3.35cm tall, where no per-panel label can ever be
legible — verified by downsampling to actual slide pixels. Its fonts were bumped
moderately (labels 30pt, ticks 20pt) so it survives being shown larger, and it
keeps corner's full stacked-quantile titles for full-size viewing. The legible
cut is `step3_corner_reduced.png`.

## Next steps

1. **Decide whether any slide should switch** from `step3_corner.png` to
   `step3_corner_reduced.png`. No `talk.tex` was edited this session — that is a
   content call. Current refs: `2026-04-macsys/talk.tex:1154` (0.85\textwidth)
   and `2026-05-brisbane-macsys/talk.tex:1133` (0.48\textwidth).
2. **Sync to the other checkouts.** The Mac tree uses `dev/talks/figures/` while
   the cluster uses `docs/talks/figures/` — these have diverged. There are also
   byte-identical copies in iCloud at
   `Work/talks/2026-08-isab/figures/`. `rsync` from
   `ozstar:/fred/oz022/tkimpson/InferCell/docs/talks/figures/`.
3. **`gamma_prot` title precision.** In the full 8x8 it reads
   `0.10^{+0.00}_{-0.00}` — pre-existing, from `title_fmt=".2f"` on a tightly
   constrained parameter. Widen to `.3f` if that figure is ever shown large.
4. **Apply the same font treatment to the Step 4 figures** if they are used at
   similar sizes: `examples/step4_plot_results.py` and
   `examples/step4_calibration_curve.py` generate `step4_calibration_curve.png`
   and `step4_corner_talk_shared.png`, which are also in the ISAB talk folder.
   They were not touched this session.
5. **Consider deleting `examples/step3_plot_results.jl`.** It is stale and it is
   actively misleading — now that it is tracked, it will mislead again. Left in
   place because deleting was out of scope.

## Open questions (carried)

- **Specific syn3A loci** for the enzyme / ribosomal-component / reporter slots —
  placeholder names (`:enzyme`, `:ribosome`, `:reporter`) in the multi-gene
  Block 2; commit to specific loci in v0.0.2 alongside the genome annotation.
- **Iterative over-tightening.** At 95% CI, 2/4 iterative params miss truth
  while single-pass covers all 4. Per auto-memory, `iterative_infer` is a
  proof-of-concept of cross-module exchange, not the production protocol — do
  not over-invest in fixing its calibration.
- **Tier B beyond Block 1.** Block-2 stochasticity and Block-3 finite-volume
  mismatches are deferred.

## Non-obvious context

- **Figure fonts must be chosen against the on-slide size, not the export size.**
  The export is ~9.4in wide; the slide shows it at 6.2cm. That is a ~4x
  reduction, and it is the reason "high resolution" did not help — 2820px at 1%
  margin was always plenty. Divide any font size by ~4 to predict legibility.
- **`examples/` scripts are now tracked but `examples/results/` and
  `examples/figures/` are not.** The plotting scripts read CSVs that exist only
  on the cluster. A fresh clone can read the scripts but cannot re-run them
  without regenerating results via the `.slurm` jobs.
- **`docs/handoff.md` is the project's handoff convention.** This file. Keep
  updating it; do not create a root-level `handoff.md` even though some skills
  suggest that path.
- **`.gitignore` uses a strict per-file allowlist for `docs/`.** Any new doc
  needs an explicit `!docs/path/to/file.md` exception. The grant Q&A under
  `docs/grant-application/` is intentionally untracked.
- **`LightMetabolism` is instance-configurable** (`mRNA_source`,
  `enzyme_source` kwargs). Defaults keep the Step-3 composition path identical.
  If you add a third coupling input, mirror this pattern rather than hardcoding
  state names in `inputs(::LightMetabolism)`.
- **`TranscriptionTranslation` is dual-constructor.** The no-arg form returns
  `gene_names=[:_default]` and the single-gene shape; the vector-arg form
  returns the multi-gene shape. Don't merge these — the `[:_default]` special
  case preserves backward compatibility with Steps 1–3.
- **`BurstyGeneExpression` sits alongside `StochasticGeneExpression`, not
  replacing it.** Step 3's regression test asserts constitutive behaviour.
- **Colour semantics across the talk:** blue `#0072B2` = unconditioned /
  base posterior, green `#009E73` = conditioned overlay, orange `#E69F00` =
  truth. Wong palette. Preserved by this change.
