# Handoff — 2026-08-07 (morning)

## Goal

Regenerate the ISAB corner figures against the August 2026 palette so the deck
is presentable. Talk is 12 Aug 2026.

## Status

**Done — the figure blocker is cleared.** All four figures regenerated from the
existing particle CSVs, `isab_ssa_alone.png` now exists, and the deck compiles
clean at 20 pages with no `FIGURE PENDING` box remaining.

The previous handoff claimed the particle CSVs "never left OzStar" and that the
job had to be re-run. That was wrong from the wrong machine: this repo lives on
OzStar at `/fred/oz022/tkimpson/InferCell`, `examples/results/isab_*.csv` (job
15152954, 6 Aug 15:23) were sitting there the whole time, and `corner` 2.2.3
imports in the default `python3`. Nothing needed re-running on the cluster.

## What changed this session

Regenerated, from the repo root, with no `--results` override:

```
python3 dev/talks/ISAB/figures/plot_bursty_boundary.py --mode ode
python3 dev/talks/ISAB/figures/plot_bursty_boundary.py --mode alone
python3 dev/talks/ISAB/figures/plot_bursty_boundary.py --mode ssa
python3 dev/talks/ISAB/figures/plot_bursty_boundary.py --mode burst
```

Note the fourth: `--mode burst` also draws from the recoloured constants
(`COL_ALONE`/`COL_COND`/`COL_TRUTH`), so backup slide 17 was stale too. The
previous handoff's command block listed only three.

Then two reduced figures, for when the 6x6 corner is too dense to read from the
back of a room:

```
python3 dev/talks/ISAB/figures/plot_bursty_boundary.py --mode alone --shared-only
python3 dev/talks/ISAB/figures/plot_bursty_boundary.py --mode ssa   --shared-only
```

`--shared-only` keeps just `k_tl`, `γ_mRNA`, `γ_prot` and writes
`isab_ssa_alone_shared` / `isab_ssa_posteriors_shared`. **Neither is referenced
by `talk.tex` yet** — they exist to be dropped in if slides 11–12 read badly in
rehearsal.

Verified rather than assumed:

- No parameter-order warning fired, and the dashed "shared with the
  differentiable block" box is drawn on both corners around the last three
  columns (`k_tl`, `γ_mRNA`, `γ_prot`). Ordering is intact.
- Slides 11 and 12 rendered to PNG and inspected. Palette matches slide 10:
  blue = differentiable, orange = regulator-alone, green = conditioned,
  black = truth.
- `pdftotext | grep -c "FIGURE PENDING"` → 0.

Committed: the four PNGs only.

## Open questions

- **Do slides 11–12 want the reduced (`_shared`) corners instead?** They are
  generated and committed but unwired. The full 6x6 carries more of the story —
  it shows the three promoter rates staying broad, which is the honest caveat —
  while the 3x3 is what actually reads at slide size. Possible compromise: 3x3
  on 11–12, full 6x6 promoted to a backup slide. Decide in rehearsal, in a room.
- **Not committed, deliberately: the `.pdf` twins** of all four figures and
  `dev/talks/ISAB/talk.pdf`. The plot script emits both formats; `talk.tex`
  still `\includegraphics` the `.png`s. Switching the deck to vector figures is
  a real improvement and a one-line-per-frame change, but it was not this
  session's ask. Decide before the talk, not during.
- **`.github/workflows/CI.yml` is untracked** and unrelated to the talk. Left
  alone. Someone should decide whether it belongs on this branch or its own.
- **PR #37** against `main` still needs review and merge.

## Next steps

1. **Look at the deck on a projector.** Specifically orange (`#D2762B`) against
   green (`#3D9970`) for anyone colour-blind — magenta truth (`#E5007D`) was
   chosen partly to stay clear of that pair, but the orange/green pair itself is
   still unverified in a room. This cannot be checked from a terminal. The five
   constants are together at the top of `plot_bursty_boundary.py`, and
   `talk.tex` carries hex copies as `figDiff`/`figAlone`/`figCond`/`figTruth`.
2. Rehearse slides 10 → 11 → 12 as a unit. Slide 10 ends on a question ("How do
   you infer θ when half the cell has no gradients?") that 11 and 12 are the two
   answers to — say it out loud before advancing.
3. Decide the `.pdf`-vs-`.png` question above.
4. Review and merge PR #37.

## Non-obvious context

- **Truth is magenta (`#E5007D`), not black, as of this session** — black was
  disappearing into the dark contour cores at slide size. Cyan was tried first
  and rejected: on the *blue* corner it reads as a tint of `COL_DIFF` rather
  than a separate quantity, which is the one panel where that matters. Magenta
  is the only hue no posterior owns, and it keeps its blue channel so it stays
  separable from the orange/green pair under deuteranopia. `thicken_truths()`
  also widens the lines to 2.6 — corner takes `truth_color` but no linewidth,
  so the artists have to be reached after the fact.
- **Say out loud on slide 12, deliberately not on the slide:** the flow is
  one-way (differentiable → regulator); nothing the regulator's data knows
  travels back. This is a cut posterior, a principled approximation to the joint
  rather than the joint itself. The iterative protocol covered truth on 35–70%
  of seeds at a nominal 95%, which is why single-pass is the default. There is a
  `SAY OUT LOUD` comment on the frame.
- **Slide 10 shows two different "3 shared" labels.** Only the boundary-crossing
  three are the problem; the intra-block three are absorbed by the joint NUTS
  fit. Be explicit when speaking it or they read as one thing.
- **The two corner panels on slide 11 are the argument, not decoration.**
  `ODE_SUBSET` reduces the left panel to `k_tx` plus exactly the three shared
  parameters, so its three tight columns are the three broad ones inside the
  dashed box on the right.
- **Colour chips in `talk.tex` are hex copies of the script's constants**
  (`figDiff`/`figAlone`/`figCond`/`figTruth`). Change one, change both. Backup
  slide 16 also *names* the truth colour in prose — it read "Orange lines mark
  truth" for two palette revisions before anyone noticed, so grep the deck for
  colour words, not just for `\definecolor`, after any recolour.
- Carried over: the promoter rates are not individually identifiable from
  replicate-mean summary statistics; only `k_tx_burst·k_on/(k_on+k_off)` is, and
  it is recovered. Disclosed on backup slide 16. The fix is distributional
  summaries in `src/summary_statistics.jl`.
- Building the deck needs `~/.TinyTeX/bin/x86_64-linux` on PATH; the system TeX
  has no `beamer`.
- Do not trust smoke-setting runs: at 100 NUTS samples the differentiable chain
  does not converge and the boundary passes the overconfidence downstream,
  making the result look better than it is.
