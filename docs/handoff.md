# Handoff — 2026-08-06 (evening)

## Goal

Restructure ISAB slides 10–12 so the deck earns its result: schematic first,
then the two blocks inferred independently, then the same regulator posterior
conditioned across the boundary.

## Status

Slides and plotting code are **done**; the deck compiles clean at 20 pages.
But **all three corner figures on disk are stale** and the deck is not
presentable until they are regenerated on OzStar. The talk is 12 Aug 2026.

## What changed this session

See `dev/log/2026-08-06-isab-slide-restructure.md`.

## Blockers / problems

**The only blocker: three figures need regenerating, and the inputs are on the
cluster.** `examples/` is gitignored, so the particle CSVs from job 15152954
never left OzStar. The Aug 2026 palette change (blue = differentiable, orange =
regulator-alone, green = conditioned, black = truth) made every existing corner
PNG stale, not just the new one.

On OzStar, either re-run `sbatch test/run_talk_bursty_boundary.slurm` or locate
the job's `examples/results/isab_*.csv`, then from the repo root:

```
python3 dev/talks/ISAB/figures/plot_bursty_boundary.py --mode ode
python3 dev/talks/ISAB/figures/plot_bursty_boundary.py --mode alone
python3 dev/talks/ISAB/figures/plot_bursty_boundary.py --mode ssa
```

`alone` and `ssa` **must** come from the same particles or slides 11 and 12
disagree. Needs `corner` installed (the laptop's default `python3` lacks it).

Until `isab_ssa_alone.png` exists, slide 11's right panel renders a red
**FIGURE PENDING** box — deliberate, so it cannot be missed in rehearsal.

Unchanged from the previous session: the promoter rates are not individually
identifiable from replicate-mean summary statistics; only the combination
`k_tx_burst·k_on/(k_on+k_off)` is, and it is recovered. Disclosed on backup
slide 16. The fix is distributional summaries in `src/summary_statistics.jl`.

## Open questions

- **`CLAUDE.md` does not exist at the repo root.** Worth adding one that points
  future sessions at `docs/handoff.md` and `dev/log/` — this repo's convention
  differs from the `/handoff` skill's default (`./handoff.md`, `./log.md`), so
  every session has to rediscover it. Not created without a say-so.
- **PR #37** was open against `main` at the last handoff. Still needs review and
  merge; this session's commits land on the same branch.
- Carried over: commit `dev/talks/ISAB/talk.pdf`? Use vector `.pdf` figures?

## Next steps

1. On OzStar: regenerate the three corner figures (command block above),
   rebuild the deck, and eyeball slides 11–12 — this is the one thing standing
   between the deck and presentable.
2. Check the new palette actually reads on a projector, especially orange
   (`#D2762B`) against green (`#3D9970`) for anyone colour-blind. If it does
   not, the four constants are in one place at the top of the plot script.
3. Rehearse slides 10 → 11 → 12 as a unit. Slide 10 now ends on a question
   ("How do you infer θ when half the cell has no gradients?") that slides 11
   and 12 are the two answers to — say it out loud before advancing.
4. Review and merge PR #37.

## Non-obvious context

- **Say out loud on slide 12, deliberately not on the slide:** the flow is
  one-way (differentiable → regulator); nothing the regulator's data knows
  travels back. This is a cut posterior, a principled approximation to the
  joint rather than the joint itself. The iterative protocol covered truth on
  35–70% of seeds at a nominal 95%, which is why single-pass is the default.
  There is a `SAY OUT LOUD` comment on the frame.
- **Slide 10 shows two different "3 shared" labels.** Only the boundary-crossing
  three are the problem; the intra-block three are absorbed by the joint NUTS
  fit. Be explicit when speaking it or they read as one thing.
- **The two corner panels on slide 11 are the argument, not decoration.**
  `ODE_SUBSET` reduces the left panel to `k_tx` plus exactly the three shared
  parameters, so its three tight columns are the three broad ones inside the
  dashed box on the right.
- **Colour chips in `talk.tex` are hex copies of the script's constants**
  (`figDiff`/`figAlone`/`figCond`). Change one, change both.
- Parameter order in the boundary figure is still load-bearing — unshared three
  first, shared three last. The script warns and drops the dashed box rather
  than drawing a wrong one; check the log after any rerun.
- Building the deck on the cluster needs `~/.TinyTeX` on PATH; the system TeX
  has no `beamer`.
- Do not trust smoke-setting runs: at 100 NUTS samples the differentiable chain
  does not converge and the boundary passes the overconfidence downstream,
  making the result look better than it is.
