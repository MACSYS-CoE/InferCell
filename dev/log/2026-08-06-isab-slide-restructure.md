# Session Log: ISAB Slides 10–12 Restructured Around the Two-Sampler Contrast

**Date:** 2026-08-06 (second session of the day)
**Branch:** `minimal-cell-example`
**Status:** Slides and plotting code done. **All three corner figures are stale**
— they must be regenerated on OzStar before the deck is presentable.

## Why this session happened

The previous session left slide 10 (the minimal-cell schematic) paired with a
single results slide. Two problems surfaced while reading it back: the schematic
sat on the right of the frame with the dynamics sketches on the left, which is
backwards for a slide whose subject is the cell; and the results slide asked the
audience to absorb a two-colour overlay without ever having seen what one colour
alone looked like.

## What the flow became

Worked out conversationally before touching the file. Six beats:

1. here is a minimal cell — modules, parameters, deliberately simple
2. we want to infer those parameters, conditioned on data
3. some blocks are ODEs (NUTS), some are stochastic (SBI)
4. and some parameters are shared across that divide
5. one way: run each block with its own sampler, independently
6. another way: let the posterior of one become the prior of the other

Beats 1–4 turned out to be *already drawn* on slide 10 — the boxes and counts,
the observed-data node, the two dynamics sketches with their sampler labels, and
the red shared-parameter annotations. So beats 1–4 became spoken text over an
unchanged figure, and only beats 5–6 needed new slides.

Two things were flagged as needing to be said out loud rather than drawn:

- Which of the two "3 shared" labels is the actual wrinkle. The intra-block
  three are absorbed by the joint NUTS fit and are a non-problem; only the three
  crossing the boundary matter. The schematic shows both.
- The missing half-beat between 4 and 5: *why not one sampler over all 11?*
  No gradient path through the SSA kills NUTS for the joint, and ABC over 11
  parameters with an ODE inside is wasteful. That is what forces block-wise
  inference in the first place.

## The two-corner argument, which was initially resisted and shouldn't have been

The first build of slide 11 used a single blue-only corner, on the reasoning
that an 8-parameter corner beside a 6-parameter one would be unreadable at slide
size. That reasoning was wrong, and the user overrode it correctly.

`ODE_SUBSET` in `plot_bursty_boundary.py` already reduces the differentiable
corner to a 4×4 — precisely because the full 8×8 is illegible on a slide. And
those four are `k_tx` plus exactly the three parameters shared across the
boundary. So the two panels not only fit, they carry the argument: the three
tight columns on the left *are* the three broad ones inside the dashed box on
the right. Point at one, then the other.

The lesson (saved to memory): open the figure before arguing from the parameter
count in the model.

## Consequence: the build became a zoom, not an overlay

The original plan had slides 11 and 12 pixel-aligned so that advancing added the
green contours and moved nothing else — the reason `--mode alone` computes its
axis ranges from *both* posteriors and reserves an invisible legend slot for the
conditioned series. Once slide 11 had to hold two corners, the SSA panel shrank
and the alignment was lost. Slide 12 is now a zoom-in rather than an overlay
reveal. The range-sharing and legend-slot machinery was kept anyway: it costs
nothing and the two panels still have to agree panel-for-panel.

## The palette collision

Putting the two corners side by side exposed an inconsistency that had been
survivable while they lived on separate slides. Both were drawn in the same
`COL_UNCOND` blue, while the slide-10 schematic codes blue = differentiable and
orange = stochastic — so the regulator was orange on slide 10 and blue on 11,
and the two adjacent corners were the same colour for different quantities.

Options weighed: recolour only the ODE corner to something neutral; leave it and
disambiguate with captions; or align the plots to the schematic. Chose the last.
Truth moved off orange (which the stochastic block now owns) to black, which is
the corner-plot convention anyway and stops truth lines competing with a
posterior for the eye.

```
COL_DIFF  #3B7EA1  blue    differentiable block's own posterior
COL_ALONE #D2762B  orange  regulator inferred alone
COL_COND  #3D9970  green   conditioned across the boundary
COL_TRUTH #1A1A1A  black   ground truth
```

The chips on slides 11–12 are `\definecolor` hex matches to these constants,
with cross-referencing comments in both files so they can't drift apart quietly.

## What could not be verified

Nothing in the plotting script was executed. `examples/` is gitignored, so the
particle CSVs from job 15152954 never came back from the cluster, and `corner`
is not installed in the laptop's default `python3`. The palette change and the
`--mode alone` code path are verified by reading and an AST parse only.

Slide 11's missing panel renders a deliberately loud red **FIGURE PENDING** box
via `\IfFileExists`, chosen over silently reusing the overlay figure — the
latter would spoil the reveal in a way that is easy to miss during rehearsal.

## Artefacts

- `dev/talks/ISAB/talk.tex` — slides 10, 11 (new), 12; figure-series colours
  defined at the top; backup slide 16 reframed as a full-size repeat
- `dev/talks/ISAB/figures/plot_bursty_boundary.py` — `--mode alone`, palette
  constants renamed and recoloured, `plot_ssa()` extracted from `main()`
