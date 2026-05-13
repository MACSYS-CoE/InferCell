# Handoff — 2026-05-12

## Goal

Scope a v0.0.1 minimal WCM (issue #6), settle the project's positioning, and tee up implementation of Steps 4–5 to be picked up on the HPC cluster.

## Status

Scoping and positioning complete and committed on branch `minimal-cell-scoping`. No code changes this session — all work is docs. Ready to start implementation (Step 4 or 5) on the cluster.

## What changed this session

- Consolidated `docs/plans/2026-04-13-v1-minimal-demo.md` into `docs/plans/2026-05-12-v0.0.1-scoping.md` — single canonical v0.0.1 plan that pairs scoping (syn3A, closed-loop three-block architecture, bursty regulator) with the 6-step engineering roadmap and the two-tier synthetic data strategy. Old file deleted; refs in `docs/handoff.md` and `docs/grant-application/questions_and_answers.md` repointed; `.gitignore` allowlist updated. Commit: `2c0a896`.
- Added `## Capabilities` section to `README.md` and wrote `docs/positioning.md` — frames inference as the *design discipline*, with parameter inference, sensitivity analysis, model reduction, and model comparison as downstream capabilities of that discipline (not separate features). Commit: `79f72b8`.
- Sharpened the Tier A vs Tier B framing in the v0.0.1 plan: Tier A (synthetic twin) shows the boundary protocol works as a *computational technique* (it has to — NUTS-ODE and ABC-SMC-SSA can't be jointly composed). Tier B (refit under mismatch) is what shows it's *honest under model error* — the actual scientific claim.

## State of the v0.0.1 plan

Steps 1–3 done on `main` (TX/TL ODE + NUTS, Stochastic GE + ABC-SMC, LightMetabolism + shared params, joint NUTS, Level-1 boundary protocol via KDE priors, 121+ tests). The new plan has 6 steps total. Next:

- **Step 4 — Level 2 boundary protocol.** Iterative ODE↔SSA message-passing with a KL-based convergence criterion. Demonstrate convergence and improved parameter recovery vs Step 3's single-pass conditioning.
- **Step 5 — Close the v0.0.1 biology.** (1) Enzyme→metabolism feedback: block 2's enzyme protein concentration enters block 1 as a rate parameter, closing the autocatalytic loop. (2) Replace constitutive `StochasticGE` with a two-state bursty regulator. (3) Assign specific syn3A genes to block 2 slots. (4) Build the Tier B "v0.0.1+ε" mismatch generator.

See `docs/plans/2026-05-12-v0.0.1-scoping.md` for full detail, plus Step 6 (end-to-end demo + figure) and the risks/AD-backend notes.

## Open questions (carried from the plan)

- Specific syn3A genes for block 2 — at minimum need ≥1 metabolic enzyme (its product is in block 1) and ≥1 ribosomal component (gates TL capacity).
- Which syn3A gene plays the bursty regulator role (sigma-factor analogue is the obvious candidate).
- Tier B mismatch design: probably one experiment each for (i) extra reactions in block 1, (ii) stochasticity in block 2, (iii) finite-volume effects on block 3.
- Should `LightMetabolism` evolve in-place, or fork to a `Syn3AMetabolism` with named species and reactions?
- Recommended start order (from this session): Step 5 first (mostly local code, light compute), then Step 4 (algorithmic, heavier compute) — but Tom hasn't committed to an order.

## Blockers / problems

None. User is moving to the cluster to start implementation.

## Next steps

1. On the cluster, pull to current HEAD (`79f72b8` on `minimal-cell-scoping`).
2. Decide Step 4 vs Step 5 start order. Suggested: Step 5 first.
3. For Step 5: wire enzyme→metabolism coupling (watch for ODE instability — start at low gain, Michaelis–Menten fallback in the risks table); two-state bursty SSA replacing `StochasticGeneExpression`; commit specific syn3A gene assignments; build the Tier B generator.
4. For Step 4: iterative loop, KL convergence criterion, integration test showing tightening vs single-pass.
5. Plan for ForwardDiff → Mooncake benchmark as parameter count climbs past ~10 (AD-backend note in the plan).

## Non-obvious context

- Branch `minimal-cell-scoping` is a docs-only scoping branch — every commit on it is planning, no code. Implementation work should be done on a fresh branch off `main` (or the cluster's equivalent), not extended off this branch unless that's the explicit intent.
- `.gitignore` uses a strict per-file allowlist for `docs/` (`docs/**` is ignored except for explicit `!docs/path/to/file.md` exceptions). Any new doc file needs to be added to that allowlist or it silently won't be tracked. The grant Q&A under `docs/grant-application/` is intentionally not allowlisted — local-only, edits there don't show in git.
- `docs/handoff.md` (this file) is the project's handoff convention, separate from the `handoff.md`-at-root pattern that the `/handoff` skill suggests. Keep updating this file rather than creating a root-level one.
- The Capabilities section in the README and the longer `docs/positioning.md` are the canonical "how to describe this project" — use them for grant text, talks, slides, etc., rather than reinventing the framing each time.
