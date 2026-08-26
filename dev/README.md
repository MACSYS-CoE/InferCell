# dev/

Development and research material for InferCell — **not** part of the published
documentation. Everything here is working material and may be dropped or
rewritten before release. The published, user-facing docs live in
[`../docs/`](../docs/) (rendered with MkDocs).

- `notes/` — research notes. The current scope lives here:
  [`reduced-syn3a-scoping.md`](notes/reduced-syn3a-scoping.md) (the Core A'
  Step 1 target) and [`well-stirred-minimal-cell.md`](notes/well-stirred-minimal-cell.md)
  (the Thornburg 2022 model it is carved out of), with
  [`modular-bayesian-inference-heterogeneous-modules.md`](notes/modular-bayesian-inference-heterogeneous-modules.md)
  and [`model-uncertainty-and-selection.md`](notes/model-uncertainty-and-selection.md)
  on the inference architecture. Background: the case for inference, the
  calibration-curve result, the 4D and spatial/RDME notes.
- `talks/` — slide decks and talk material.
- `council/`, `references/` — local-only (git-ignored): council transcripts and
  large reference PDFs.

Forward-looking planning has moved to spec-driven development in
[`../openspec/`](../openspec/); session state lives in
[`../docs/handoff.md`](../docs/handoff.md). The former `plans/` and `log/`
directories were removed on 2026-08-26 — they are in git history if needed
(`git log --diff-filter=D --stat`).
