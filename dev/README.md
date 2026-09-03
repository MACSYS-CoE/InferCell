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
- `plans/` — how the work is sequenced, not what it is:
  [`reduced-syn3a-wave-plan.md`](plans/reduced-syn3a-wave-plan.md) holds the
  Core A' dependency graph. **Superseded** by `../spec/spec.md`, which reorders
  it; kept for the reasoning behind the decomposition.
- `talks/` — slide decks and talk material.
- `council/`, `references/` — local-only (git-ignored): council transcripts and
  large reference PDFs.

The authoritative spec is [`../spec/spec.md`](../spec/spec.md); session state lives in
[`../docs/handoff.md`](../docs/handoff.md). `plans/` was cleared on 2026-08-26
and now holds only cross-change sequencing. The old `plans/` and `log/` contents
are in git history (`git log --diff-filter=D --stat`).

- `archive/openspec/` — the retired OpenSpec root, read-only. Its
  `specs/corea-interface/` holds the frozen wave-0 interface requirements the
  code implements, and `changes/` holds four drafted module designs whose derived
  findings are absorbed into §4 of the spec. Nothing here is live.
