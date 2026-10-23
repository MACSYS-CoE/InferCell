# Paper 1 draft: status

**Branch:** `paper-1-draft` (from `main` at `37e680c`). **Build:** `latexmk -pdf main.tex`
in `paper/` (system TeX Live; natbib, cleveref, xcolor present).
**Working title:** *InferCell: inference-first whole-cell modelling*.
**Venue:** journal-agnostic; Cell Systems or PLOS Computational Biology. Two-column ICML
2026 template in `[preprint]` mode, copied from `FNO_for_turing/paper/`.
**Numbers:** every quoted number is a `\pn{key}` in `numbers.tex`, with its source file.
**Placeholders:** red `\placeholder{...}` text and boxed `\placeholderfig` figures.
`grep -n placeholder sections/*.tex` lists them all.

## What is drafted

| Section | State | Sources |
|---|---|---|
| Abstract | Drafted; 3 placeholders (B7, B8, budget) | — |
| 1 Introduction | Drafted; Table 1 cells all placeholders | `dev/notes/the-case-for-inference.md`, outline §1 |
| 2 Design | Drafted; Fig 1 placeholder; B5 and B6 placeholders | `src/`, docs, agent survey |
| 3 Model | Drafted; Fig 2 is `fig1c_state_graph_corea.pdf` | spec §0, §2; `src/organisms/coreA/m0.jl` |
| 4 The boundary | Drafted; Fig 3 placeholder | discussion note §2–§6 |
| 5.1 Exact sampler | Drafted from V6, V7, F10 | `dev/scripts/v6_*`, `v7_*`, `f10_*` result files |
| 5.2 Jensen gap | Drafted from B3, the rethink, the bias report, route (ii) | `jensen_*`, `route_ii_89_result.md` |
| 5.3 Approx vs exact | Placeholder (B7) | — |
| 5.4 Full scale | Placeholder (B8, needs O1) | — |
| 5.5 Uncertainty budget | Placeholder (not built) | outline §6 |
| 5.6 Cost | Drafted; approximate row placeholder | 16-recovery §12 |
| 6 Discussion | Drafted; future work placeholder | outline §7 |
| Methods (appendix) | Drafted; B5 details, panel, σ_b prior, compute placeholders | spec, `src/` |

## Flags (to decide or fix)

1. **No biological finding yet.** Reviewers at either venue will ask what the method
   reveals. The uncertainty budget (§5.5) is the best candidate and is not built
   (≈50 CPU-h on M0). Recommend scoping it as its own task.
2. **AD and the "single compiled graph" are scoped.** Only route (ii) and all-ODE
   compositions are one `ODEProblem` with ForwardDiff through them. The hybrid is two
   stepped integrators and is not differentiable (clamped drain, integer counts). §2.4
   says so plainly; keep it that way.
3. **The dispatcher (B6) is unbuilt.** "Inference attaches to each data stream" is a
   design principle, not code. O3 is open: build a generic registry, or narrow the
   framework claim to one organism.
4. **The spec is stale in places** (do not fix from this branch; `spec/` belongs to track B):
   - `spec.md` §1 says 21 metabolic reactions; §0, §2, README and docs say 22. The
     paper uses 22.
   - `spec.md` §1 still gives the 60 s rebuild as the "likely reason" for the ENO
     null; A3 refuted it. The paper gives saturation as the reason for the weak gain,
     which is the measured part.
   - The 2026-10-23 framing is only in `spec/phases/approximate-method.md` §12.
   - `README.md` says "calibration has not yet passed"; V7 passed (`dev/scripts/v7_16a9b_result.md`).
5. **Caveats on the calibration claim:** M0 only; 23 of 50 V7 replicates have
   min bulk ESS < 100; F10 is one chain; FBA resolved by 0.020 SD. All stated in §5.1.
6. **The GTP mechanism is an interpretation.** Saturable regeneration (PGK3/PYK3)
   is not measured; §5.2 says so.
7. **Citations.** `references.bib` has 55 entries, each checked against a live source
   (`references_report.md`); about 20 unverified extras were dropped.
   - The 4D model's first author is Zane R. Thornburg (CrossRef,
     doi 10.1016/j.cell.2026.02.009). `README.md` and the outline say
     "Luthey-Schulten et al." and should be corrected on `main`.
   - The syn3a-ensemble companion is a `@misc` "in preparation"; update it when it
     is on arXiv.
   - Check Liu, Bayarri and Berger (2009) for author order against the PDF.
   - The RDME note uses Bayes factors; the model-selection note rejects them; the
     paper follows the latter.
8. **Companion paper.** syn3a-ensemble (*A whole minimal cell in a minute*) holds the
   speed claims. It is cited once, in §5.6 (cost).
9. **Template.** ICML preprint branding. Both venues reformat; PLOS CB needs an
   Author Summary, Cell Systems needs Highlights and an eTOC blurb (placeholder after
   the abstract).
10. **Figures.** Figs 2–4 reuse note figures; they need a common style for print.
    Figs 1, 3, 5, 6, 7 are not drawn.
11. **Table 1** (framework comparison) has no checked cells. Each needs a source.

## Next

- Fill §5.3 and the B5/B6 placeholders as track B lands (B5 → B6 → B7).
- Decide the uncertainty-budget task and O3.
- Draw Fig 1 (architecture) and Fig 3 (three graphs).
