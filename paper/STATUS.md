# Paper 1 draft: status

**Branch:** `paper-1-draft`, rebased onto `main` at `5ec56d3` on 2026-10-31 (KIM-44;
it first branched from `37e680c`). **Build:** `latexmk -pdf main.tex` in `paper/`
(TeX Live 2023; natbib, cleveref, xcolor present). At `5ec56d3` it builds to 12 pages
with no undefined references or citations.
**Working title:** *InferCell: inference-first whole-cell modelling*. It follows the
v0.3.0 claim (KIM-61).
**Venue:** PLOS Computational Biology, as a methods paper (KIM-15; `spec/spec.md` §12
2026-10-29, "v0.3.0 launch decisions"). The venue material is an Author Summary
(KIM-49). The two-column ICML 2026 template in `[preprint]` mode, copied from
`FNO_for_turing/paper/`, stays for drafting.
**Tracking:** Linear, project InferCell, milestone "M6 Paper 1 complete" (KIM-44 to
KIM-50, KIM-61).
**Numbers:** every quoted number is a `\pn{key}` in `numbers.tex`, with its source file.
Every source path it cites exists on main at `5ec56d3`.
**Placeholders:** red `\placeholder{...}` text and boxed `\placeholderfig` figures.
`grep -n placeholder sections/*.tex` lists them all.

## What is drafted

The text is unchanged by the rebase; only main moved under it. The last column names
the issue that clears each section's placeholders.

| Section | State | Sources | Issue |
|---|---|---|---|
| Front matter | Co-authors, affiliations and email are placeholders | — | KIM-49 |
| Abstract | Drafted; placeholders for B7, B8, the budget and the Author Summary | — | KIM-61 (claim, B7), KIM-46 (budget; B8 withdrawn), KIM-49 (Author Summary) |
| 1 Introduction | Drafted; Table 1 cells all placeholders; the approximate-route contribution and the last paragraph wait on B7 and B8 | `dev/notes/the-case-for-inference.md`, outline §1 | KIM-48 (Table 1), KIM-61 |
| 2 Design | Drafted; Fig 1 placeholder; B5 and B6 placeholders, though both are built (flag 3) | `src/`, docs, agent survey | KIM-48 (Fig 1); B5 and B6: none (flag 15) |
| 3 Model | Drafted; Fig 2 is `fig1c_state_graph_corea.pdf`; its caption is a placeholder | spec §0, §2; `src/organisms/coreA/m0.jl` | none (flag 15) |
| 4 The boundary | Drafted; Fig 3 placeholder; the FBA placeholder is stale (flag 12) | discussion note §2–§6 | KIM-48 (Fig 3); FBA: none (flag 15) |
| 5.1 Exact sampler | Drafted from V6, V7, F10; its F8 and F10 panels need one style | `dev/scripts/v6_*`, `v7_*`, `f10_*` result files | none (flag 15) |
| 5.2 Jensen gap | Drafted from B3, the rethink, the bias report, route (ii); one sentence waits on AV6 (B7) | `jensen_*`, `route_ii_89_result.md` | KIM-45 |
| 5.3 Approx vs exact | Placeholder. B7 has run | `dev/scripts/av6_89_b7_result.md`; `spec/phases/approximate-method.md` §12 2026-10-29, "the mixed route's status in one sentence" (KIM-19) | KIM-45 |
| 5.4 Full scale | Placeholder. B8 is parked outside paper 1, so the section is withdrawn | `approximate-method.md` §12 2026-10-29, "B8 parked" (KIM-17) | KIM-46 |
| 5.5 Uncertainty budget | Placeholder. M5's evidence run fills it | the KIM-40, KIM-41 and KIM-42 result notes (not yet written) | KIM-46 |
| 5.6 Cost | Drafted; the approximate row is a placeholder | 16-recovery §12; `av6_89_b7_result.md`, "KC (rule 7)" | KIM-45 (M0 cell), KIM-46 (B8 cell withdrawn) |
| 6 Discussion | Drafted; opening, B7, O3, asserted-priors and future-work placeholders | outline §7 | KIM-61 (opening, O3), KIM-45 (B7), KIM-46 (priors) |
| Methods (appendix) | Drafted; the venue's form, the panel and σ_b prior, B5's S and settings, and compute placeholders | spec, `src/` | none (flag 15) |
| Availability | Release DOI, the draws archive and funding are placeholders | — | KIM-54 (DOI), KIM-43 (archive); funding: none (flag 15) |

## Flags (to decide or fix)

1. **No biological finding: decided.** The v0.3.0 claim (O2, spec §1) does not claim
   one. The predictive and intervention demonstration on M0 (KIM-40, KIM-41) and the
   unscored reduced-model budget (KIM-42) are the paper's predictive result, in §5.5.
   Reviewers may still ask what the method reveals; the budget's sentence is the
   answer to aim for.
2. **AD and the "single compiled graph" are scoped.** Only route (ii) and all-ODE
   compositions are one `ODEProblem` with ForwardDiff through them. The hybrid is two
   stepped integrators and is not differentiable (clamped drain, integer counts). §2.4
   says so plainly; keep it that way.
3. **B5 and B6 are built; O3 is decided.** B5 is the mixed likelihood
   (`approximate-method.md` §12 2026-10-24, "B5 as run", #101). B6 is the dispatcher,
   `inference_route` in `src/inference.jl` (§12 2026-10-24, "B6 as built", #102). O3:
   no generic registry, and the wording is "a compositional inference framework
   demonstrated on reduced syn3A" (KIM-13). The text still calls B5 and B6 unbuilt in
   §2 and Methods (flag 15), and O3 open in the Discussion (KIM-61).
4. ~~**The spec is stale in places.**~~ Resolved on main:
   - §1's reaction count (22) and the ENO null's reason (saturation, not the rebuild)
     are corrected (KIM-24; `spec/spec.md` §12 2026-10-29, "§1's reaction count and
     the ENO null's reason"). The paper already matched.
   - The 2026-10-23 framing is superseded by the v0.3.0 claim, now in spec §1 (KIM-12).
   - `README.md` says M0 calibration has passed (KIM-23).
5. **Caveats on the calibration claim:** M0 only; 23 of 50 V7 replicates have
   min bulk ESS < 100; F10 is one chain; FBA resolved by 0.020 SD. All stated in §5.1.
6. **The GTP mechanism is an interpretation.** Saturable regeneration (PGK3/PYK3)
   is not measured; §5.2 says so.
7. **Citations.** `references.bib` has 55 entries, each checked against a live source
   (`references_report.md`); about 20 unverified extras were dropped.
   - The 4D model's first author is Zane R. Thornburg (CrossRef,
     doi 10.1016/j.cell.2026.02.009). `README.md` is corrected (KIM-23). The outline
     (`dev/notes/2026-10-19-paper-1-outline.md`, its references) still says
     "Luthey-Schulten et al. (2026)", and its "superseded in part" note, added by this
     branch, still calls the paper journal-agnostic.
   - The syn3a-ensemble companion is a `@misc` "in preparation"; update it when it
     is on arXiv.
   - Check Liu, Bayarri and Berger (2009) for author order against the PDF.
   - The RDME note uses Bayes factors; the model-selection note rejects them; the
     paper follows the latter.
8. **Companion paper.** syn3a-ensemble (*A whole minimal cell in a minute*) holds the
   speed claims. It is cited once, in §5.6 (cost). Full-cell stage 1 through its port
   is paper 2 and does not block this paper (KIM-15; KIM-55 to KIM-58).
9. **Template and venue material.** ICML preprint branding; PLOS CB reformats. The
   venue is chosen, so the abstract's "Author Summary or Highlights" placeholder becomes
   the Author Summary (KIM-49), and Methods goes after the Discussion, not as STAR
   Methods.
10. **Figures.** Figs 2–4 reuse note figures; they need a common style for print.
    Figs 1, 3, 5, 6 and 7 are not drawn (KIM-48). Fig 6 (F6) sits in §5.4, which
    KIM-46 withdraws with B8, so KIM-48 may not need it.
11. **Table 1** (framework comparison) has no checked cells. Each needs a source
    (KIM-48).
12. **§4's FBA placeholder is stale.** It says "K2 at 1,000 cells decides it
    (16b.5)". 16b.5 is budget-bound and was not run; at M0, K2's divergence is a
    50-cell read-off labelled as M0's, not K2's verdict (`spec/phases/16-recovery.md`
    §11 16b.5 and §12 "T3 at M0").
13. **The claim moved (O2).** The title, abstract, the introduction's contributions
    and the Discussion's opening follow spec §1 (KIM-61). A grep for "Core A′" in
    user-facing text is already clean at `5ec56d3`; only file names and `\pn` keys
    contain `corea`.
14. **B7.8's matched KS null may be degenerate** (KIM-67): its null chains are
    byte-identical to the base fit's. KIM-19's sentence already leaves B7.8 out. §5.3
    (KIM-45) should not lean on B7.8 until KIM-67 is settled.
15. **Placeholders no issue owns:** §2's B5 and B6 items; Methods' venue form, panel
    and σ_b prior, B5's S and sampler settings, and the compute table; Fig 2's caption
    and the restyle of Figs 2–4 (§5.1's F8 and F10 panels among them); §4's FBA sentence
    (flag 12); and funding in Availability.

## Next

- Can start now: KIM-61 (claim), KIM-47 (model-form table), KIM-48 (figures, Table 1),
  KIM-49 (front matter, Author Summary), and KIM-45 (§5.3) once KIM-67 is settled.
- Waits on M5's evidence run (KIM-40 to KIM-42): KIM-46 (§5.5, and withdrawing §5.4).
- Give flag 15's placeholders an owner.
- Last: KIM-50 (review and compliance passes).
