# Paper 1 framing: the decision

> **Reopened on 2026-11-26. This file is now the record of the old decision, not the current plan.**
> Tom's Decision log entry of 2026-11-26, "One paper, on the full syn3A model; the reduced-model
> preprint is held. Reopens the entry of 2026-11-25", names this decision and reopens it. That is
> condition 3 of [What would reopen this](#what-would-reopen-this). There is now one paper, on the
> full syn3A model. M0 and the reduced model are its verification tier. The manuscript in `paper/`
> is not posted as a preprint. It is the source for the new paper's framework and verification
> sections. v0.3.0 is tagged as software only. The technical risks in §3 below still hold, and they
> became the new plan's first gates (milestone M11). The plan is the Linear project InferCell. The
> text below is unchanged from 2026-11-25.

**Date:** 2026-11-25
**Decided by:** Tom
**Status:** REOPENED 2026-11-26 (see the note above). Was: DECIDED. Recorded in the Linear Decision log, entry of 2026-11-25
("Paper 1 claims the space on the reduced model; the full cell and the upgrades go to a later
flagship paper"). That entry is the record. This file is the reasoning behind it.

**Do not reopen this without new evidence.** The only things that would reopen it are listed in
[What would reopen this](#what-would-reopen-this). A new review asking for more scope is not one
of them.

---

## 1. The decision

1. **Paper 1 stays on the reduced syn3A model.** The full syn3A model is not in paper 1.
2. **Paper 1's job is to claim the space.** It puts the InferCell idea on record: composition with
   declared priors and reductions, inference across the ODE–stochastic boundary with a measured
   error, and scored probabilistic prediction.
3. **What claims the space is the arXiv preprint and the v0.3.0 release,** not acceptance at a
   journal. The bar for paper 1 is "sound and clearly scoped". It is not "survives a referee at a
   flagship venue".
4. **A later flagship paper takes the upgrades.** The full cell, real data, the degradation ladder
   (M10), the wider metabolic block and a large knockdown belong to that paper. Keeping them out
   of paper 1 also keeps that paper's novelty intact.
5. **The degradation ladder leaves paper 1.** The contribution sentence loses its last clause.
   This settles the parked release-vs-M10 question: the release and the preprint follow KIM-116,
   and do not wait for M10.

### The contribution sentence

**Was** (`story.md`, signed off 2026-11-23):

> InferCell infers a reduced whole-cell model of JCVI-syn3A across a coupled ODE–stochastic
> boundary and turns the posterior into predictions scored against held-out truth, then shows how
> both degrade as the synthetic data become realistic.

**Now** (the ladder clause is gone; the final wording waits for KIM-116, see §4):

> InferCell infers a reduced whole-cell model of JCVI-syn3A across a coupled ODE–stochastic
> boundary and turns the posterior into predictions scored against held-out truth.

## 2. The documents behind it

Read in this order. All are in `paper/` on the KIM-98 branch.

| File | What it is |
|---|---|
| [`manuscript_review.md`](manuscript_review.md) | The first review of the draft at `32b8496`. Six issues, and two strategic proposals: Point A (targeted inference) and Point B (drop the reduced model for the full cell) |
| [`checklist.md`](checklist.md) | The response, item by item: 35 items checked against the repository, with verdicts and tiers |
| [`response_to_review.md`](response_to_review.md) | The synthesis of that response. It scopes a full-cell paper and recommends three tiers |
| [`reply_to_response.md`](reply_to_response.md) | The reviewer's reply: five pushbacks and a four-gate protocol |
| this file | The decision, and what follows from it |

Also relevant:

- [`story.md`](story.md): the signed-off story of 2026-11-23, amended by this decision (local
  only, not in the repository).
- [`STATUS.md`](STATUS.md): the draft's open items.
- [`../dev/notes/2026-10-26-infercell-framing.md`](../dev/notes/2026-10-26-infercell-framing.md):
  the project's goal and its four stages. It calls the reduced model "validation infrastructure,
  not the main road" and sends the full cell through the port as a black box.
- [`../spec/prereg/kim-116.md`](../spec/prereg/kim-116.md): the modular arm's pre-registration
  (on the KIM-116 branch until it merges).
- Linear: project InferCell, the Decision log, and KIM-97 (Plan M10), which holds the three
  options of the release-vs-M10 question.

## 3. Why

### Why not the full cell in paper 1

- **It would make the least verified tier the headline.** The full cell's mean field does not
  close: 363 second-order reactions, 493 replication reactions, and rate rules with `round`, `min`
  and `max`. M0's exact check covers none of that, and no joint posterior is affordable there to
  referee a modular route.
- **The cheap route drops the paper's subject.** The data plug-in is cheap because it never
  crosses the ODE–stochastic boundary. On the full cell it also meets gene copy number as a new
  latent variable, because the reduced model cuts replication.
- **The full cell with no real data is the weakest combination.** `reply_to_response.md` is right
  about that. The fix chosen here is not to pay the full-cell cost in paper 1, not to add real
  data to it.
- **It contradicts the project's own framing.** The note of 2026-10-26 puts the full cell through
  the port, in stages that need no new inference method. A doubling-time check belongs there.
- **Cost.** One to two months of engineering before the first full-cell posterior, on the
  responder's own guess, added to a release held since 2026-11-09.
- **The scope has moved twice this month already.** The paper was centred on M0 (2026-11-07), then
  on the reduced model (2026-11-09). Each agent review adds scope. A paper that keeps moving
  claims nothing.
- **Portfolio.** PreCell already covers amortised inference on the full cell from bulk data.
  InferCell's distinct ground is the mechanistic likelihood, verified against an exact sampler.

### Why "claim the space" is enough for now

- The idea is new and the software works. A preprint with the release fixes priority on both.
- The upgrades in the review are a coherent second paper. Spending them on paper 1 delays the
  claim and leaves less for the flagship.

## 4. What is still open

These are open **inside** the decision. None of them reopens it.

| Item | Whose | Note |
|---|---|---|
| The submission date | Tom | "Claim the space" needs a date. Fix one, and treat whatever cannot meet it as flagship material. Not set yet |
| The framing sentence around the contribution | Tom, after KIM-116 | See the table below |
| Ten replicate truths by MAP and Laplace, 130–200 CPU-h | Tom | Optional. It removes "one dataset, one truth". Run it only after the promoter question, or it replicates a known bias. It lifts the non-goal "one truth per rung" |
| Paper 1's venue after the preprint | Tom | Not a flagship venue. PLOS Comput. Biol. was last named in Linear |
| Whether the 2026-11-30 target for M1–M7 moves | Tom | Follows from the date above |

### How KIM-116 sets the framing

The plug-in arm has passed its gate and is scored. The factor-cut arm is on its one extension
(job 21462423). The pre-registration reports each arm's standardised error against the truth.

| Outcome | Paper 1 says |
|---|---|
| The plug-in is close to the joint | Targeted inference is the practical route, at about 30% of the joint's gradient cost in stage 1 (11.3 s against 37.9 s, KIM-116 profile) |
| The plug-in differs and sits nearer the truth | Cutting protects a module from the other module's approximation error |
| The plug-in is worse | The joint was needed here. The current story stays |

A caution for the first row. On the reduced model the joint posterior has no exact reference, and
it carries the promoter shift. "Close to the joint" is not the same as "close to the truth".
Report both.

## 5. The minimum for paper 1

Three things, and no more without a reason that meets "What would reopen this". They are tracked
in Linear, milestone M6, as KIM-160 to KIM-169, with KIM-170 for the optional replicate truths.
`checklist.md` now opens with the map from its items to those issues.

1. **The promoter shift is diagnosed, or stated accurately.**
   - The 17 promoters sit about 2% low in log rate (*login-node* read of
     `pilot_coordinates.tsv`: mean −2.2%, precision-weighted −1.7%, 3.2 SE). Not for quoting
     until a tracked script produces it.
   - The suspect from `reply_to_response.md`, Pushback 3: the transcript factor reads the
     mean-field pools with no correction. It is a hypothesis. On M0 the mean-field GTP sits
     10–16% below the population mean, which pushes promoters up. At the saturated elasticity of
     0.044–0.051 that is also under 1%. The mechanism needs cells below K_d inside the fit window
     (saves to 3,120 s), and the 99.3% clip figure is from the 6,300 s scored runs.
   - The test already paid for is KIM-116's plug-in arm, which rebuilds the rates from observed
     pools. If its promoters are not low, the mean-field pools are the cause.
2. **KIM-116's result is filled in.** It decides the framing (§4).
3. **The text fixes** (tier A of `checklist.md`):
   - reword the "MAP" statement in §5: θ_ref is the MAP at δ = S = 0, not of the fitted
     likelihood;
   - fill the three KIM-90 todos (merged, PR #189, no cause named);
   - resolve the contradiction between the Discussion and `sec:info`;
   - add the positioning paragraphs and the missing citations (Gneiting and Raftery 2007 for
     CRPS; `lindsten2014particle`);
   - state the first-order restriction as a limit, and say why ENO and FBA;
   - move or replace the §2 code listing; fix the two schematics;
   - remove §7, F7, the ladder appendix and every ladder clause in the abstract, introduction and
     discussion.

## 6. What goes to the flagship paper

Not in paper 1, by this decision:

- the full syn3A model, by any route;
- real data, as a fit or as a check;
- the degradation ladder (M10), including the snapshot scRNA-seq rung;
- a wider block of free metabolic constants, and the Hessian eigen-spectrum probe;
- a second, resolvable intervention;
- fitting the reduced model to data from the full simulator;
- the four-gate protocol's Gate 4 (full-cell timing measurements).

Paper 1's Discussion can name these as what comes next. The list of record is the Linear
document "Paper 2: deferred scope".

## What would reopen this

Only these:

1. KIM-116 or the promoter diagnosis shows the mixed route is wrong in a way that makes §5–§6 of
   the paper unsound. That would change what paper 1 can claim. It would not move the full cell
   into it.
2. Another group publishes compositional Bayesian inference on a whole-cell model first.
3. Tom decides to, and writes a new Decision log entry that names this one.

## Done since, and not done

Done on 2026-11-26:

- Linear matches the decision: the project summary and description, milestones M6, M7 and M10,
  KIM-97 (deferred), KIM-98 (decomposed), new issues KIM-160 to KIM-170, and the document
  "Paper 2: deferred scope".
- `CLAUDE.md` (local, not tracked) carries the scope and points here.

Not done:

- No edit to the manuscript. KIM-160 to KIM-169 are all open.
- No submission date is set.
- The older paper issues from the M0 draft (KIM-45 to KIM-49, KIM-61, KIM-85) are not re-scoped.
- This file, `checklist.md`, `response_to_review.md` and `reply_to_response.md` are uncommitted.
  Whether they go into the public repository is not decided.
