# Review checklist: paper 1

**Source:** [`manuscript_review.md`](manuscript_review.md) (2026-11-25), turned into one item per point.
**Written:** 2026-11-25, against the draft at `32b8496` on the KIM-98 branch. Nothing here is done yet,
and nothing in the manuscript, Linear or the Decision log was changed in writing it.

**How to read an item.** *Review* quotes the point with its line in the review. *Checked* says what I
verified in the repository, and what I could not. *Verdict* is mine: **agree**, **partly**, or
**disagree**. *Tier* is one of:

- **A, fix regardless:** text, figures or analysis of data we already hold. No new runs.
- **B, compute:** needs new runs. Worth it only if the venue asks for it.
- **C, direction:** changes the claim, the scope or the plan. Your decision, then a Decision log entry.

Numbers marked *login-node* are my own quick reads of archived draws. They are indicative. Recompute
them in a tracked script before quoting any of them.

> **Status, 2026-11-26. Read this before anything below.** The framing is decided
> ([`framing_decision.md`](framing_decision.md); Decision log 2026-11-25). This file is now the
> evidence behind the Linear issues, not a list of open decisions.
>
> - **Tier A is tracked in Linear,** milestone M6: 1a → KIM-160; 1b, 4a, H2 → KIM-161; A1, 6a
>   (promotion), 6c → KIM-162; the ladder's removal → KIM-163; 3a, 3b, 6a → KIM-164; 4b → KIM-165;
>   6d → KIM-166; F1, F2 → KIM-167; 6b → KIM-168; M1, H1, H3 → KIM-169. Work from the issues.
> - **Settled:** V1 (preprint first; the journal after it is open), 5a (the ladder leaves paper 1),
>   5b, 5c, B1 (the reduced model stays), RC1–RC3. "Decisions only you can make" and "What I would
>   do, in order" below are history.
> - **Deferred to paper 2:** 2b, the Hessian probe in 3a, 5a's rungs, 5b, 5c, the promotion in 6a,
>   R1, RC1–RC3. They are listed in the Linear document "Paper 2: deferred scope".
> - **Still open:** 2a, replicate truths, as KIM-170.
> - **Changed:** M1's thirteen ladder todos are deleted, not filled. 4a step 4 goes with the
>   ladder. 6c waits on KIM-116 only.

---

## Read this first

### 1. Where I disagree with the review

| Item | The review says | My view |
|---|---|---|
| [B1](#b1) | Drop the reduced model; go to M0 + full syn3A | **Disagree.** The full cell is not an InferCell composition, so the mixed route cannot run on it. PreCell already calibrates the full model, though only through 28 group multipliers. It also reverses the 2026-11-09 decision |
| [5c](#5c) | Add a real-data case study | **Disagree for paper 1.** It is a stated non-goal, and the reduced model cannot be compared with real data by construction |
| [5b](#5b) | Fit the reduced model to data from the full simulator | **Disagree for paper 1.** A good idea in the wrong paper: it is the transfer test that paper 2 owns |
| [6a](#6a) | Promote the three graphs to the main text as a flagship result | **Partly disagree.** Its key direction has no exact result, and its saturation premise fails on the reduced model. Fix the contradiction; promote only if KIM-116 supports it |
| [5a](#5a) | M10 is mandatory before submission | **Partly disagree.** True for *Nat. Comput. Sci.*; not for a methods venue. This is the parked release-vs-M10 question, not a settled fact |
| [6c](#6c) | End the abstract on a forward-looking advance | **Partly disagree.** The review misquotes the ending, and its replacement would claim things the paper does not show |
| [6b](#6b) | Delete "We did not read the supplements" | **Partly disagree.** Only after someone reads them. I would also keep "the comparison favours it" |
| [4b](#4b) | Say neural SBI is prohibitive here; integrate the 27 uncited references | **Partly disagree.** PreCell's draft does amortised inference on the full model from bulk data (28 group multipliers), so a flat "prohibitive" is unsafe. Cite what the positioning needs, not 27 entries |
| [3a](#3a) | Run a NUTS stress test with 10–30 metabolic constants free | **Partly disagree.** A fair question. Answer it with a cheap Hessian probe and a stated limit, not a new posterior |
| [V1](#v1) | Target *Nat. Comput. Sci.* | **Partly disagree.** The review never weighs the venue last named in Linear (PLOS Comput. Biol., methods) or a preprint. Most of "will block acceptance" is calibrated to venues the project has not chosen |
| [A1](#a1) | The paper "accidentally argues" for monolithic joint inference | **Partly disagree on the diagnosis, agree on the direction.** The joint posterior is there on purpose, as the reference for the modular arm |

### 2. Facts in the review that are stale or wrong

- **KIM-90 is finished.** PR #189 merged on 2026-11-23 (`ae731c6`). The review treats it as pending
  ("ensuring KIM-90 *diagnoses* the M0 gap"). It ended with no cause named. Its three `\todo`s can be
  filled once main is merged into this branch. See [1a](#1a).
- **The disjoint-cohort test was considered and declined.** It needs a new exact run of 508–678 CPU-h
  (your call, 2026-11-23, in the KIM-90 note). See [R1](#r1).
- **The abstract quote is spliced.** "The data are synthetic, from one truth" is not the last
  sentence. See [6c](#6c).
- **`fig_model` already prints the 1 s and 60 s clocks** on its arrows. See [F2](#f2).
- **`tab:reducedcost`'s caption already itemises the 24.3 CPU-h probe.** See [H2](#h2).
- **"99.3% GTP clipping is a pathology the reduction likely introduces"** is unverified. The clip rule
  is the published model's. Nobody has checked how often the full cell clips. See [B1](#b1).
- **"~5,000+ parameters"** for a whole-cell model has no source in the repository.

Everything else I checked is accurate: the todo counts (22, 2, 1, 1), 27 of 56 references uncited, and
every quoted number (3.3 SE, −0.81, 0.4%, 1.09, 0.96–1.08, 2.3 and 1.4 SD, 99.3%, 93 reductions, the
40,000 CPU-h of exact runs).

### 3. Two things I found while checking

1. **The promoter shift is harder to explain as noise than the draft says** ([1b](#1b)). The 3.3 SE
   assumes independent errors. With the posterior's own correlations it is 3.1 SE, and 4 of the 1,600 posterior
   draws put the promoters' mean standardised error at or above zero (*login-node*, from
   `pilot_draws.tsv`). Noise shared across promoters would show as posterior correlation between
   them, and their mean pairwise correlation is 0.01. So either this dataset is roughly a 1-in-400 draw (four draws is a noisy
   tail), or the posterior is miscalibrated. "One dataset cannot separate the two possible causes" is too even.
2. **The "MAP" is the MAP of a different target, and ENO and FBA are not skewed** ([4a](#4a)).
   - The pilot script computes θ_ref as "the bounded MAP at δ = S = 0"
     (`dev/scripts/corea_pilot_kim92.jl:18`), as Methods says (`app_methods.tex:103`). The pilot
     note calls it "the MAP of the same likelihood", and §5 follows the note. The fitted likelihood
     has δ and S in it; θ_ref's does not.
   - ENO's and FBA's marginals have skewness 0.05 and 0.06 (*login-node*). So "skewed" is out. The
     likeliest reading of their 2.3 and 1.4 SD offsets is that δ and S move them, not that the
     optimiser stopped short. This is my inference; a MAP of the fitted likelihood would settle it.
   - "The MAP shares the shift, so it comes from the likelihood, not the sampler"
     (`05-inference.tex:52–54`) needs rewording. What the numbers show is that the promoters sit
     0.86 SD low before δ and S enter, and 0.81 after. So the values of δ and S are not what puts
     the promoters low.

### 4. What I would do, in order

1. **Decide the venue** ([V1](#v1)). It sets which tier-B items matter.
2. **Do tier A now.** Fill the KIM-90 todos, resolve the Discussion contradiction, reword the
   promoter paragraph, add the positioning paragraphs and the two missing citations, fix the code
   listing and the schematics.
3. **Run the cheap replicate check** ([2a](#2a)): MAP and Laplace on about ten new truths, roughly
   16 CPU-h each. It answers Issues 1 and 2 together for under 200 CPU-h.
4. **Add one resolvable intervention** ([2b](#2b)), after a dose scan at θ\*. This is the review's
   best point. It is the one large run here: of order 2,000–5,000 CPU-h.
5. **Wait for KIM-116,** then settle Point A's reframing ([A1](#a1)) and the M10 question
   ([5a](#5a)) together. Both hang on how close the modular arm sits to the joint.
6. **Do not pivot to the full cell in paper 1** ([B1](#b1)).

### 5. Decisions only you can make

- [ ] **Venue:** PLOS Comput. Biol. (last named in Linear), *Nat. Comput. Sci.*, or preprint first.
- [ ] **M10:** wait, narrow the claim, or one probe rung (the three options on KIM-97).
- [ ] **Replicate truths:** lifts the non-goal "one truth per rung". Yes or no, and in which form.
- [ ] **A second intervention:** it is post hoc, since the dose was frozen in KIM-39.
- [ ] **Point A:** change the contribution sentence, or only the framing around it.
- [ ] **Point B:** confirm the reduced model stays.

---

## §1 Venue

<a id="v1"></a>
- [ ] **V1. Choose the target venue.**
  - **Review:** "Make *Nature Computational Science* the primary target unless a real-data
    application or biological discovery… is added specifically for *Cell Systems*."
    ([L18](manuscript_review.md#L18); table at [L13–16](manuscript_review.md#L13))
  - **Checked:** `STATUS.md` says "Venue: open" and `story.md` says "none (preprint)". Linear's M6
    milestone still says "PLOS Comp Biol, methods", which is older text: the PLOS template was
    dropped on 2026-11-07. A keyword search of the Decision log finds no venue entry (I did not
    read the log in full). The review weighs neither a methods journal nor a preprint.
  - **Verdict: partly.** I agree *Cell Systems* is a poor fit: "no real data, and no biological
    finding about syn3A" are written non-goals. The *Nat. Comput. Sci.* verdict ("reject with
    encouragement") is the reviewer's guess, and reaching it costs most of tier B. Two of your
    earlier positions pull apart here. On 2026-11-05 the framework launch was the deliverable and
    the paper an add-on. On 2026-11-09 you held the release for a stronger paper, with no cost cap
    on the pilot. The venue choice says which of the two governs now.
  - **Tier C.** Decide this first. If the venue stays a methods journal or a preprint, items 2a, 2b
    and 5a become choices and stop being requirements.

## §2 What to protect

- [ ] **P1. Keep the problem statement as the opening.**
  - **Review:** "A genuinely important, longstanding open problem." ([L24](manuscript_review.md#L24))
  - **Verdict: agree.** No action. One caution: the review says the work "directly answers" the
    DREAM8 challenge. The draft does not claim that, and should not start to.

- [ ] **P2. Keep the rigour visible when cutting for length.**
  - **Review:** "Exceptional statistical rigor and reproducibility." ([L25](manuscript_review.md#L25))
  - **Verdict: agree.** No action. Keep the sentences that bound it: SBC ran on M0 only, and 17 of 50
    mixed and 23 of 50 exact replicates fall short of the convergence targets (`04-m0.tex:32–38`).

- [ ] **P3. Keep declared reductions and the cut-versus-reduction distinction prominent.**
  - **Review:** "Sharp, memorable conceptual contributions." ([L26–29](manuscript_review.md#L26))
  - **Verdict: agree on the first two, partly on the third.** The three graphs are judged in
    [6a](#6a). Note that Point B ([B1](#b1)) would remove the model these reductions describe.

- [ ] **P4. Leave Figs. 3–6 and `fig_trajectories` alone.**
  - **Review:** "Publication-grade data figures." ([L30](manuscript_review.md#L30))
  - **Verdict: agree.** No action.

## §3 Issue 1: biases in the mixed route

<a id="1a"></a>
- [ ] **1a. State the M0 gap with what KIM-90 found.**
  - **Review:** "3 of 6 parameters sit 0.42–0.56 SD from the exact posterior, and the text states
    'We have not found the gap's cause.'" ([L39](manuscript_review.md#L39))
  - **Checked:** accurate (`04-m0.tex:63–78`). But KIM-90 merged on 2026-11-23. Its verdict: no cause
    named. Target identity and MCMC error are ruled out (the gap is 7–11 z beyond both). Varying the
    reference closes about 0.17 SD of the ptsG promoter and up to 0.16 of the decay rate. Shared
    cells are indeterminate. Double counting is real as a width effect (the bulk carries 0.19–0.66
    of the information the route credits it with), and its effect on location is unmeasured.
  - **Verdict: agree it is a weak point; the remedy has been tried.** The review asks for a
    diagnosis. A full issue of diagnosis (nine checks, four refit arms) did not produce one.
  - **Tier A.** This branch was cut at `24273c7`, before #189, so merge main first. Then fill the
    three KIM-90 todos (`04-m0.tex:74`, `app_m0.tex` "The remaining gap", `08-discussion.tex:61`)
    from `gap_diag_kim90_result.md`. `story.md` already records the §4 wording change as agreed on
    2026-11-23; Linear's story line still says "the remaining gap is diagnosed" and needs the same
    change. Say that the gap is a location shift beyond sampler error, so the text must not call
    it noise.

<a id="1b"></a>
- [ ] **1b. Rewrite the promoter-shift paragraph.**
  - **Review:** "Leaving an unexplained 3.3 SE systematic bias in the headline fit, and explicitly
    stating 'a second dataset would tell them apart, and we have not run one'… reads as unfinished
    rather than thorough." ([L40](manuscript_review.md#L40), [L43](manuscript_review.md#L43))
  - **Checked:** accurate (`05-inference.tex:47–58`). See finding 1 above: 3.1 SE with correlations,
    mean promoter–promoter correlation 0.01, 4 of 1,600 draws at or above zero (*login-node*).
  - **Verdict: agree, and it is somewhat worse than the review says.** The shared-noise reading has
    little support in the posterior itself. KIM-90's double counting predicts a posterior that is
    too narrow, which makes any offset look larger in SD units; the paragraph should say so.
  - **Tier A** for the rewording: give the correlation-adjusted figure, drop the even-handed "two
    possible causes", and link to the M0 double counting. **Tier B** to settle it: see [2a](#2a).

- [ ] **1c. Tie the pool misses to the same question.**
  - **Review:** "10 of 17 metabolite pools miss by 0.9–1.4 SD (G3P +1.39, DHAP +1.23, plus
    systematic energy-charge splitting)." ([L41](manuscript_review.md#L41))
  - **Checked:** accurate (`06-prediction.tex:56–72`). The draft already reports this as its finding
    and already says the link to the promoter shift is post hoc and untested
    (`08-discussion.tex:24–27`).
  - **Verdict: agree with the facts, no new action beyond 1b and 2a.** I would not soften it. A
    scored miss is the paper's point.

<a id="1d"></a>
- [ ] **1d. Decide whether to diagnose the reduced-model bias: noise, fixed δ and clipping, or double counting.**
  - **Review:** "Diagnosing whether the −0.81 SD promoter shift and 10/17 pool shifts come from (a)
    a single unlucky noise draw, (b) fixed δ / GTP clipping below K_d, or (c) double-counting across
    streams would turn a vulnerability into one of the strongest methodological findings."
    ([L45](manuscript_review.md#L45))
  - **Checked:** the premises hold. 99.3% of baseline cells clip (`kim95.clip`). The mean field does
    not clip (`app_methods.tex:93`). δ and S are fixed at one reference value. The route ignores
    that both factors read the same cells (`app_methods.tex:101`).
  - **Verdict: partly.** The three hypotheses are the right ones. "One of the strongest findings" is
    optimistic: on M0 the same exercise ended with no named cause. Each hypothesis has a cheap test:
    - (a) replicate truths by MAP and Laplace, [2a](#2a);
    - (b) recompute δ and S at the posterior mean and refit, the M0 "postmean" arm. One MAP is about
      3 CPU-h; the chains of a full posterior are about 450. Finding 2 already weakens the "fixed δ"
      half of (b) for the promoters: they sit as low with δ = S = 0. The clipping half is untouched,
      because the mean field does not clip whatever δ is;
    - (c) fit transcripts from cells 1–100 and bulk from cells 101–200. On the reduced model there is
      no exact reference, so this shows only whether the posterior moves or widens.
  - **Tier B.** I would run (a) first. If the shift repeats across truths, (b) and (c) at MAP level
    are worth their few CPU-h each. If it does not, stop.

## §3 Issue 2: one dataset, and the knockdown

<a id="2a"></a>
- [ ] **2a. Replicate truths on the reduced model.**
  - **Review:** "Running even 3–5 replicate truths on the reduced model (or using the modular/data
    plug-in route if KIM-116 validates it as cheap and accurate) would separate random dataset noise
    from systematic bias." ([L48](manuscript_review.md#L48))
  - **Checked:** one joint posterior costs 492 CPU-h and 5.5 days on the critical path. Chains are
    serial, so replicates run side by side. Scoring one at KIM-95's sizes cost 5,447 CPU-h of cell
    work. "No SBC or coverage rate on the reduced model; one truth per rung" is a written non-goal.
  - **Verdict: agree this is the draft's largest weakness; partly on the remedy.** Full NUTS
    replicates are the expensive form. The cheap form uses the review's own observation in
    [4a](#4a): Laplace matches NUTS widths to 4–8%. The pilot's prep stage cost 7.3 CPU-h (MAP 2.7 h,
    reference ensemble 4.5 h), the Hessian 0.7, the dataset 4.5. So MAP plus Laplace on ten new
    truths is roughly 130–200 CPU-h, with no scoring. That gives a recovery rate for the promoters
    over eleven truths, which is exactly what 1b needs.
  - **Caveat:** by finding 2, the pilot's MAP is of the δ = S = 0 target. The replicates need the
    MAP of the fitted likelihood, with δ and S in place: a second optimisation per dataset, about
    3 CPU-h more. The Laplace-to-NUTS width match was measured at θ_ref, so check it holds at the
    new MAP on the one dataset that has both.
  - **Tier B, recommended in the cheap form.** It lifts a non-goal, so it needs a Decision log
    entry and a short pre-registration of what counts as "the shift repeats".

<a id="2b"></a>
- [ ] **2b. Add an intervention the predictions can resolve.**
  - **Review:** "The headline out-of-sample intervention (30% ptsG knockdown) changes Q by only
    0.4%… completely buried in cell noise (MC share = 1.09)." Fix: "70–90% CRISPRi-style knockdown
    of ptsG, eno, or gapd, or a step change in extracellular glucose… 15–40% shift."
    ([L49–50](manuscript_review.md#L49))
  - **Checked:** accurate (`06-prediction.tex:29–44`; `kim95.dqshare` 0.4, `kim95.dqmc` 1.09). The
    draft concedes it: "covering Q_kd tests nearly the same thing as covering Q." The dose was
    frozen in KIM-39 from M0 pilots, before evidence.
  - **Verdict: agree. This is the review's strongest point and the fix is bounded.** As it stands,
    the paper's only intervention tests nothing the baseline does not.
  - **How I would do it:**
    1. Keep the frozen 30% result as it is. Do not replace it.
    2. Scan dose and target on the true simulator at θ\* only. This is forward simulation at 32 s a
       cell. Find a knockdown whose ΔQ clears cell noise.
    3. Pre-register the new intervention, then run one scored arm. A larger effect needs fewer
       cells than KIM-95's 973 per draw.
    4. Label it post hoc in the text, with a Decision log entry.
  - **Cost:** the scan is small. The scored arm is of order 2,000–5,000 CPU-h at KIM-95's sizes
    (KIM-95's cell work was 5,447 CPU-h for baseline and knockdown together), and less if the
    larger effect lets the sizing rules pass with fewer cells. Not costed properly yet.
  - **Cautions on the review's specifics:** the knockdown scales a promoter strength, so any of the
    17 genes works, and a glucose step needs a new intervention type. A strong knockdown may starve
    cells, and the scored runs apply no viability rule. "15–40%" is the reviewer's number. The model
    makes no claim about syn3A, so the test is statistical resolution, not biology.
  - **Tier B, recommended.**

## §3 Issue 3: parameter scope and first-order generality

<a id="3a"></a>
- [ ] **3a. Justify ENO and FBA, and answer the sloppiness question.**
  - **Review:** "18 of 21 free parameters are constitutive gene-expression rates; only 2 (ENO, FBA
    k_cat) are nonlinear metabolic constants… inferring 17 linear promoter rates and 1 decay rate is
    nearly a decoupled textbook birth–death problem. Does the mixed route still work when 10–30
    sloppy, correlated metabolic parameters are free?" ([L55–56](manuscript_review.md#L55))
  - **Checked:** accurate, and "~62 other constants" matches the 62 releasable constants (Spec D1,
    `traj.nreleased`). The draft states the free set but never says why ENO and FBA
    (`03-testsystem.tex:52–55`). The reason exists: Spec D1 took the union of the budget's inferred
    and held sets, and ENO and FBA came from the phase-16 targets. Spec D1 lists "all 62 releasable
    constants" as "a later rung".
  - **Verdict: agree on the justification; partly on the stress test.** The "decoupled birth–death"
    remark is fair, and KIM-116's factor cut measures exactly that. If the cut lands close to the
    joint, the paper should say plainly that the promoters are a near-decoupled problem. That is an
    honest result and it supports Point A.
  - **Tier A:** add two sentences on why these two constants, and name the restriction in Limits.
  - **Tier B, optional:** rather than a new posterior, free 10–15 glycolytic constants and report
    the Hessian's eigenvalue spectrum at the MAP. That shows the sloppy directions for a few CPU-h.
    Cost scales with parameter count under forward-mode AD, and the unit-metric condition number is
    already 528 at 21 parameters. I would not run NUTS on a wider set for paper 1.

<a id="3b"></a>
- [ ] **3b. State the first-order restriction as a limit, and say what would replace it.**
  - **Review:** "Exact closure of the mean-field ODEs and the exact birth–death transcript factor
    rely on gene expression being strictly first-order… Fix: show how the compositional interface
    accommodates nonlinear gene regulation (FSP, LNA/moment closure, or a local neural likelihood
    surrogate)." ([L57–58](manuscript_review.md#L57))
  - **Checked:** accurate (`02-infercell.tex:104–109`, `:189`). The Limits list in the Discussion
    does not include it. "One organism, one reduction" is the nearest entry.
  - **Verdict: agree.** A referee from *E. coli* or yeast will ask at once.
  - **Tier A.** Add a Limits bullet and one paragraph. Describe which parts are specific (the
    closed mean field, the birth–death factor) and which are general (the typed edges, the route
    dispatch). Name FSP, LNA and surrogates as candidates. Do not write that the interface
    "accommodates" them, because none is implemented. One correction to the review: first order is a
    property of the published model, not of the organism.

## §3 Issue 4: positioning against alternatives

<a id="4a"></a>
- [ ] **4a. Say why NUTS, given that Laplace matches it.**
  - **Review:** "Why pay 492 CPU-h (5.5 days wall-clock) for NUTS if Laplace matches it within
    4–8%?… (a) the MAP optimizer missed the ENO and FBA posterior means by 2.3 SD and 1.4 SD… (b)
    along the degradation ladder, sparser/noisier data will likely produce non-Gaussian ridges."
    ([L63](manuscript_review.md#L63))
  - **Checked:** accurate (`05-inference.tex:71–75`). See finding 2: ENO and FBA are symmetric, so
    skew does not explain the offsets, and the "MAP" is θ_ref, the MAP at δ = S = 0. The review's
    answer (a), "the optimizer missed", inherits the draft's mislabel.
  - **Verdict: agree the question must be answered; partly on the answer.** The honest answer is:
    on this generous dataset Laplace gives the right widths, and we know that only because NUTS
    ran. That is useful to a user, and it is what makes [2a](#2a) cheap. Suggestion (b) is a
    prediction, not a result.
  - **Tier A.**
    1. Correct `05-inference.tex:52–54` and `:71–75`, and the pilot note, to say θ_ref is the MAP at
       δ = S = 0 (finding 2). The Laplace SDs are also taken at θ_ref.
    2. Add the skewness numbers from a tracked script.
    3. Write the answer above into the cost paragraph.
    4. If M10 runs, the MAP and Hessian exist at every rung for the warm metric, so
       Laplace-versus-NUTS along the ladder costs nothing extra. Note it on KIM-97.

<a id="4b"></a>
- [ ] **4b. Add two positioning paragraphs and the citations they need.**
  - **Review:** "Where does InferCell sit relative to neural SBI, LNA, and linked emulators? 27 of
    the 56 entries in `references.bib` are never cited." ([L64–65](manuscript_review.md#L64))
  - **Checked:** exactly 27 of 56 are uncited, and every key the review names is among them.
    "Bayesian balancing of metabolic parameters" appears in the Introduction with no citation to
    `lubitz2010parameter` (`01-intro.tex:40`, `:63–66`). Ancestor sampling is uncited at
    `app_methods.tex:67`. Gneiting and Raftery is not in the bibliography at all.
  - **Verdict: agree on the gap; partly on the content.**
    - Do not write flatly that neural SBI is "prohibitive" here. PreCell's draft abstract describes
      one diffusion transformer, trained on simulations from the fast port, that returns posteriors
      over the full model's kinetic parameters from bulk measurements. Its θ are 28 multipliers,
      one per parameter group, within a factor of 2 of the published values, and it uses no
      single-cell trajectories. So it does not contradict the review's narrower claim (21+
      individual parameters with 200-cell trajectories). State the difference that is
      safe: a likelihood with an exact reference on M0, against amortisation with none. Cite
      PreCell as a companion once it can be cited.
    - "No neural SBI comparison" is a non-goal for experiments. A positioning paragraph does not
      break it.
    - Cite what the two paragraphs need. Uncited entries do not print, so there is nothing to
      "integrate". Check `references_report.md` before adding any.
  - **Tier A.** Two paragraphs (SBI and amortisation; LNA and moment closure against the exact
    birth–death factor; linked emulators and multimodel inference), plus `lubitz2010parameter`,
    `lindsten2014particle`, and a new Gneiting and Raftery (2007) entry for both CRPS sentences.
    Stage 7 (`/check-refs`) is still to run and covers the rest.

## §3 Issue 5: the synthetic, well-specified gap

<a id="5a"></a>
- [ ] **5a. Decide M10, and whether its rungs are the right ones.**
  - **Review:** "Completing M10 is mandatory before submission to either journal. Realistic rungs
    should include destructive/snapshot scRNA-seq…, capture dropout / binomial thinning, sparser
    sampling intervals, and subsetted metabolite panels." ([L70](manuscript_review.md#L70))
  - **Checked:** M10 has not started and KIM-97 is blocked by KIM-116. Spec D6 (draft) has five
    rungs: (a) generous, (b) transcripts for two genes, (c) bulk only, (d) capture efficiency, (e)
    unknown cell number. Thirteen of the 22 todos wait on M10. Your estimate on KIM-97: about two
    weeks without scored predictions, three to five with.
  - **Verdict: partly disagree on "mandatory"; agree on the rungs.**
    - "Mandatory" is the parked release-vs-M10 question. Option 2 on KIM-97, a narrower paper on
      the generous rung, is an honest scope for a methods venue. For *Nat. Comput. Sci.* I agree
      some ladder is needed.
    - The snapshot rung is the review's best addition and D6 lacks it. Real scRNA-seq is
      destructive. The draft itself says time-resolved records are what separate decay from
      synthesis (`app_m0.tex:171–174`), so a snapshot rung tests the generous design's most
      unrealistic feature. It needs a new transcript factor on marginals, not transitions.
    - One correction: the 98.9% is the decay rate's variance reduction under transcripts only on
      M0, not a measured share "driven by" autocorrelation.
  - **Tier C.** Decide with KIM-116 in hand. If M10 runs, add a snapshot rung to D6 and consider
    dropping (b) for it.

<a id="5b"></a>
- [ ] **5b. Structural misspecification test with the full simulator.**
  - **Review:** "Generate synthetic data from the *full* cell model and fit the *reduced* model."
    ([L72](manuscript_review.md#L72))
  - **Checked:** SpeedyCell.jl simulates the full well-stirred model and records counts per minute,
    so the data could be made. But the reduced model removes about 92% of the cell and "cannot
    reproduce the doubling time, the metabolite levels or the proteome" (`03-testsystem.tex:43–48`).
    The fit would be grossly misspecified, not subtly. `STATUS.md`'s standing flag says the transfer
    "is assumed, not tested, and a test is likely paper 2". I did not check how the two models'
    observables map onto each other.
  - **Verdict: disagree for paper 1.** It is a good experiment and it belongs to the full-cell
    milestone. In paper 1 it would most likely produce a failed predictive check and a long
    explanation.
  - **Tier C.** If a misspecification result is wanted in paper 1, a smaller one is available:
    generate data with one declared reduction switched and refit. Otherwise record it under paper 2.

<a id="5c"></a>
- [ ] **5c. Real-data case study.**
  - **Review:** "Fit a real experimental dataset (e.g., published JCVI-syn3A
    proteomics/metabolomics/growth constraints…, or a well-characterized… dataset in E. coli /
    yeast)." ([L73](manuscript_review.md#L73))
  - **Verdict: disagree for paper 1.** It is a written non-goal. The reduced model cannot be set
    against syn3A data by design, and the code refuses to report a doubling time. A second organism
    is a second implementation, which the Limits list already names as future work. This option is
    only live if *Cell Systems* is the target, and I would not make it the target.
  - **Tier C.** No action unless V1 says otherwise.

## §3 Issue 6: buried results and tone

<a id="6a"></a>
- [ ] **6a. Resolve the three-graphs contradiction; decide on promotion.**
  - **Review:** "Promote the 'Three Graphs' and Elasticity Working Rule from Appendix B to the Main
    Text… a modeller can diagnose cross-module identifiability analytically before running a single
    MCMC chain. Burying this… is both a logical contradiction… and a missed opportunity."
    ([L78](manuscript_review.md#L78))
  - **Checked:** the contradiction is real and already in `STATUS.md`. `app_m0.tex:102` says "No
    claim in the main text rests on this subsection", and `08-discussion.tex:29–33` draws a working
    rule from it. But the evidence under the rule is thin:
    - metabolites to gene expression is shown by a linearised bound and the mixed route only. The
      exact fit did not mix (ESS 39), and the draft says "we do not claim it there";
    - transcripts to metabolism is a null read-off on M0. The planned 1,000-cell test did not run;
    - the rule rests on saturated rate laws. On the reduced model 99.3% of cells clip GTP, below
      K_d, and the Discussion says the rule "is not tested there".
  - **Verdict: partly disagree.** Fix the contradiction; do not promote it as a headline. The
    review's summary of the result is stronger than what was measured. If it leads the paper, a
    referee will ask for the exact-route result that is on record as missing.
  - **Tier A:** make one side give way. I would keep the appendix sentence and reword the Discussion
    paragraph as a hypothesis from M0.
  - **Tier C:** KIM-116's factor cut is the rule's test on the reduced model. If it lands close to
    the joint, a short main-text paragraph with `tab:gain` is earned, presented as the reason the
    cut works. Decide then.

<a id="6b"></a>
- [ ] **6b. Read the supplements, then fix Table 1's footnote.**
  - **Review:** "Remove 'We did not read the supplements' (verify the supplements of those 4 papers
    so Table 1 is bulletproof) and reframe neutrally… Update 'for the 2026 syn3A model we read the
    preprint' once checked against the published *Cell* 2026 paper." ([L80](manuscript_review.md#L80))
  - **Checked:** accurate (`01-intro.tex:93–101`). `tables/frameworks_sources.md` lists the unread
    supplements and says the *Cell* text returned 403, so the bioRxiv version was read.
    `STATUS.md` already has "Spot-check Table 1, read the supplements" against your name. Karr's and
    Macklin's supplements hold their fitting methods, so two cells could change.
  - **Verdict: partly.** The sentence goes only when the reading is done. I would keep "These
    columns are InferCell's own design aims, so the comparison favours it." It answers the first
    objection a referee will raise, and it costs one line.
  - **Tier A,** but it is reading time, yours or an agent's with your spot-check.

<a id="6c"></a>
- [ ] **6c. Revisit the abstract's ending once the M10 call is made.**
  - **Review:** "Ending on 'The data are synthetic, from one truth... show where it is wrong' lands
    flat… close on the forward-looking scientific advance (what compositional inference enables for
    experimental design, uncertainty-aware perturbation prediction, and scaling to full whole-cell
    models)." ([L81](manuscript_review.md#L81))
  - **Checked:** the quote splices two sentences. The abstract ends: "A whole-cell posterior can
    thus be carried to scored predictions that show where it is wrong. Because every departure from
    the published model is declared, each score states what it does not cover."
    (`00-abstract.tex:21–23`)
  - **Verdict: partly disagree.** The scope sentence can move earlier. The proposed ending claims
    experimental design, perturbation prediction and full-cell scaling, none of which the paper
    shows; the knockdown is the one perturbation and it is unresolved. The present ending is the
    paper's actual claim.
  - **Tier A, later.** Run `/write-abstract` once M10 and KIM-116 settle the last clause.

<a id="6d"></a>
- [ ] **6d. Replace the §2 code listing.**
  - **Review:** "Currently demonstrates releasing fixed constants for an appendix-only sensitivity
    check… Replace with a clean 6-line snippet of the main paper's workflow, or move the release
    snippet to `app_m0predict.tex`." ([L82](manuscript_review.md#L82))
  - **Checked:** accurate (`02-infercell.tex:225–238`), and already open in `STATUS.md`. The listing
    also shows the internal name `COREA_CYCLE_S`.
  - **Verdict: agree.**
  - **Tier A.** Move the release listing to `app_m0predict.tex`. Put `infer` then `predict`, as
    §5–§6 call them, in §2. Check the snippet runs.

## §4 Presentation and polish

<a id="f1"></a>
- [ ] **F1. Redraw `fig_architecture`.**
  - **Review:** "the `predict` box crowds the short vertical arrows emerging from `Mean field` and
    `Mixed route`." ([L89](manuscript_review.md#L89))
  - **Checked:** I looked at the figure. The four arrows into `predict` are a few points long and
    the lower pair nearly vanish.
  - **Verdict: agree.** **Tier A.** Use `/publication-figure` or `scientific-schematics`.

<a id="f2"></a>
- [ ] **F2. Redraw `fig_model` so arrows join modules.**
  - **Review:** "the four horizontal arrows sit evenly between columns rather than connecting to
    specific modules, forcing the caption to warn 'Arrows run between the blocks, not the boxes they
    touch.'… explicit 1 s / 60 s clocks." ([L90](manuscript_review.md#L90))
  - **Checked:** accurate on the arrows and the caption (`03-testsystem.tex:39`). The clocks are
    already on the arrow labels.
  - **Verdict: agree, apart from the clocks.** A caption that warns the reader about its figure is
    a sign the figure is wrong. **Tier A.** Draw block enclosures and run each arrow from its real
    source module, then delete the warning sentence.

<a id="m1"></a>
- [ ] **M1. Clear the unfinished markers.**
  - **Review:** "22 `\todo{}`, 2 `\placeholder{}`, 1 `\DRAFT{}`…, and 1 `[CITE NEEDED]`."
    ([L91–92](manuscript_review.md#L91))
  - **Checked:** counts are exact, plus one `\placeholderfig` for F7. By dependency: 13 todos wait
    on M10, 6 on KIM-116, and 3 on KIM-90, which can be filled now ([1a](#1a)). The `\DRAFT{?}` is
    the M0 mixed route's wall-clock days (`05-inference.tex:137`), with no source found. The
    placeholders are the Zenodo DOI and the thanks.
  - **Verdict: agree.** **Tier A** for KIM-90's three, the `\DRAFT` and the citation. The rest
    follow their runs. `check_paper --final` does not see `\todo`, so grep before submission.

<a id="h1"></a>
- [ ] **H1. Close the open items in `STATUS.md`.**
  - **Review:** "naming of Q vs. 'lactate exported'; clarifying that the prior state in Fig. 6 draws
    20 parameters rather than 21; reconciling gradient timing 37.9 s vs. 43–45 s."
    ([L94](manuscript_review.md#L94))
  - **Checked:** all three are in `STATUS.md`. `budget.ndrawn` is 20 and the text does say "all 20
    drawn parameters" (`06-prediction.tex:98`), but never says why it is not 21. The 37.9 s figure
    appears only inside a todo.
  - **Verdict: agree.** **Tier A.** One clause that σ_b is not drawn; one sentence that the two
    gradient costs are different measures on different nodes; your call on Q's name.

<a id="h2"></a>
- [ ] **H2. Give the production cost beside the whole-run cost.**
  - **Review:** "separate the production posterior cost from the exploratory 24.3 CPU-h unit-metric
    probe and PPC." ([L95](manuscript_review.md#L95))
  - **Checked:** the caption already itemises the probe (`05-inference.tex:122–128`). The pilot note
    has the split: the four chains took 452 CPU-h of the 492 (299 for the first 200 + 200, and
    153 for the extension).
  - **Verdict: mostly already done.** The headline figure overstates the cost, which is the safe
    side. **Tier A, small:** add "of which the chains are about 452" to the caption or the text.

<a id="h3"></a>
- [ ] **H3. Remove the dead label in `03-testsystem.tex`.**
  - **Review:** "Remove redundant `\label{sec:results}` / `\label{sec:model}` pair."
    ([L96](manuscript_review.md#L96))
  - **Checked:** `sec:results` (line 2) is referenced nowhere. `sec:model` (line 12) is used throughout
    and resolves to the same section.
  - **Verdict: agree.** **Tier A.** Delete `sec:results` and move `sec:model` up to the heading.

## §5 Roadmap

Steps 1–6 of the roadmap ([L104–114](manuscript_review.md#L104)) restate §3 and map to 1a, 1d and
2a, 2b, 3a and 4a, 6a and 4b and F1–F2, and 5b–5c. Two points are new.

<a id="r1"></a>
- [ ] **R1. The disjoint-cohort test of double counting on M0.**
  - **Review:** "splitting cells into disjoint transcript vs. bulk cohorts cleanly tests
    double-counting." ([L104](manuscript_review.md#L104))
  - **Checked:** the gap is mixed against exact, so a clean test needs an exact posterior on the
    disjoint design too: 508–678 CPU-h. The KIM-90 note records that it was not run, by your call
    of 2026-11-23.
  - **Verdict: agree it is the clean test; the decision is already made.** No action unless
    [2a](#2a) shows the reduced-model shift repeats, which would raise its value.

<a id="r2"></a>
- [ ] **R2. Order of work.**
  - **Review:** finish M10, KIM-116 and KIM-90 first, then diagnose, then the perturbation, then
    scaling, then the restructure. ([L104–113](manuscript_review.md#L104))
  - **Verdict: partly.** I would not start with M10. It is the most expensive item, and both its
    design and whether it runs at all depend on KIM-116. The order in "What I would do" above puts
    the venue and tier A first, then the two cheap experiments.

## §6 Strategy

<a id="a1"></a>
- [ ] **A1. Reframe around targeted and per-module inference.**
  - **Review:** "Positioning InferCell around compositional, targeted, and per-module inference,
    rather than monolithic joint inference over all parameters, is a significantly stronger, more
    original, and more defensible thesis. Currently, the manuscript accidentally argues the
    opposite." ([L122–133](manuscript_review.md#L122))
  - **Checked:**
    - The direction matches the project's own record. The framing note of 2026-10-26 dropped
      "posteriors for every parameter" for "goal-directed inference". Spec D7 says "per-module or
      super-module inference is the practical route on a full whole-cell model."
    - The joint posterior is not an accident. D7 keeps it as the reference: "the reduced model is
      the one place where a joint posterior is affordable, so it is where the cut's loss and saving
      can be measured."
    - "Superseded" (`app_m0predict.tex:12`) refers to the budget transferred from M0, not to
      modular inference.
    - "Not a joint posterior… a cut" is a precise statement of the estimand, not defensiveness. It
      stays.
    - The budget as an experimental-design tool has a snag: KIM-96 found no held set left on the
      reduced model, so the only budget with something to split is KIM-42's, now in the appendix.
  - **Verdict: agree on the direction, partly on the diagnosis, and it is conditional.** Everything
    turns on KIM-116. If the cut lands close to the joint at a fraction of the cost, Point A is the
    paper's best story and the joint becomes its yardstick. If it does not, the result is "the
    joint was needed here", and Point A would oversell.
  - **Tier C.** Wait for KIM-116's note. A full reframe changes the contribution sentence and pulls
    appendix material back into the main text, two days after the last restructure; it needs your
    sign-off and a Decision log entry. Meanwhile, two sentences are safe now: the Introduction and
    the Discussion can say the joint posterior is the reference, and targeted inference the
    practical route.

<a id="b1"></a>
- [ ] **B1. Drop the reduced model for M0 plus full syn3A.**
  - **Review:** "replacing the 3-tier hierarchy… with a 2-tier architecture (M0 + Full JCVI-syn3A)
    would make the manuscript substantially cleaner and much higher impact."
    ([L137–153](manuscript_review.md#L137))
  - **Checked:**
    - SpeedyCell.jl is a standalone simulator of the 2022 model, and InferCell does not depend on
      it. A grep of its `src` finds ForwardDiff only as the ODE solver's sparse Jacobian backend
      (`src/ODE/ODE.jl:330`), and no mean-field gene expression and no parameter gradient. The mixed
      route needs both. I grepped; I did not read the package. The decision of 2026-10-19 was "the
      full cell as a black box first".
    - Calibrating the full model already has a paper. PreCell's draft (26 pages, handoff dated
      today) is "Amortised Bayesian calibration of a whole-cell model with a diffusion
      transformer", on JCVI-syn3A, from bulk measurements including real data. It infers 28
      group multipliers, not individual parameters, so the overlap with Point B is in the organism
      and the real-data constraints, not in the method or the parameter set. This is a weaker
      objection than the first one.
    - It reverses the decision of 2026-11-09 that made the reduced model the main result, and it
      moves about 6,000 CPU-h of finished §5–§6 results out of the main text.
    - Joint-versus-modular "inside the full cell" has no reference unless a joint posterior is
      affordable there. That is the reason the reduced model exists.
    - It contradicts the review's own §2, which praises the declared reductions.
    - 93 reductions is accurate (`red.total`). "40,000+ CPU-h" is accurate for the exact long run
      and calibration (6,810 + 33,900). The claim that GTP clipping is an artefact of the reduction
      is unverified. The clip rule is the published model's, and no one has counted clips in the
      full cell.
  - **Verdict: disagree.** The review's own third reality check has two branches, and the one that
    applies is "months of core solver engineering: finish the reduced-model manuscript… and reserve
    full-cell compositional inference for a follow-up". So the review reaches this conclusion too,
    given the facts above.
  - **What is right in it:** an editor may see a toy. The answer is a clear sentence on what the
    reduced model is for (the one place a joint reference is affordable), and the Full-cell stage 1
    milestone, not a pivot.
  - **Tier C.** Confirm the reduced model stays, and I will record nothing further.

<a id="rc1"></a>
- [ ] **RC1. Can the full cell take a mean-field swap, or only module-local inference?**
  - **Review:** "Discrete full-cell processes… do not close as linear first-order ODEs. Verify
    whether the full cell supports a smooth mean-field swap… or whether InferCell should perform
    module-local inference." ([L159–161](manuscript_review.md#L159))
  - **Verdict: agree it is the right question, for paper 2.** Today's answer is that neither
    exists: the full cell is not composed in InferCell. KIM-116's data plug-in, which needs no ODE
    solve, is the first evidence on the module-local route. Not checked further.
  - **Tier C.** Note it on the Full-cell milestone when KIM-116 reports.

<a id="rc2"></a>
- [ ] **RC2. Per-gradient cost on the full cell.**
  - **Review:** "If full-cell gradients run in ≤60–90 s…, full-cell targeted NUTS is tractable; if
    gradients take 10–30 min, module-local inference + full-cell composition is the clear path."
    ([L162–163](manuscript_review.md#L162))
  - **Checked:** unmeasured, and unmeasurable until RC1 is answered. For scale, the reduced model's
    gradient was 277 s before KIM-99 and 38–45 s after, at about 8% of the cell.
  - **Verdict: agree as a gate for paper 2.** No action for paper 1.

<a id="rc3"></a>
- [ ] **RC3. Timeline and publication strategy.**
  - **Review:** "If… integrates cleanly in weeks: pivoting… is a clear upgrade… If full-cell
    integration requires months…: finish the reduced-model manuscript for *Nature Computational
    Science* (adopting the compositional framing from Point A plus the diagnostic fixes)."
    ([L164–166](manuscript_review.md#L164))
  - **Verdict: agree with the second branch.** It is months, for the reasons in [B1](#b1). That
    leaves the plan this checklist already describes: tier A, the two cheap experiments, then Point
    A and M10 decided on KIM-116.
  - **Tier C.** Settled by B1 and V1.
