# The cut predictive distribution and p(Q | D, do(a))

**Date:** 2026-10-29
**Status:** a methods note for KIM-20. It fixes what `predict` (KIM-27) computes
and seeds paper 1's Methods. It decides nothing that spec §12 has not already
recorded. The choices it leaves open are named, and are fixed in M4's freeze entry
(KIM-39).

**Based on:**
- the completion plan (`dev/notes/2026-10-29-completion-plan.md`) §3.1–3.4 and §4;
- the launch-readiness note
  (`dev/notes/2026-10-29-launch-readiness-and-predictive-demonstration.md`) §5–§8;
- spec §12 2026-10-29, "v0.3.0 launch decisions".

## 1. Notation

- θ_I: the inferred dynamical parameters, M0's five targets (`M0_TARGETS`).
- σ_b: the bulk observation noise. It is inferred beside θ_I, but it is not a cell
  parameter and never enters a forward simulation.
- θ_H: the held parameters, the rest of the reduced model's informed constants.
  They are completed from priors, not inferred.
- ξ: intrinsic noise, meaning one cell's stochastic trajectory, fixed by its seed.
- D: the data, M0's cell and bulk record.
- V: the viability event, a declared rule on simulated summaries (KIM-36).
- Q(θ, ξ): the predicted quantity for one cell. O4 is provisional: lactate
  export, as `carbon_accounts`' `lactate_out` or `exported`.
- a: an intervention. g_a(θ) is θ with the intervention applied.

## 2. The aim and the estimand

The aim is the predictive distribution of Q given the data:

```text
p(Q | D) = ∫ p(Q | θ) p(θ | D) dθ,     p(Q | θ) = ∫ δ(Q − Q(θ, ξ)) p(ξ) dξ.
```

The estimand is the population mean (spec §12 2026-10-29, O4):

```text
μ(θ) = E_ξ[ Q(θ, ξ) ],
```

and the prediction is the distribution of μ(θ) under the parameter distribution.
The intrinsic spread σ²(θ) = Var_ξ[Q(θ, ξ) | θ] is a separate output. Only the
hybrid forward model has it; the mean field has no single-cell distribution.

The hybrid estimates μ(θ) from n_cells cells per draw:

```text
μ̂(θ) = (1/n) Σ_c Q(θ, ξ_c),     Var[μ̂(θ) | θ] = σ²(θ) / n.
```

So the spread of μ̂ over the parameter draws overstates the spread of μ by
E_θ[σ²(θ)] / n. `predict` reports that Monte Carlo term beside the result. The
parameter variance in §5 subtracts it, and does not leave it inside.

## 3. The cut distribution

### 3.1 On M0: no held parameters

M0's 50 datasets come from M0 with every non-target parameter at its published
value, and the fit holds them there. So θ_H is empty and

```text
p(μ | D) = ∫ δ(μ − μ(θ_I)) p(θ_I | D) dθ_I,
```

with σ_b marginalised out. This is the **scored** case. The truth θ_I* is stored
beside each archived posterior (`dev/notes/2026-10-29-archived-posteriors.md`).
So μ(θ_I*), estimated with many more cells than any posterior draw uses, is the
value whose coverage is measured.

On M0, a viability filter can act only on θ_I. That would be extra data, and it
would change the posterior being scored (completion plan §3.3). So on M0 there are
two options:
- apply no filter, which is the default. For scale: 4% of prior draws over D11's
  six targets starve on Core A′ (14c.2). That is a Core A′ prior figure, not M0's
  posterior rate, which has not been measured;
- condition jointly, as in §3.3 below, and report it as a different, conditioned
  posterior.

### 3.2 On the reduced model: the modular transfer

θ_I's posterior comes from M0. θ_H gets no data, only its prior and V. The cut
distribution in the factorised form is

```text
p_cut(θ_I, θ_H | D, V) = p(θ_I | D, M0) · p(θ_H | V).
```

It is a deliberate modular transfer, **not** the joint posterior under the reduced
model. Three things are cut:
1. **D does not update θ_H.** The M0 likelihood does not involve θ_H.
2. **θ_H does not feed back into θ_I.** M0 fixed the counterparts of θ_H at their
   published values. If θ_I's posterior would move when they vary, the transfer
   hides it. This is the dependence assumption to state. **It is untested in
   v0.3.0.** The prior-width runs (½, 1, 2) do not test it: they change how far
   θ_H spreads the prediction, but θ_I's posterior is the same in every run. A
   test would refit M0 with some held counterparts moved off their published
   values and compare θ_I's posteriors. That is outside v0.3.0's budget, and the
   paper says so.
3. **V does not act as data on θ_I** in the factorised form.

The reduced model is **unscored**. No truth matches these data under the reduced
model, so this case shows the method's scale and its uncertainty budget. It is not
a calibration result.

### 3.3 When V depends on both groups

V is read from simulated cells, so it depends on θ_I as well as θ_H. The
factorised form is only exact if it does not. There are two well-defined variants:

```text
(a) joint:        p_V(θ_I, θ_H) ∝ p(θ_I | D) · p(θ_H) · P(V | θ_I, θ_H)
(b) per-θ_I cut:  p_V(θ_I, θ_H) = p(θ_I | D) · p(θ_H | θ_I, V),
                  p(θ_H | θ_I, V) ∝ p(θ_H) · P(V | θ_I, θ_H)
```

- **(a)** reweights θ_I by P(V | θ_I) = ∫ P(V | θ_I, θ_H) p(θ_H) dθ_H. To sample
  it, draw a pair and accept it if viable.
- **(b)** keeps θ_I's marginal equal to the M0 posterior, which is what "modular"
  means. For each θ_I draw, resample θ_H until a viable one is found. Record each
  θ_I's acceptance rate, which is an estimate of P(V | θ_I).

Both are reported with:
- the acceptance rate;
- the ESS;
- the correlations between θ_I and θ_H that the filter induces;
- the rejected draws, each with its reason.

"Clipped", "starved" and "viable" stay separate categories. Which variant is used
is fixed with V in M4 (KIM-36, KIM-39).

P(V | θ) is itself a probability over ξ: one parameter draw can give a viable cell
and a starved one. So V has to say whether it is a per-cell event or a property of
θ. A property of θ could be, for example, the fraction of n_cells that are viable
being at least a threshold. A per-cell V conditions the cells, not the parameters,
and changes μ's definition.

If no starvation is part of V, starvation cannot also be Q (note §6).

## 4. The intervention, p(Q | D, do(a))

The intervention is a ptsG promoter knockdown at dose d. It scales `S_0779` (ptsG)
by (1 − d) after the parameter draw:

```text
p(Q | D, do(a)) = ∫ p(Q | g_a(θ)) p(θ | D) dθ.
```

The reported effect is paired, with common cell seeds:

```text
Δ(θ, ξ) = Q(g_a(θ), ξ) − Q(θ, ξ),     Δμ(θ) = E_ξ[ Δ(θ, ξ) ].
```

Common seeds remove most of ξ's variance from Δ. The `predict` test for this
(KIM-32) is that d = 0 gives Δ ≡ 0 exactly.

**Assumptions** (stated in the paper's Methods):
1. **The intervention acts only through the named parameter.** g_a changes
   `S_0779` and nothing else.
2. **The posterior carries over.** It is fitted on unperturbed data, and p(θ | D)
   is used unchanged under do(a). No parameter is assumed to shift with the
   perturbation.
3. **The model stays valid at the dose.** The one-way nucleotide reduction rests on
   elasticities of 0.04–0.09. The dose screen (KIM-37) keeps the pools away from
   K_d. Its criterion concerns only whether the regime is valid, never the effect's
   size or coverage.
4. **V is read at baseline.** On the reduced model, θ is conditioned on baseline
   viability, and the intervened cell may then fail. That failure is part of the
   prediction, not a reason to reject the draw. This is proposed here, and fixed in
   the freeze.

On M0 the true effect Δμ(θ_I*) is computable with the same seeds. So coverage of
the true effect is scored over the 50 datasets, for exact vs mixed posterior and
for hybrid vs mean-field forward.

## 5. The variance split

### 5.1 Parameters against intrinsic noise

For one cell's Q, the law of total variance is

```text
Var(Q | D, V) = Var_θ[ μ(θ) ] + E_θ[ σ²(θ) ].
```

The two terms are parameter uncertainty and intrinsic noise. For the
population-mean estimand, only the first term is the estimand's uncertainty. The
second is the intrinsic spread, which is reported separately (§2).

### 5.2 Inferred against held, in both orders

Viability conditioning correlates θ_I and θ_H, so an independent ANOVA split is
not valid. The split is nested, and it is reported in both orders:

```text
I first:  Var_θ[μ] = Var_{θ_I}[ E_{θ_H | θ_I}[μ] ] + E_{θ_I}[ Var_{θ_H | θ_I}[μ] ]
H first:  Var_θ[μ] = Var_{θ_H}[ E_{θ_I | θ_H}[μ] ] + E_{θ_H}[ Var_{θ_I | θ_H}[μ] ]
```

- With independent groups, as in variant (b) when acceptance does not depend on
  θ_I, the I-first first term is θ_I's first-order effect. The H-first second term
  is θ_I's total effect.
- The gap between the two orders measures the interaction, plus the dependence the
  filter induced.
- Shapley effects are not used (completion plan §3.4).
- **How both orders are computed.** Simulate a crossed matrix: M posterior draws
  θ_I^m × K prior draws θ_H^k, with μ̂_mk from n_cells cells and the viability
  indicator v_mk for each pair. The cost is M·K·n_cells cells, and KIM-38 costs it.
  Each variant is then a set of cell weights on that one matrix:
  - (b): w_mk = v_mk / (M · Σ_k v_mk);
  - (a): w_mk ∝ v_mk.

  Both orders are weighted two-way splits of the same matrix, using conditional
  means along rows (θ_I fixed) and along columns (θ_H fixed). Under (b), θ_I given
  θ_H is not something one can sample directly, because it needs P(V | θ_I) to
  normalise. The matrix sidesteps that, since the row sums of v estimate P(V | θ_I).
  - **Bias:** the variance of the row means over m overstates
    Var_{θ_I}[E_{θ_H|θ_I} μ] by about E[within-row variance] / K_eff. Each μ̂_mk adds
    σ²/n on top. Both are subtracted with the usual unbiased two-way estimator,
    and an estimate that comes out negative is reported as zero, flagged.
  - **Fallback:** if M·K·n_cells is too dear, report the I-first order only, and
    say so. The H-first order is not then claimed.

The sentence to aim for (plan §4.3): *x% of this prediction's uncertainty comes
from parameters no measurement supports.*

### 5.3 Measurement noise

σ_b is observation noise. It is never pushed into a latent biological quantity.
It enters only when the target is a future noisy **bulk** observation. The
model's bulk noise is multiplicative lognormal (`observe_bulk`): a reading of the
mean over N cells, x̄, is Y = max(x̄, floor) · exp(σ_b ε), with ε ~ N(0, 1). So
the split is on the log scale. Away from the floor, it is

```text
Var(log Y | D, V) = Var_θ[ E(log x̄ | θ) ]  +  E_θ[ Var(log x̄ | θ) ]  +  E[σ_b²]
                    parameters             intrinsic, over N cells     measurement
```

By the delta method, Var(log x̄ | θ) ≈ σ²(θ) / (N μ(θ)²). Because a bulk reading
averages N cells, the intrinsic term shrinks as 1/N. It is not the single-cell
E[σ²]. The parameter term is split as in §5.2. Otherwise, for the
population-mean estimand, σ_b is dropped once the posterior is read in.

### 5.4 Model form

Model form is not a variance term. The formalism swap is a bias, measured on M0 by
the crossed design's columns (spec §12 2026-10-29, KIM-16):

```text
E_θ[ μ_hybrid(θ) ] − E_θ[ m(θ) ]          (uncorrected)
E_θ[ μ_hybrid(θ) ] − E_θ[ m(θ) + δ ]      (corrected)
```

- m(θ) is the mean-field Q.
- δ = μ̂_hybrid(θ_ref) − m(θ_ref), from a hybrid ensemble at θ_ref that is computed
  once.
- The corrected form assumes δ is constant over the posterior's support. The
  crossed design tests that assumption on M0.
- The mixed route's own δ covers only the 17 panel pools. It does not cover Q
  (completion plan §3.1), so the δ for Q is a new ensemble.

Every other reduction is declared, not bounded (the model-form table, KIM-47).

## 6. For `predict` (KIM-27 to KIM-34)

| Object here | `predict` | Issue |
|---|---|---|
| p(θ_I, σ_b \| D): V7 replicate files or B6 `Chains` | `posterior` | KIM-28 |
| σ_b carried, not simulated | a separate field of the draw type | KIM-28 |
| θ_H ~ p(θ_H), Haldane-consistent, width scale ½/1/2 | `completion = :prior` | KIM-29 |
| p(θ_H \| θ_I, V) or the joint variant | `completion = :viability_conditioned_prior`, V passed as data | KIM-30 |
| μ̂(θ) from n_cells, or m(θ), or m(θ) + δ | `forward = :hybrid \| :mean_field`, `delta = :none \| :ensemble` | KIM-31 |
| g_a, Δ with common seeds | `intervention` | KIM-32 |
| seeds, sources, V's decisions, failures | the result's provenance; the rng is forked and the global state restored | KIM-33 |
| μ(θ), σ²(θ) and the E[σ²]/n term | `estimand = :population_mean`; the result's fields | KIM-27 |
