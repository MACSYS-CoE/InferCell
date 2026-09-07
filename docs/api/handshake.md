# The 1 s handshake and the 60 s rebuild

Source: [`src/handshake.jl`](https://github.com/MACSYS-CoE/InferCell/blob/main/src/handshake.jl).

A mixed ODE/jump composition is advanced by a **split-operator exchange**, not by a single solver call. Every simulated second, counts from the stochastic block become rate-law parameters of the metabolic block, the metabolic block integrates one second, and the costs the stochastic block accrued are debited back against its pools. This mirrors the published whole-cell model's `hookSimulation`, whose communication interval is `delt = 1.0` s.

The two blocks keep **separate state vectors** — the ODE block in mM over `Float64`, the jump block in particles over `Int` — which is what makes the conversion below a real boundary rather than a units comment.

## `HandshakeDriver`

Returned by [`build_problem`](orchestrator.md) for a mixed composition. It holds two live integrators and the exchange records lowered from the declared [coupling edges](corea-interface.md).

```julia
driver = build_problem([metabolism, expression]; tspan = (0.0, 600.0))
record = run_handshake!(driver, 600)
record.census        # (; handshakes, drains, clipped, fraction, deficits, pending)
record.rebuilds      # one row per rebuilding module: params, pools, interval, refreshes
record.growth        # per handshake: area_nm2, radius_nm, volume_litres, factor, fractional, capped
```

`run_handshake!` refuses to run past the `tspan` the driver was built with — `n_steps * interval` must fit inside it — and refuses a run that would not end on a drain, since that would leave accrued cost in the counters that nothing debits and nothing names. It returns `(; t, ode, jump, jump_p, growth, census, rebuilds)`: the handshake times, the ODE state in mM, the jump state in particles, the jump block's parameter vector and the cell's geometry at each of them — per handshake rather than per solver step, because the handshake is the only instant at which the two blocks agree on a state. `jump_p` and `growth` are what make the two channels with no state of their own observable from outside: a refresh count is read off `jump_p` and the growth law off `growth`, rather than restated from the schedule or the counts that produced them. `handshake_step!` runs one exchange.

### Policy keywords

| keyword | default | what it fixes |
|---|---|---|
| `interval` | `1.0` | the exchange period, in simulated seconds |
| `drain_interval` | `interval` | how often the deferred counters are debited; must be a whole number of exchanges, and is a labelled reduction when coarser |
| `rounding` | `:fractional_carry` | how a continuous pool is written back as whole particles |
| `radius_nm` | `COREA_INITIAL_RADIUS_NM` (200.0) | the cell radius behind the count↔concentration conversion, for a composition with a **fixed** cell; refused alongside a growing one |
| `initial_surface_area_nm2` | `COREA_INITIAL_SURFACE_AREA_NM2` (502831.0) | the published initial membrane area, the primitive for a **growing** cell; refused on a fixed one |
| `footprint_nm2` | `MEMBRANE_PROTEIN_FOOTPRINT_NM2` (28.0) | nm² of surface area per membrane protein; refused on a fixed cell |

The three geometry keywords all default to `nothing` in the signature and resolve to the constants above, so that a value passed to the composition it does not apply to can be refused rather than accepted and never read.
| `ode_solver` | `Rodas5P()` | the stiff integrator for the metabolic block |
| `abstol`, `reltol` | `1e-10`, `1e-8` | pinned tolerances |

## Counts and concentrations

`corea_particles_per_mM(radius_nm)` derives the conversion factor from Avogadro's constant and the volume of a sphere, rather than carrying a transcribed constant: at the published 200 nm radius it returns 20180.39 particles per mM, the 20,180 the scoping note records. `counts_to_mM(n, factor)` is exact division; the other direction needs a policy, because particles are whole.

The factor is **live** where the composition grows. A growing cell starts from the published *area* rather than a radius, so its initial factor is 20190.998 rather than 20180.39 — a 0.05% difference, and the reason the driver refuses to be handed both.

### Rounding policies

Writing a pool back as whole particles once per handshake, ~6,300 times per cell cycle, gives the three policies three different signatures:

- `:fractional_carry` — the remainder is carried to the next handshake, so the emitted counts track the exact running total and the residual **does not accumulate**. The default.
- `:stochastic` — round up with probability equal to the fraction. Unbiased, so the residual accumulates as the **square root** of the handshake count.
- `:deterministic` — round to nearest, discarding the remainder. Biased, so the residual accumulates **linearly**, reaching up to 0.156 mM per species over a cycle — about 140× the integrator's own bound, and measured at 0.125 mM on the phase-3 toy. Implemented only so that this signature can be demonstrated; it is not a policy this project runs under.

The policy is a field on the driver rather than an argument to the arithmetic, so a trajectory carries the policy that produced it.

## Deferred debits

A [`DeferredCounterEdge`](corea-interface.md) declares that the stochastic block accrues a cost in a named counter and the hook debits it against a pool one step later. Under the published `:clamped_deficit_carried` policy the pool floors at zero and the shortfall is carried to the next handshake — the interface's `max(0, ·)`, and the reason a clamped counter is reported as a gradient obstruction. The two labelled departures also execute: `:unclamped` lets the pool go negative, and `:smoothed` replaces the floor with a softplus of the declared width, which converges to the clamped answer as that width narrows.

One counter may feed several pools — the published charged-tRNA transfer debits the charged pool and credits the uncharged one from a single accrual. The counter is read and cleared once per handshake, so every pool it feeds sees the same accrued value, and a credit matches what the debited pool actually paid rather than what was asked of it.

`clipping_census(driver)` reports how many handshakes ran, how many of them applied a debit, and at how many of *those* a counter carried a shortfall material against its cost. **The fraction is per drain, not per handshake** — only a drain can clip, so at a coarse `drain_interval` a per-handshake denominator would dilute it by exactly `steps_per_drain`. `pending` reports the accrual still sitting in the counters, which between drains belongs to no `deficit`. The comparison is relative rather than against zero: under `:smoothed` the softplus leaves a strictly positive residue at every step however deep the pool, so an absolute test would report every handshake as clipping. A non-zero count at published parameters means the non-smooth drain is in the operating regime.

## The 60 s rebuild

A [`RateConstantEdge`](corea-interface.md) declares that a pool re-enters the stochastic block as a recomputed rate constant. A jump module names the parameters the rebuild fills with `rebuilt_params(m)` and computes them in `rate_constants(p, t, m, pools)`, receiving the pools its inbound edges name, in mM. The edge carries no parameter slot and one pool feeds many constants, so the edges name the inputs and the module names the outputs.

The values go into the jump block's **parameter vector**, written by the hook after the debit and before the SSA step — so `remake(prob; p = θ)` and `model_free_params` still see them, where a constant held on the sub-model struct would be invisible to both. A rebuild fires at the handshake whose end time is a multiple of the edge's declared interval, and an interval that is not a whole number of handshakes is refused rather than rounded. A `cadence = :continuous` edge is refused too: the outer loop holds a constant between refreshes by construction, so declare a shorter piecewise-constant interval instead.

Note that a rebuilt parameter is a *derived* quantity that the inference entry points do not yet know is derived — they sample every `model_free_params` entry, and the hook overwrites the draw at the first refresh.

`rate_constant_elasticity(driver)` measures the channel's gain, `d ln k / d ln pool`, at the driver's current pools, by a central difference in log space. It carries the pool concentration with the number, because an elasticity is a function of the pools and not a constant of the model.

`rebuild_census(driver)` reports one row per rebuilding module: its parameters, its pools, its interval and how many times it has refreshed.

## Growth and volume

A [`VolumeEdge`](corea-interface.md) declares that a state's count sets the cell's membrane surface area. The module owning that state flags it with `membrane_protein_states(m)` and declares one outbound edge per flagged state; `membrane_protein_states(models)` is the sweep that finds them all, so a composition is asked rather than a module type named. Either block may declare it — a jump state is already a count, an ODE state's count is `u * factor` at the factor in force at that hook, which is the published semantics.

The chain is the published model's, from `in_out.py:107`:

```
area   = baseline + 28.0 nm² × Σ(membrane protein counts)
radius = sqrt(area / 4π)
volume = (4/3) π radius³,   capped at twice the initial volume
```

The baseline is derived at build time, not written down: it is whatever makes the composition's own initial counts come to `initial_surface_area_nm2`. For Core A′ that is 502,831 − 831 × 28.0 = 479,563 nm², and it is **frozen** — the lipid leaflet and the membrane-protein loci outside the reduction are exogenous, which is why growth reaches ~1.07× rather than approaching the cap.

It runs at **step 0 of every handshake**, on neither the 1 s nor the 60 s clock, so every conversion in that exchange uses one volume. A change in the factor rescales every ODE state at constant particle count, which is what makes growth a coupling channel rather than a reported diagnostic. A state whose millimolar is referred to another volume — external lactate, in the medium — is exempted with `extracellular_states(m)` and is not diluted.

`growth_census(driver)` reports the flagged states and which block each is read from, the membrane count, the baseline and footprint, and the live area, radius, volume, factor and fractional growth. `run_handshake!` records the same geometry per handshake, so the law is read off the trajectory rather than asserted from the counts that produced it.

**Growth may be reported as fractional growth or time-to-threshold, and as nothing else.** `reporting_constraints(driver)` returns that constraint as retrievable data and `doubling_time(driver)` throws: Core A′ removes roughly 92% of the cell, so a comparison to syn3A's 105-minute doubling time would be meaningless, and spec §7 requires the code rather than the prose to refuse it. `growth_report(driver)` and `time_to_threshold(record, x)` are the two admissible reportings.

## Driver-level departures

`driver_declarations(driver)` enumerates the departures carried by the driver's *policy* rather than by any module's declarations — a `drain_interval` coarser than the exchange, a rounding policy other than fractional carry, and, where the cell grows, the 28.0 nm² footprint being the published model's *calibrated* constant (its own docstring says 35 nm²; the code governs) and the frozen exogenous membrane baseline. `reduction_declarations(models)` structurally cannot see any of them, so pass the driver as a second argument to `reduction_declarations` or `reduction_report` when reporting beside a result.

## What this layer does not yet do

Volume re-entering an ODE rate law — the other half of the volume channel, which spec §11 task 7.2 wants for the lactate exporter's `3P/r` — is not built: `VolumeEdge` carries no parameter slot to write it into, and an inbound edge is refused rather than left to resolve and never run. The clamped edge's held value still travels as a fixed parameter. Nor does the driver plug into the inference entry points, which build a single SciML problem and call `remake` on it.
