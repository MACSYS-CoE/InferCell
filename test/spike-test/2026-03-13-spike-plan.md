InferCell's core bet is that Bayesian inference can flow through an ODE solver in pure Julia — NUTS sampling parameters, AD computing gradients through the solve, posteriors recovering known ground truth. 

There could be some potential friction between Turing.jl + DifferentialEquations.jl, e.g. 

- **AD backend compatibility.** Not all solver/AD combinations work. In-place dynamics functions fail with some backends. Solver choice matters.
- **Sensitivity algorithm selection.** `ForwardDiffSensitivity()`, `InterpolatingAdjoint()`, and `BacksolveAdjoint()` have different trade-offs and failure modes. The right choice depends on problem size and structure.
- **Numerical stability under sampling.** NUTS proposes parameter values across a wide range. The ODE solve must not fail, produce NaNs, or silently return garbage for extreme parameter proposals. Solver tolerances and error handling matter.
- **Performance.** If a single NUTS sample takes seconds (due to recompilation, allocations, or slow sensitivity computation), inference on even a toy model becomes impractical.

A spike that validates all of this in ~50 lines — before any abstraction layer exists — saves us from building architecture on top of a broken foundation.

## What the spike proves

1. **Forward solve works.** TX/TL ODEs produce plausible mRNA/protein trajectories.
2. **AD flows through the solver.** Turing.jl + NUTS can differentiate through `solve()` with `Tsit5()` and `ForwardDiffSensitivity()`.
3. **Parameter recovery.** Posteriors concentrate around true values used to generate synthetic data.
4. **Observation noise recovery.** Inferred `σ_obs` matches the noise level used to generate data.
5. **Reasonable performance.** Sampling completes in minutes, not hours.

If all five hold, we proceed to implement the full architecture with confidence. If any fail, we learn exactly what needs to change before investing in abstractions.

## Implementation plan

A single self-contained script: `spike/spike_txl_inference.jl`.

### Step 1: Define the ODE (out-of-place)

```julia
function txl_dynamics(u, p, t)
    mRNA, protein = u
    k_tx, k_tl, γ_mRNA, γ_protein = p
    SA[k_tx - γ_mRNA * mRNA, k_tl * mRNA - γ_protein * protein]
end
```

Using `StaticArrays.SA` for the return value — small fixed-size, allocation-free.

### Step 2: Generate synthetic data

- True parameters: `k_tx=1.0, k_tl=2.0, γ_mRNA=0.5, γ_protein=0.1`
- Initial conditions: `mRNA₀=0.0, protein₀=0.0`
- Timespan: `(0.0, 50.0)`, observe at `t = 0:2.5:50`
- Additive Gaussian noise: `σ_true = 0.3`
- Solve with `Tsit5()`

### Step 3: Define the Turing model

```julia
@model function txl_inference(data, times, prob)
    # Priors
    k_tx ~ LogNormal(0, 1)
    k_tl ~ LogNormal(0, 1)
    γ_mRNA ~ LogNormal(0, 1)
    γ_protein ~ LogNormal(0, 1)
    σ ~ HalfNormal(1.0)

    # Solve ODE
    p = SA[k_tx, k_tl, γ_mRNA, γ_protein]
    sol = solve(remake(prob, p=p), Tsit5(); saveat=times,
                sensealg=ForwardDiffSensitivity())

    # Bail if solver failed
    if sol.retcode !== ReturnCode.Success
        Turing.@addlogprob! -Inf
        return
    end

    # Likelihood (per-timepoint)
    for i in eachindex(times)
        data[:, i] ~ MvNormal(sol[:, i], σ)
    end
end
```

### Step 4: Run NUTS

```julia
model = txl_inference(data, times, prob)
chain = sample(model, NUTS(), 1000)
```

### Step 5: Validate

- Print summary statistics — do posterior means/medians bracket the true values?
- Check that `σ` posterior concentrates near `σ_true = 0.3`.
- Print elapsed time.
- If available, plot posterior predictive trajectories against data (optional, not required for the spike to pass).

## Success criteria

| Criterion | Pass condition |
|-----------|---------------|
| Sampling completes | No errors, no crashes |
| `k_tx` recovery | True value (1.0) within 90% credible interval |
| `k_tl` recovery | True value (2.0) within 90% credible interval |
| `γ_mRNA` recovery | True value (0.5) within 90% credible interval |
| `γ_protein` recovery | True value (0.1) within 90% credible interval |
| `σ_obs` recovery | True value (0.3) within 90% credible interval |
| Performance | 1000 samples in under 10 minutes |


