# Contributing to InferCell

Thanks for your interest in contributing. This guide covers the local
development workflow, the checks CI runs on every PR, and the repo-level
settings that keep `main` healthy.

## Local setup

```bash
git clone <repo-url>
cd InferCell
julia --project -e 'using Pkg; Pkg.instantiate()'
```

`Manifest.toml` is committed, so `instantiate` reproduces the exact
dependency versions CI uses.

## Running tests

```bash
# Unit tests only (~1 min, runs on every PR)
julia --project -e 'using Pkg; Pkg.test()'

# Unit + integration tests (~3 min, also runs on every PR)
INFERCELL_INTEGRATION_TESTS=true julia --project -e 'using Pkg; Pkg.test()'
```

Aqua quality checks (stale deps, missing compat bounds, unbound type
parameters, undefined exports) run as part of the test suite — they need no
extra invocation.

## Formatting

InferCell uses [JuliaFormatter](https://domluna.github.io/JuliaFormatter.jl/)
with the [SciML style](https://github.com/SciML/SciMLStyle), configured in
`.JuliaFormatter.toml`. To format the whole tree:

```bash
julia --project -e 'using Pkg; Pkg.add("JuliaFormatter"); using JuliaFormatter; format(".")'
```

The `Format` workflow runs the same check on every PR. It is currently
**advisory** (it surfaces diffs as a warning but does not block merges).
After the codebase has been baseline-formatted, it will be flipped to a
required check.

## Branch workflow

- Branch off `main` for every change. Don't commit to `main` directly.
- Branch naming convention: `issue-<n>-short-slug` (e.g. `issue-14-ci-cd-setup`).
- Open a PR back into `main`. Squash-merge is the default.
- Link the issue in the PR description (`Closes #<n>`).

## Commit messages

Short imperative subject (50 chars), then a body that explains the *why*.
References to issues live in the body or PR description. Example:

```
Add Aqua quality checks to the test suite

Closes #14. Catches stale compat entries, missing compat bounds, and
undefined exports as part of `Pkg.test()`. `ambiguities` and `piracies`
are disabled for now (SciML/Turing extensions produce false positives);
re-enable after an audit.
```

## CI checks on every PR

| Check                  | Blocking? | Notes                                              |
| ---------------------- | --------- | -------------------------------------------------- |
| `Test`                 | Yes       | Unit + integration tests, Aqua, coverage summary.  |
| `Format`               | No (yet)  | JuliaFormatter diff; advisory until Stage 2.       |
| `Documentation/build`  | Yes       | mkdocs strict build of the `docs/` site.           |

Coverage is reported on every CI run (job summary line and `lcov.info`
artifact) but is **not** gated — there is no minimum coverage threshold.

## Repo-level settings (admin checklist)

These are set in the GitHub UI, not via PRs. Apply once per repo:

1. **Settings → Branches → Add rule for `main`**:
   - Require a pull request before merging
   - Require status checks to pass: `Test`, `Documentation/build`
     (add `Format` after Stage 2)
   - Require linear history
   - Require conversation resolution before merging
   - Do not allow administrators to bypass (once the team has grown)

2. **Settings → Actions → General**:
   - Allow GitHub Actions
   - Workflow permissions: read-only (workflows that need writes opt in
     explicitly via `permissions:` blocks)

3. **Settings → Pages** (when ready to publish docs publicly):
   - Source: GitHub Actions
   - Set the `DOCS_DEPLOY_ENABLED` repository variable to `true`
     (Settings → Variables → Actions).
