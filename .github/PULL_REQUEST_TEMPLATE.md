<!--
Thanks for contributing to InferCell. Please fill in the sections below.
Delete anything that doesn't apply.
-->

## Summary

<!-- 1-2 lines: what changes, and the *why* (not the *what* — the diff shows that). -->

## Linked issue

Closes #

## Checklist

- [ ] Unit tests pass locally: `julia --project -e 'using Pkg; Pkg.test()'`
- [ ] If this PR touches inference / boundary protocols, the integration suite passes on Slurm: `sbatch test/run_integration.slurm` — quote the job id below. GitHub does not run these on the merge path, so the author is the gate.
- [ ] Docs updated if behaviour or public API changed
- [ ] `Manifest.toml` changes (if any) were regenerated via `Pkg`, not hand-edited
- [ ] If scientific results change (CIs, posterior summaries, calibration), the PR description quantifies the change

## Notes for reviewers

<!-- Anything reviewers should focus on, alternatives you considered, or things you're unsure about. -->
