# Regular-geometry screen for timed previews

This development smoke check uses the exposed COD development input in
`../2026-11-27-selenium-crystal-stereo/source-development-input.json`. Selection
is deterministic and outcome-independent: ascending source index, first 30
SMILES containing `@` whose parsed molecule has at most 22 heavy atoms. There
were no selection errors. All molecules were parsed and had explicit hydrogens
added before using the timed worker's prepared corrected embedder once with
seed 42 and its default iteration setting.

All 30 embedding attempts returned coordinates. On these same coordinates,
the existing preview stereo check and the new opt-in regular-geometry check
both returned true for all 30; each was cross-checked against its corresponding
public `stereo_3d_diagnostics` result. Exact coordinates, atom-center reports,
selection and all outcomes are in `results.json`. The selected collection is
small and is not a molecular holdout, a representative rate, an independent
stereo oracle, or proof of conformer quality. It does not measure timeout or
minimization behavior. A separate distorted-coordinate regression fixture
exercises the important difference: the ordinary preview check accepts it and
the regular-geometry preview check rejects it.

The local macOS arm64 debug wheel used CPython 3.13.6, was built from commit
`baada70f` plus the accompanying source changes with `maturin build --locked`,
and has SHA-256
`51ee2466457abbdfc350446c09abe714b70f25af7cc2dbdbf5aabef13b0b78c2`. It is not
a published wheel. Reproduce with that candidate installed:

```sh
python validation/results/2026-11-28-timed-preview-regular/run.py > output.json
```

Worker-interruption, checkpoint metadata validation, ordinary-mode
compatibility, and a completed-ensemble control are tested separately in
`crates/chematic-py/tests/test_timed_conformers.py`. The option affects only
provisional previews: it does not add the regular-geometry screen to completed
conformers, energies, or force-field optimization.
