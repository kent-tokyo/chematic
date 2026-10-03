# A6 MMFF94 source-wheel termination and adoption boundary

This is a **local source candidate**, not the published v1.0.31 wheel. It
contains the post-v1.0.31 compact-ring typing change as well as the diagnostic
termination field; the effects cannot be attributed to the new field alone.
The macOS arm64/Python 3.13 dev-profile wheel SHA-256 is
`e662bda42a6360c292946711e4b0a02f2e4c23ed09e1ec55faa2ed87efff1c76`.
Base revision: `91ccdcb04c4b538e90c4da667f76b2d289629975`;
the wheel-relevant source diff SHA-256 at build time is
`3aa9cb9c89ebf688f0a023b3a043f17905a5ab4a8ec757642a1f72753c735ab3`.

The fixed Tier A+B 265-row manifest, seed 20260801, eight embedding attempts,
stereo-safe strict MMFF94, diagnostic-only ring torsions and a 200-iteration
limit match the [published baseline](2026-10-03-a6-published-v1031-mmff94-quality.md).
The runner was `scripts/public_package_3d_chematic.py` with
`--package-kind source_candidate --tiers AB --arms
chematic_pipeline_v2_mmff94_strict_stereo_safe`; it imported the extracted
wheel, not the checkout. The raw JSONL SHA-256 is
`ef0d0052dd921fa40c1a4ffec0b7e400468215586edbdead8597886ed1a635f0`;
metadata SHA-256 is
`717758b60c66d35cd523ebe93e7a129e41da3323c8f6a7c829d673b46e1464dd`.
The diagnostic `scripts/summarize_a6_mmff94_termination.py` produced
summary SHA-256
`cfa373304cf2f3a82f5b1e088e497b80851c2db5d8bb1ef39887f67bddfc86e8`;
its fixed-265 success-retention gate is **false**. These raw temporary files
are not checked into Git.

| 200-iteration source outcome | Rows |
|---|---:|
| Successful, gradient converged | 101 |
| Successful, iteration cap | 160 |
| Successful, geometry constraint rejected fallback | 1 (row 233, iteration 58) |
| Typed pipeline failure | 3 (rows 240, 246, 248) |

All 262 successful rows report `final_validation.sound=true`; the
`converged` flag agrees with the typed termination for every successful row.
The published v1.0.31 baseline returned 265/265 successes and 100 converged
at 200 iterations. The profiles differ, so this is an **adoption warning**,
not proof of a release-build regression:

| Row | Source failure | Published v1.0.31 |
|---|---|---|
| 240 | Final stereo verification violation | Success, non-converged |
| 246 | Force-field refusal: excessive residual force, iteration 20, constraint-rejected fallback, max component 557.58 kcal/mol/Å | Success, non-converged |
| 248 | Final stereo verification violation | Success, non-converged |

For the four published early-stop rows (170, 191, 227, 234), this source
wheel reports `iteration_limit` at 200 for each. Under an exploratory
400-iteration run, all four report `gradient_converged` at iterations 397,
319, 311 and 357. The underlying source typing and trajectories have changed:
this **does not identify** the reason the published v1.0.31 wheel stopped early.
Nor does it show that raising the default limit is safe. Geometry, stereo,
clash, same-coordinate per-term energy and independent conformer quality were
not scored for this source candidate by the external scorer.

An optimized macOS source wheel was built locally (SHA-256
`a9f0185ab92f51f89d10807f384bc37a7f3557802668a6cb5259561840292883`)
but could not be imported on this host: macOS reported a misaligned Mach-O
LINKEDIT string pool. Rebuilding with `maturin --strip true` did not fix
the loader error. No optimized-wheel 265-row result was obtained; the
independent Linux/macOS CI wheel lanes or another working host must decide
whether rows 240/246/248 also fail under release settings.
The new CI diagnostic records the corresponding optimized source-wheel rows,
metadata and summary on both hosts; until those artifacts are read, the
release-profile outcome is unknown.

The row-246 Python `PipelineV2Error.diagnostics.cause.force_field_bridge_error`
contains `reason=excessive_residual_force` and
`mmff94_termination=constraint_rejected_fallback`; these are separate
failure and stopping conditions. The WASM error JSON exposes the analogous
`mmff94Termination` field.

Next: diagnose rows 240/246/248 against matched starting coordinates and
MMFF atom types; rerun the full independent A6 quality matrix before
considering publication.
Keep the 200/400 published diagnosis separate from these source results.
