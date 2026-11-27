# Reasons for timed-conformer geometry-screen rejections

The optional completed-conformer screen now reports a disjoint per-conformer
breakdown: `stereo_violated`, `distorted_tetrahedral_geometry`, and
`stereo_unevaluable`, in that precedence order. The new breakdown describes
rows returned by the backend and then removed by the optional screen. It is
separate from backend rejection counters; it is not an atom-count histogram.
Malformed or inconsistent diagnostic results fail the worker rather than being
assigned a guessed category. Rejected coordinates are not exposed.

## Public timed-API comparison

`after.json` repeats the 12 exposed development molecules and both normal and
partial timed modes from `../2026-11-28-timed-conformer-regular/`, with the same
seed, MMFF iteration limit and 20-second process deadline. The previous output
was produced with local wheel SHA-256
`7f0c988b4d1a537542e88df914291f742b005c252a719c1088b8e790c1b3625b`; the
candidate local wheel SHA-256 is
`46715cfd22e52931ed1f105f46a1d218bc5e016ccb34aa3ce36ca52dfc093fb0`.
The candidate was built from commit `92045c21` plus the accompanying source
changes using macOS arm64, CPython 3.13.6 and `maturin build --locked`. Neither
wheel is a published distribution. The exact selection and protocol are in the
prior directory. Reproduce using the candidate wheel:

```sh
python validation/results/2026-11-28-timed-conformer-regular/run.py after > after.json
```

All 24 pairs completed without exceptions, with no geometry-screen rejections.
The prior and candidate results match after removing elapsed time and the new
reason metadata. Every new reason map sums to its rejection count. This small
sample does not establish rejection rates or a speed, accuracy or parity claim.

## Simulated-candidate diagnostic

`diagnose.py` wraps one rule-based `generate_3d_etkdg()` geometry from each of
the preceding 30 exposed development inputs as a **simulated** candidate row.
It does not run the timed backend or imply that these rows would pass that
backend. Of these 30 rows, 3 were retained, 17 classified as stereo violations
and 10 as distorted tetrahedral geometry; none raised a processing error. The
report also includes separate controls for each of the three reason labels,
including an unsupported square-planar case. The controls and these
unoptimized geometries do not establish chemical or conformer quality.

```sh
python validation/results/2026-11-28-timed-rejection-reasons/diagnose.py > diagnose.json
```

`crates/chematic-py/tests/test_timed_conformers.py` additionally checks the
three labels, precedence, force-field fallback, counter reconciliation,
interrupted-worker recovery and legacy checkpoints. A legacy checkpoint with
an unrecorded nonzero breakdown yields `None`, not guessed zeroes; a recorded
zero rejection count has an exact all-zero breakdown.
