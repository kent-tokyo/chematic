# Stereo repair on the exposed boron development subset

This is the 24-molecule exposed COD development selection from
`../2026-11-27-boron-initial-geometry/before.json`, applied to the same source
input. The selection contains 14 molecules for which the verifier reports no
stereo; no reserved holdout was used. `run.py` records every result or exception,
per-element statuses and failures, and the number of nonbonded atom pairs under
0.8 Å. It uses the molecule's own bond table to compare bond lengths before and
after successful repair. This is a diagnostic screen, not an independent
chemical geometry or steric-quality oracle.

A local debug wheel was built from commit `3b090dc0` plus the accompanying source
changes, with CPython 3.13.6 on macOS arm64 and `maturin build --locked`.
Its SHA-256 is `4741ca14f243cf970e8637edf0b0431f01ba8f9070715eca859cf82ad79c103b`.
It is not a published wheel. To repeat against a candidate installed wheel:

```sh
python validation/results/2026-11-27-boron-stereo-repair/run.py > results.json
```

All 24 inputs ran without an exception. Twenty returned without repairs (14 of
these had no reported stereo), two were repaired, and two returned failure
because required ring substituents were not bridge-isolatable. The successful
repairs addressed one center in record 2877 and two in record 9498. For 2877,
one atom moved up to 3.472 Å; for 9498, the largest movement was 1.518 Å. No
successful repair changed any bond length by more than 1e-12 Å or changed the
count of nonbonded pairs closer than 0.8 Å. Record 9498 retained four such pairs,
but pairwise reinspection found that one original pair was replaced by a new
close pair (atoms 3 and 32). Unchanged counts must not be read as unchanged
contacts or clash-free geometry. All successful
outputs have zero violations and zero unevaluable elements according to the
same verifier used by the repair implementation.

A separate Python regression test asks RDKit to reassign stereochemistry from
both the original and repaired record 2877 coordinates: the original center
has the opposite CIP code from the declared molecule, while the repaired center
matches it. That independent check is limited to this one center. Neither this
sample nor this test establishes general RDKit parity, conformer accuracy, or
quality after force-field minimization.
