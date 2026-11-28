# Square-planar planarity and stereo screening: development probe

This probe reruns the 40 fixed cases from the earlier
[bounded-retry development probe](2026-11-27-square-planar-retry-development.md):
16 previously exposed development attempts and 24 authored synthetic attempts.
The input and threshold are not an independent holdout or a representative
population. The checked method uses the same angle policy and bounded retry,
then an explicitly selected maximum unit-vector triple product of 0.3 and the
library's other-declared-stereo verifier. The 0.3 value is an example, not a
validated physical cutoff.

| Group | Attempts | Returned | Embedding failed | Planarity limit failed | Other stereo failed |
| --- | ---: | ---: | ---: | ---: | ---: |
| Exposed development | 16 | 5 | 4 | 6 | 1 |
| Authored synthetic | 24 | 8 | 7 | 9 | 0 |

The checks run in that order; a structure rejected earlier is not assessed by
later screens. The other-stereo failure was the exposed selenium case at seed
42. Its second selenium center was reported as unevaluable because its
coordination is unsupported by the verifier. This does not establish that its
stereochemistry is wrong. All 40 unmodified retry results reproduced the prior
wheel's coordinate/attempt or error output, and all 13 checked results retained
the corresponding retry coordinates exactly.

For all 29 returned retry structures, an independent NumPy determinant
calculation remeasured the maximum absolute unit-vector triple product. The
largest numerical difference was 2.78×10⁻¹⁶. This checks the calculation on
these coordinates; it does not establish a valid threshold, accuracy of the
conformer, or general stereochemistry. The stereo check itself is limited to
the elements reported by the current verifier, and the checked operation does
not validate bond lengths, energy, interfragment placement or subsequent
optimization. The screen performs no new attempt after a planarity or
other-stereo rejection.

The host ran macOS 27.0.1 arm64, Python 3.13.6, NumPy 2.5.3,
`RAYON_NUM_THREADS=1` and `OPENBLAS_NUM_THREADS=1`. The nonpublished source
wheel was built in release mode with LTO disabled and 16 codegen units. Wheel
SHA-256: `65740f3f6d3d596c26001977d55f0078f98135f88c70b70979702ae75cd0ac1d`;
extension SHA-256:
`7482766512e0cd6fb2ed20f2ff640bce55649d8400d9990a827499b15d0770c1`.
The [probe](../scripts/validation/square_planar_checked_development.py)
SHA-256 is `6f7b1f5c1fe93a8672fe69af736afc120594ed4519e473bf026a05df06e584f8`;
the previous input SHA-256 is
`7a75e3ba32b8c0237e106b9b20e8847d6e4ae1f27541e081012b8bc31b0070ae`.
The raw coordinates, comparisons and errors are in
`validation/results/2026-11-27-square-planar-checked/development.json`.
