# Square-planar bounded retry: development probe

This local source-wheel probe compares `rdkit_embed_square_planar` with the
new opt-in `rdkit_embed_square_planar_with_retry`, using the same input,
seed and original final-angle policy. The retry leaves an initial success
unchanged. After an angle rejection only, it makes one additional whole
embedding call with the same seed, trans bounds tightened by 20° and cis
bounds tightened by 10°, subject to the documented caps. Both outputs must
satisfy the original policy: trans ≥150° and cis 60°–120°. The reported
successes therefore satisfy this screen by construction.

The first group is the eight previously exposed development molecules and two
seeds from the [earlier probe](2026-11-27-square-planar-embedding-development.md).
The second group contains 12 authored synthetic Pt, Pd and Ni structures with
seeds 7 and 23. The second group was authored after exploratory tuning of the
retry rule. Neither group is an independent holdout or a random chemical sample.

| Group | Calls | Single attempt returned | Retry API returned | Returned on second call | Both calls failed |
| --- | ---: | ---: | ---: | ---: | ---: |
| Exposed development | 16 | 5 | 12 | 7 | 4 |
| Authored synthetic | 24 | 1 | 17 | 16 | 7 |

Every one of the six initial successes retained exactly the same coordinates
and reported one embedding call. All first failures in this sample were angle
rejections and led to a second call; no setup or initial embedding errors were
observed. The two exposed iodine seeds went from zero to one accepted result;
the selenium case retained its one initial success and one failure. This is an
angle-policy check, not independent evidence of correct stereochemistry, bond
geometry, energy or crystal-structure accuracy. A later force-field
optimization may invalidate the checked angles. The retry can require two
whole embedding calls, and the probe does not establish a speed improvement or
a general success-rate improvement.

An additional local timing probe ran the same 40 calls in single, retry,
retry, single block order, including failed calls and excluding molecule
parsing and hydrogen addition. The four block medians were 0.331, 0.580,
0.563 and 0.295 ms per call respectively. In each retry block, 34 calls made a
second embedding call and six did not. These are small molecules on one
non-isolated host while Rust tests were running, with no warmup or general
runtime guarantee. The
[timing script](../scripts/validation/square_planar_retry_timing.py) SHA-256
is `8f3a5b0baa5f20caec5d7535b1764bb14a95bb090c91d1e2707fbec9b9243d35`;
the per-call timings are retained in `timing.json`.

As a separate measurement check, RDKit 2026.09.2 re-evaluated the 29 returned
retry structures at the exact stored coordinates after atom-order and bond
connectivity checks. All 126 measured pair angles agreed within 1.14×10⁻¹³°.
This checks the angle calculation, not the correctness of the selected trans
pairs or whether either engine would generate the same geometry. The
[verification script](../scripts/validation/square_planar_retry_rdkit_angles.py)
and `rdkit-angles.json` retain the individual comparisons; its script
SHA-256 is `c84dc6a04ef43b7aa5aff6b359f6c7275e04f013f6a7b3c6afb424c2fa8f39b8`
and the input JSON SHA-256 is
`7a75e3ba32b8c0237e106b9b20e8847d6e4ae1f27541e081012b8bc31b0070ae`.

The host ran macOS 27.0.1 arm64 and Python 3.13.6 with `RAYON_NUM_THREADS=1`.
The nonpublished source wheel was built in release mode with LTO disabled and
16 codegen units. Wheel SHA-256:
`2e46dc1e0c1ce15be9a71243e4eb29d91428e0b718643eb779346c0e00897f41`;
extension SHA-256:
`13899a39006757d03244bca678a210bc16710bf87ff20f19cab1dc280c62b24d`.
The [probe](../scripts/validation/square_planar_retry_development.py)
SHA-256 is `eb91441d9fd7b4c7c0b4ae3aceff757a11ad9d68e2939010652a5e35b5f2d619`.
Raw input SMILES, seeds, coordinates, measured angles and error messages are
in `validation/results/2026-11-27-square-planar-retry/development.json`.
Run the probe from the repository root with that wheel installed and specify
`--wheel PATH --output PATH` to reproduce the format.
