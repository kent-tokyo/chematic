# Published v1.0.31 MMFF94 atom-type census

This is an atom-type comparison, **not** a force-field energy, 3D geometry,
or convergence result. The published PyPI `chematic==1.0.31` macOS arm64
wheel used for the [265-row quality rerun](2026-10-03-a6-published-v1031-mmff94-quality.md)
was installed in an isolated Python 3.13 environment. Its pinned wheel
SHA-256 is `b481316bb5a49cb572b035c805316e38a381f4f0e6fa141a05d1da3dc6f05c1f`.
The comparator was RDKit **2026.03.6**. Both engines parsed the same
`rdkit-js-browser-10k-v1.smi` bytes (SHA-256
`f48777e96f74738336f47f682fbff5af6775a3474fe25daa28b595340feee95f`).
Heavy-atom order and explicit hydrogens were aligned by the existing census
script before numeric MMFF94 type comparison.

| Outcome | Count |
|---|---:|
| Input molecules | 10,000 |
| Comparable molecules | 9,774 |
| RDKit MMFF94 unsupported | 204 |
| CheMatic typing error/refusal | 22 |
| Heavy atoms compared / differing | 214,990 / **2,908** |
| Explicit H compared / differing | 190,737 / **56** |

The entire heavy-type, hydrogen-type and row-status breakdown is **identical**
to the published v1.0.25 census on these same bytes. None of the 56 hydrogen
type differences has an agreeing parent heavy-atom type. The largest heavy
bucket is carbon RDKit type 37 versus CheMatic type 2 (1,322 atoms); next are
carbon 64 versus 37 and carbon 63 versus 37 (203 atoms each). These are
prioritization buckets, not yet an adjudication of which implementation is
chemically correct. A separate pass over the same 10,000 inputs found the
1,322 type-37/type-2 atoms in 317 molecules; RDKit marks all 1,322 aromatic.
That makes aromaticity/typing interaction the first investigation target,
not proof that changing aromaticity alone would fix these types. The
[machine census](../validation/results/mmff94-atom-type-census-v1.0.31-2026-10-03.json)
contains every bucket and row example; SHA-256
`aa32fa942bb5f3a7b9d39f363dc2e1c60f64838df5133c2d0cbce77188c6a28c`.

Reproduce with the pinned wheel and RDKit imports in an isolated Python 3.13
environment:

```sh
python scripts/mmff94_atom_type_census.py \
  --corpus validation/benchmark_corpora/rdkit-js-browser-10k-v1.smi \
  --output /tmp/mmff94-atom-type-census-v1031.json
```

The A6 gate remains open: 2,908 current heavy-atom numeric type differences
and 22 CheMatic error/refusal rows require category-specific adjudication.
The 265-row geometry/stereo/clash result and same-coordinate total-energy
result are separate measurements; neither makes these atom types equivalent.
