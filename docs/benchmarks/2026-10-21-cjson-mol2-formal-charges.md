# CJSON → MOL2 formal-charge preservation, 2026-10-21

The checked converter previously lost nonzero formal charges except +1 nitrogen
encoded by `N.4`. It now emits explicit `UNITY_ATOM_ATTR` charge records and
requires exact charge equality through verified atom correspondence. Partial
charges remain independent floating-point data, with their own method label.

UNITY precedes BOND: Open Babel 3.2.1 reads this order and skips trailing UNITY
in the previously measured [section-order diagnostic](../mol2-section-order.md).
Its [MOL2 documentation](https://openbabel.org/docs/FileFormats/Sybyl_Mol2_format.html)
also identifies UNITY records as formal-charge information. Element-only atom
types remain the adapter's limited typing profile. Open Babel warns about
oxygen types in the oxide and nitrate readbacks; those warnings are archived.
No general atom typing, hydrogen assignment or stereochemistry claim follows.

## Same-input development measurement

Six original synthetic inputs cover a helium neutral control, NaCl, CaCl2,
oxide, zinc and a connected nitrate ion. These are supplied coordinates for
file conversion, not optimized geometries. Both commands receive identical
canonical JSON bytes. An independent small Python reader scores the indexed
element/bond graph, formal charges, coordinate binary64 bits and any declared
partial-charge bits/method. Source arrays define expected values; neither
toolkit is used as the expected-value oracle. No tolerance is allowed.

| Required field | Previous adapter | Fixed adapter | Open Babel 3.2.1 default |
|---|---:|---:|---:|
| Indexed elements/bonds and formal charges | 1/6 | 6/6 | 1/6 |
| Coordinate bits | 6/6 | 6/6 | 5/6 |
| Declared partial-charge bits | 3/3 | 3/3 | 0/3 |
| Declared charge method | 3/3 | 3/3 | 0/3 |
| All fields required by this packet | 1/6 | 6/6 | 1/6 |

All twelve paired commands completed with output; no failure or refusal is
excluded from the denominator. An additional six commands read the fixed
adapter's MOL2 using Open Babel and write SDF. Independent inspection of SDF
elements, explicit bonds and formal charges agrees with all six emitted
graphs. That readback verifies interoperability of emitted charges, not SDF
coordinate precision or partial-charge preservation.

The before-fix production source is `7d773b72` (1.2.0 development main).
The baseline was captured before the production edit, using the same example
adapter and six inputs. The final harness additionally archives external
readbacks. Raw inputs, hashes, outputs, warnings, terminal outcomes and
per-field scores are retained in the
[baseline](https://github.com/kent-tokyo/chematic/blob/main/validation/results/cjson-mol2-charges-baseline-2026-10-21.json),
[candidate](https://github.com/kent-tokyo/chematic/blob/main/validation/results/cjson-mol2-charges-candidate-2026-10-21.json)
and [validation index](https://github.com/kent-tokyo/chematic/blob/main/validation/results/cjson-mol2-charges-validation-2026-10-21.json).

These inputs were authored during development and include the motivating
failure. This is a bounded improvement on the declared CJSON profile, not an
independent holdout, general CJSON-schema conformance, released-package proof
or toolkit-wide superiority. Names, arbitrary metadata, isotopes, implicit
hydrogens, generated charge columns and stereochemistry are outside this
packet's scoring. The API continues reporting source-text loss, inferred
MOL2 metadata and unverified stereo; strict policy still refuses.

## Regression checks and reproduction

Rust regressions cover every i8 charge value for three elements, distinct
formal/partial charges, sparse charged atom IDs, absence of partial-charge
columns, exact output-byte limits and withheld correspondence after isotope
loss. Extreme charge values test the representation range, not chemical
validity. Existing finite-coordinate and partial-charge bit tests remain.
The packet oracle has negative controls for missing charges, altered bonds,
signed zero, one-ULP differences and malformed UNITY/SDF references.

```sh
cargo test --locked -p chematic-mol --test cjson_mol2_conversion
cargo build --locked -p chematic-mol --example cjson_mol2_report
python3 -m scripts.check_cjson_mol2_formal_charges \
  --probe target/debug/examples/cjson_mol2_report \
  --obabel /opt/homebrew/bin/obabel --output /tmp/cjson-mol2-charge-new.json
python3 -m pytest -q scripts/tests/test_cjson_mol2_formal_charges.py \
  scripts/tests/test_interchange_preservation.py
```

The comparator is pinned to Open Babel 3.2.1; its binary identity and version
are in the result. Output paths must be new. Fixtures are original project
contributions under `MIT OR Apache-2.0`; the checked-in recipe and result
record provenance, minimization, expected fields and normalized input hashes.
