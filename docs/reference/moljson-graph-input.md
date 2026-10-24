# MolJSON graph input integrity

The existing chematic atom-record MolJSON reader now retains quadruple bonds
and checks known graph fields before constructing a molecule. Previously,
`charge: 128` became `-128`, `isotope: 65537` became `1`, and a non-aromatic
`order: 4` became aromatic `1.5`. Duplicate IDs could redirect bonds and failed
bond insertions could be ignored. These interpretations are now corrected or
explicitly refused.

```python
import json
import chematic

source = json.dumps({
    "atoms": [{"id": "left", "element": "Cr", "hydrogens": 0},
              {"id": "right", "element": "Cr", "hydrogens": 0}],
    "bonds": [{"source_id": "left", "target_id": "right",
               "order": 4, "aromatic": False}],
})
result = json.loads(chematic.from_moljson(source).to_moljson())
assert result["bonds"][0]["order"] == 4
assert result["bonds"][0]["aromatic"] is False
```

## Known-field contract

- A present atom `charge` must be an integer fitting signed 8-bit storage
  (`-128..=127`). Omission retains the zero default.
- A present `isotope` is either null or an integer in `0..=65535`. Omission
  retains absence. Zero is retained as supplied; storage bounds do not certify
  an isotope's physical existence or a charge's chemical plausibility.
- Atom and bond `aromatic` fields must be booleans when present. Null, strings
  and numeric flags do not become false defaults.
- Present numeric bond `order` is exactly `1`, `1.5`, `2`, `3` or `4` after
  JSON number decoding. These mean single, aromatic, double, triple and
  quadruple. Approximate values are not rounded or matched with a tolerance.
  Explicit aromatic flags must agree with the order. When order is absent,
  `aromatic: true` means aromatic; otherwise the existing single-bond default
  remains. When the flag is absent, `1.5` determines aromatic bond order.
- Atom IDs must be nonempty, unique strings. Present bond IDs must also be
  nonempty, unique strings in their separate namespace; omission is still
  accepted. Endpoint references must resolve to existing atoms. Self bonds,
  repeated endpoint pairs and reversed duplicates are refused. Failed bond
  insertion propagates as an error.
- A present `bonds` field must be an array. An omitted field remains an empty
  graph edge list. Duplicate object keys are refused recursively, including
  escaped spellings that decode to the same key. The existing shared unique-key
  JSON decoder is reused; MolJSON's own resource and graph checks still apply.

Known present fields with the wrong type are errors rather than neutral
defaults. Integer fields do not coerce floating/exponent forms or booleans.
JSON number decoding still uses binary64 for fractional values; this is not an
arbitrary-precision decimal contract. The existing authoritative
[non-node hydrogen counts](moljson-hydrogens.md) remain unchanged.

Rust `parse_moljson` and its bounded variant, Python `from_moljson`, actual
WASM `mol_from_moljson`, their common conversion routes and CLI `convert` use
this reader. Missing required fields retain the existing `MissingField` error
contract; graph validation errors use existing error variants. No new binding
export, dependency or format alias is introduced.

## Representation boundary

This is chematic's existing **atom-record dialect**, not the
[oxpig MolJSON schema](https://github.com/oxpig/MolJSON#format-notes), whose
charge and aromatic-H lists have a different layout. IDs resolve source
references and are reassigned by the writer; they are not retained as source
labels. Source atom order and indexed bond endpoints remain observable.

The validated subset is element, charge, isotope, non-node H count, atom/bond
aromatic flags, connectivity and the five supported numeric bond orders.
Unknown fields are still outside the graph projection. Source bytes, comments,
coordinates, atom maps, stereo and spin multiplicity are not preserved by this
contract. This is not a valence/aromaticity sanitizer or a general lossless
MolJSON writer. The legacy native writer still projects zero/dative/stereo/query
bond orders to single bonds; this change does not broaden that writer domain.
Use a format/API that represents those fields when they must survive conversion.

## Local evidence, 2026-10-25

The original packet `validation/moljson_graph_input/cases.json` contains 13
positive graphs and 51 malformed or ambiguous graph inputs. Authored graph
oracles include all declared atom fields, source-order indices, bond endpoints,
exact orders and flags. Cases include single/double/triple/quadruple/aromatic
controls, isotope/charge bounds, omitted-field defaults, reversed references,
Unicode IDs and an empty graph. Bounds cases test representation, not chemical
validity. Native and binding tests perform three successive read/write cycles.

| Route | Earlier local graph preservation | Candidate graph preservation | Earlier / candidate explicit refusals |
|---|---:|---:|---:|
| Python direct / common, each | 11/13 | 13/13 | 2/51 → 51/51 |
| Actual Node/WASM direct / common, each | 11/13 | 13/13 | 2/51 → 51/51 |
| CLI conversion | 11/13 | 13/13 | 2/51 → 51/51 |

The two changed positive graphs contain quadruple bonds. Refusals are counted
separately from preserved graphs. Observer controls distinguish validation errors from panics/timeouts and detect changes in bond
order/aromaticity, charge/isotope/H values and types, missing edges, endpoint
changes and duplicate emitted JSON keys. Native tests also inspect the molecule
directly, independently of writer normalization.

The baseline uses preserved local runtimes with the native MolJSON reader source
at `d6de9d6e`, verified byte-identical to comparison parent main `7704ee20`.
It is not a release build. Receipts retain fixture/source/runtime hashes,
complete outputs/errors and bounded CLI/Node process records. Processes have a
10-second and 1-MiB combined-output bound.

Open Babel 3.2.1's captured installed format inventory has no MolJSON entry.
It receives **no same-format score**. These results close a local graph-integrity
gap; they do not establish Open Babel-wide superiority, upstream oxpig schema
support, speed, an independent holdout or a published-package result.

Receipts are `validation/moljson_graph_input/baseline-2026-10-25.json` and
`validation/moljson_graph_input/candidate-2026-10-25.json`. Reproduce with a fresh
Python extension on `PYTHONPATH`, the CLI binary and Node-target WASM:

```sh
python scripts/validation/check_moljson_graph_input.py \
  --cli /path/to/chematic --wasm-module /path/to/chematic_wasm.js \
  --output /tmp/moljson-graph-input.json
python -m unittest scripts.tests.test_moljson_graph_input
```
