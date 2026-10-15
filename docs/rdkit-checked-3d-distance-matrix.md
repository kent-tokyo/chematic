# Checked 3D distance matrices

The existing RDKit-model distance matrix already applies the allocation guard
motivated by [RDKit #9647](https://github.com/rdkit/rdkit/issues/9647).
The new `rdkit_distance_matrix_3d_checked` adds an opt-in finite-result contract.
The numeric cases below are local follow-up findings, not additional upstream
heap-corruption claims. The issue's oversized crash reproducer was not run.

## Contract

- Original and prepared inputs are limited to 1,024 atoms / 1,048,576 entries.
  The native size check precedes model preparation and finite-value validation.
  Flat and row buffers use the existing fallible reservation helpers.
- Coordinates must match prepared RDKit atom order, including explicit-H
  removal and added-H state, and every component must be finite.
- `use_atom_wts=True` sets diagonals to `6 / atomic number`. Atomic number
  zero has no defined weighted diagonal and returns an error. With weights
  disabled, dummy atoms remain usable and their diagonals are zero.
- Ordinary finite squared norms use the legacy expression unchanged. When
  squares overflow, are subnormal, or underflow to zero for distinct points,
  the kernel uses a scaled `hypot` norm. A representable distance of `1e308`
  or `1e-300` therefore remains finite and nonzero.
- A genuinely unrepresentable norm, for example between `1e308` and `-1e308`
  on one axis, returns an error. No partial matrix is exposed.
- An empty graph with empty coordinates returns an empty matrix; a nonempty
  coordinate array for an empty graph rejects.

This inherits the legacy RDKit-model preparation contract. It is not an
unprepared native-graph distance routine. The cap bounds native f64 storage
at approximately 16 MiB while flat and nested buffers coexist, excluding
binding objects and JSON. Coincident distinct atoms have a valid zero distance;
the API does not reject a conformer merely for coincident coordinates.

## Use

```python
import chematic
mol = chematic.from_smiles("CC")
coords = [[0, 0, 0], [1e308, 0, 0]]
matrix = mol.rdkit_distance_matrix_3d_checked(coords)
assert matrix[0][1] == 1e308
```

Rust exports `chematic_smiles::rdkit_distance_matrix_3d_checked`, returning
`Result<Vec<Vec<f64>>, RdkitSmilesError>`. Python exposes the corresponding
`Mol` method and maps errors to `ValueError`. Coordinate rows must have exactly
three components.

WASM exposes `MolHandle.rdkit_distance_matrix_3d_checked_json(coords_json,
use_atom_wts)`. It also applies the shared 1,000,000-byte JSON-input and 16 MiB
JSON-output bounds. Nonfinite JSON numbers and malformed rows refuse during
JSON decoding. An otherwise valid native matrix can exceed the WASM JSON-output
budget, in which case that binding returns an error.

The original `rdkit_distance_matrix_3d` remains available with its pinned
compatibility behavior, including extreme/nonfinite numeric outcomes. Normal
coordinates do not require selecting a different preparation or atom ordering.

## Evidence

The original recipes in `validation/rdkit_issues/checked-3d-matrix/` include
finite, extreme, undefined-diagonal and nonfinite-coordinate cases. Nonfinite
values use string spellings in the recipe JSON; the Python runner converts
them to actual floats. Those spellings are not valid numeric input for WASM.

The reproducible runner compares 128 ordinary configurations against RDKit
2026.09.1 and the legacy chematic API. The observed maximum RDKit difference
is **0.0**, with **zero** legacy disagreements on this packet. Eight numeric
boundary recipes additionally verify checked acceptance/refusal and retain
legacy/RDKit outcomes as categories so the evidence JSON contains no NaN or
Infinity literals.

```sh
cargo test -p chematic-smiles --lib --test rdkit_matrix_3d_checked --test rdkit_matrix_limits --offline --quiet
python3 -m pytest -q crates/chematic-py/tests/test_rdkit_matrix_3d_checked.py
node crates/chematic-wasm/tests/checked_3d_matrix.test.mjs
python3 scripts/validate_checked_3d_matrix.py
```

Local validation on main `2c3ff6ff` plus this change passed **406 Rust tests**
(9 new; one existing ignored), **26 Python tests** (15 new), and the actual
Node/WASM numeric/input/allocation contracts. WASM also matched all 128
ordinary packet matrices exactly across runtimes. Clippy with warnings denied,
formatting, binding inventories, refactor budgets and strict documentation
builds passed.

Python needs a source-built extension and the pinned RDKit runtime. Node uses
`pkg-node` by default or an explicit `CHEMATIC_WASM_MODULE`. Source/runtime
hashes, test counts and the full comparison packet are stored under
`validation/results/rdkit-checked-3d-matrix-2026-10-16.*`. These are source
candidate checks, not evidence of package publication or general RDKit parity.
