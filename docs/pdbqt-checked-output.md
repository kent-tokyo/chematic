# Checked rigid PDBQT output

[Open Babel #2101](https://github.com/openbabel/openbabel/issues/2101) reports
coordinates overflowing adjacent PDBQT columns. Rust formatting widths are
minimum widths too: the legacy `write_pdbqt` can reproduce the same defect.
The new `write_pdbqt_checked` refuses an unrepresentable block before returning
any output. Existing writers and parser error enums retain their behavior.

```python
import chematic

mol = chematic.from_smiles("CCO")
pdbqt = chematic.write_pdbqt_checked(
    mol,
    [(0.0, 0.0, 0.0), (1.54, 0.0, 0.0), (2.5, 1.2, 0.0)],
    [-0.1, 0.05, -0.4],
    name="ETH",
)
```

Coordinates are in Å and partial charges in elementary charge units. Supply
exactly one tuple and charge per **graph atom**, including explicit hydrogen
atoms. Empty or excess arrays fail; no zero defaults, implicit-H materialization,
atom removal, hydrogen merging or charge recalculation are applied.

Rust exports `write_pdbqt_checked(&Molecule, &[(f64, f64, f64)], &[f64], &str)`
returning `Result<String, PdbqtWriteError>`. The new error enum is non-exhaustive;
atom indices are zero based. Python exposes the free function above and raises
`ValueError`. JavaScript/WASM exposes
`smiles_to_pdbqt_checked(smiles, coords_json, charges_json, name)` returning a
string or throwing on failure. Its coordinate JSON contains `[x,y,z]` arrays.
Malformed Python argument shapes/types fail at the binding boundary.

## Output contract

The layout follows the
[AutoDock 4.2.6 user guide, PDBQT format specification](https://autodock.scripps.edu/wp-content/uploads/sites/56/2022/04/AutoDock4.2.6_UserGuide.pdf).

| Field | One-based columns | Validation |
| --- | --- | --- |
| Atom serial | 7–11 | Sequential 1–99999; empty molecules and larger counts fail |
| x, y, z | 31–38, 39–46, 47–54 | Finite; actual fixed-decimal text must fit 8 columns |
| Partial charge | 71–76 | Finite; actual fixed-decimal text must fit 6 columns |
| AutoDock type | 78–79 | Existing heuristic; wildcard atoms and unknown `X` mappings fail |
| REMARK name | Separate header | Empty or up to 1024 printable ASCII bytes including spaces; controls and non-ASCII fail |

Coordinates and charges are rounded to **three decimals**. Checks run after
rounding and include the sign: `9999.9999` becomes `10000.000` and is refused;
`-999.9999` becomes `-1000.000` and is refused. Charge values `99.9999` and
`-9.9999` likewise overflow. Representable rounded coordinates range from
`-999.999` to `9999.999`, and charges from `-9.999` to `99.999`.
Negative zero remains representable. No clipping, recentering, column widening
or silent partial output is performed. Callers retain original numeric values
if they need more precision than this format provides.

Accepted input produces exactly the legacy writer's bytes, with all explicit
atoms in source order, synthetic residue identity `LIG`/`A`/`1`, occupancy 1 and
B-factor 0. Output is rigid and has no ligand torsion tree. AutoDock types use
the existing heuristic; this change does **not** establish chemical correctness
of acceptor/donor typing, protonation or docking preparation. PDBQT output is
not a lossless molecular graph/isotope/stereochemistry round trip.

WASM additionally applies existing input budgets: 1 MB SMILES, 100 KB per JSON
string and 10000 parsed atoms. The 1024-byte name cap also applies. JSON rejects
NaN, Infinity and out-of-range numeric literals; the Rust kernel separately
checks all native and Python floating-point values.

## Reproduction and validation

The original fixture `validation/openbabel_issues/pdbqt-checked/coordinate-overflow.xyz`
contains one N atom at `(2261608.729, 1334600.251, 0)` Å. It deliberately tests
serialization limits, without gen3D randomness or the issue's attachment.
Open Babel 3.2.1 on this macOS arm64 host exits successfully but emits the
concatenated coordinate text `2261608.7291334600.251`. The legacy native writer
also emits this pattern. Checked APIs reject it. This reproduces the field-width
failure, not the historical Windows coordinate-generation cause.

The evidence runner records 27 accepted cases and 33 refusal cases. Accepted
cases include explicit H, aromatic C, two-letter elements, all supported default
elements, translated coordinates and rounded width boundaries. Python and real
Node/WASM outputs agree byte for byte. A separate fixed-column decoder checks
coordinates/charges/serials, and a small original C++ probe linked to Open Babel
reads every accepted block with automatic bonding disabled, checking all indexed
elements, coordinates and stored partial charges. This tests field preservation,
not molecular connectivity or AutoDock type quality.

Eight Rust contract tests also exercise all numeric components, NaN/Infinity,
exact counts, header injection, wildcard/unknown types, typed errors and the
99999/100000 atom serial boundary. Python regression tests and the real
Node/WASM contract include malformed binding input and WASM limits. The existing
public API compatibility test checks that older exhaustive error matches still
compile. CI includes the new Rust integration target.

Build the diagnostic against the installed Open Babel, then run with a Python
extension and Node module built from this source:

```sh
clang++ -std=c++17 \
  validation/openbabel_issues/pdbqt-checked/openbabel_read_probe.cpp \
  $(pkg-config --cflags --libs openbabel-3) -o /tmp/pdbqt-read-probe
python scripts/validation/validate_pdbqt_checked.py \
  --wasm-module /path/to/node/chematic_wasm.js \
  --openbabel-probe /tmp/pdbqt-read-probe
cargo test -p chematic-mol --test pdbqt_checked_contract \
  --test public_api_compat_v1_1_3 --offline
python -m pytest crates/chematic-py/tests/test_pdbqt_checked.py -q
node crates/chematic-wasm/tests/pdbqt_checked.test.mjs /path/to/node/chematic_wasm.js
```

Raw inputs, outputs, refusal diagnostics, the independent cross-read and runtime/
source hashes are retained in `validation/results/pdbqt-checked-2026-10-16`.
Local checks and hosted CI are separate. No publication, docking-quality,
performance or general Open Babel superiority claim follows from these cases.
