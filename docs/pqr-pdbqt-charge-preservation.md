# PQR-to-PDBQT charge preservation

This local source candidate addresses the charge-preservation workflow in
[Open Babel #2845](https://github.com/openbabel/openbabel/issues/2845).
The existing PQR parser already keeps charge and radius in `PqrAtomRecord`;
its plain molecule projection has no storage for those source fields.
The new adapter reads directly from those records.

## Use

```python
import chematic

source = """HETATM 1 ZN ZN A 1 0 0 0 2.0000 1.0900
HETATM 2 CL CL A 2 20 0 0 -1.0000 1.7500
"""
report = chematic.pqr_to_pdbqt(
    source, ["Zn", "Cl"], allow_radius_loss=True
)
pdbqt = report["pdbqt"]
assert report["source_charges"] == [2.0, -1.0]
assert report["source_radii"] == [1.09, 1.75]
assert report["radii_omitted"]
```

Rust exports `pqr_to_pdbqt(&PqrResult, &[&str], bool)` returning
`Result<PqrPdbqtReport, PqrPdbqtError>` from `chematic_mol`.
Python accepts the PQR text and returns a dictionary containing the same
report fields, plus `charge_provenance="pqr_atom_records"` and
`atom_type_provenance="caller_supplied"`. Invalid input raises `ValueError`.

## Contract

- Preserve source atom order, serials, atom and residue names, residue numbers,
  optional chain and insertion code, record kind, and every explicit H.
- Copy partial charges directly; do not infer bonds, recompute charges, change
  protonation or merge hydrogen charges. Values are in elementary charge units.
- Emit fixed-width PDBQT with three decimal places for charges and coordinates.
  Preserve original charges, radii and coordinates in source-order sidecar
  arrays. Report maximum absolute charge and coordinate rounding errors;
  coordinate and radius units are Å.
- Require explicit `allow_radius_loss=True`; PDBQT has no radius field.
  A default call fails rather than silently discarding radii. The sidecar
  preserves them; the PDBQT payload alone does not.
- Require caller-supplied AutoDock types in source atom order. Validate the
  supported type and element match. The caller supplies aromaticity and donor/
  acceptor chemistry; this adapter cannot derive them from a bondless PQR file.
- Reject non-finite numbers, negative radii, duplicate or out-of-range serials,
  field overflow after rounding, non-ASCII/whitespace metadata, and names too
  long for their fields. A five-character residue name is rejected even after
  radius-loss consent. Never silently truncate or renumber identity.
- Bound conversion to 99,999 source atoms; the PQR text parser also applies its
  ordinary input, line and atom limits.

The supported types are `C`, `A`, `N`, `NA`, `O`, `OA`, `OS`, `S`, `SA`, `H`,
`HD`, `P`, `F`, `Cl`, `Br`, `I`, `Mg`, `Mn`, `Zn`, `Ca` and `Fe`.
PQR element inference remains the existing documented atom-name heuristic.
This is rigid record adaptation. It supplies no flexible-ligand torsion tree,
force-field validation, pH preparation, nonpolar-H suppression or docking
accuracy guarantee. Keep the report together with the payload when radii or
full precision matter downstream.

The field contracts follow the
[APBS PQR documentation](https://apbs.readthedocs.io/en/latest/formats/pqr.html)
and [AutoDock 4.2.6 manual, Appendix I](https://autodock.scripps.edu/wp-content/uploads/sites/56/2021/10/AutoDock4.2.6_UserGuide.pdf).
No upstream source code or issue attachment was copied.

## Reproduction and evidence

The independent two-ion fixture is
`validation/openbabel_issues/pqr/charged-ions.pqr`. Open Babel 3.2.1 on this
host reads both atoms with the expected elements and coordinates, then writes
both charges as zero. This happens both with its default options and with
`-xr -xh -xp -xn` (rigid output, preserve H, indices and names).
The new adapter writes charges `2.000` and `-1.000`.
The separate five-record fixture exercises mixed record kinds, missing chain,
non-contiguous serials, an insertion code, explicit H and precision rounding.
It is not used as the Open Babel charge-loss comparison because its loose
PQR field spacing also triggers that reader's parsing errors.

```sh
cargo test -p chematic-mol --lib --test pqr_pdbqt_contract --test mmcif_pqr_integration --offline
python -m pytest crates/chematic-py/tests/test_pqr_pdbqt.py crates/chematic-py/tests/test_module_functions.py
python scripts/validate_pqr_pdbqt.py
```

The script requires Open Babel 3.2.1 and an extension built from this source.
It checks source charges against an independent fixed-column decoder, captures
upstream results, and records fixture provenance, license, minimization,
expected values, tolerances, source hashes and extension hash in
`validation/results/openbabel-pqr-pdbqt-2026-10-14.json`.
The fixtures are original project data under `MIT OR Apache-2.0`.

The initial broad `cargo test -p chematic-mol` attempt reached four existing
fuzz-corpus regression failures because this isolated sparse source snapshot
omits `fuzz/corpus`. These are missing-data errors, not passed tests. The
bounded validation command above covers the complete library suite, new
adapter contract and existing mmCIF/PQR integration without that absent corpus.
The focused Rust run passed 767 tests (753 library, 10 new adapter, 4 existing
mmCIF/PQR integration); Python passed 85 tests (17 new, 68 existing format
functions). Clippy with warnings denied, formatting, refactor budget, static
binding inventory, naming and strict MkDocs passed.

This record concerns local source, not a published package or upstream fix.
