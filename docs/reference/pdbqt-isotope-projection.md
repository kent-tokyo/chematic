# PDBQT isotope projection

[Open Babel #3005](https://github.com/openbabel/openbabel/issues/3005) reports
different retained hydrogens, hydrogen types and torsion counts for ordinary
and deuterated acetophenone in Open Babel 3.1.0. The report asks whether this is
intended. This feature addresses the associated isotope-information boundary;
it does not implement flexible-ligand preparation or claim to fix that torsion
count in Open Babel.

PDBQT has no isotope field. Previously, chematic's rigid checked writer accepted
`[2H]O[13CH2][3H]` and emitted ordinary H/O/C/H atoms without notifying the caller
that all three labels were omitted. `write_pdbqt_checked` now refuses any isotope
annotation with a typed `PdbqtWriteError::IsotopeLoss` (Python `ValueError`, WASM
exception), identifying the zero-based atom index and mass number. This includes
heavy-atom labels and explicitly constructed `Some(0)` metadata.

When an isotope-free docking projection is deliberate, use the separate API:

```python
import chematic

mol = chematic.from_smiles("[2H]O[13CH2][3H]")
report = chematic.project_pdbqt_isotopes(
    mol,
    [[0, 0, 0], [1, 0, 0], [2, 0, 0], [3, 0, 0]],
    [-0.2, 0.3, -0.4, 0.1],
)
print(report["pdbqt"])
print(report["isotopes_omitted"])
```

The sidecar is:

```json
[
  {"source_atom_index": 0, "pdbqt_serial": 1, "atomic_number": 1, "mass_number": 2},
  {"source_atom_index": 2, "pdbqt_serial": 3, "atomic_number": 6, "mass_number": 13},
  {"source_atom_index": 3, "pdbqt_serial": 4, "atomic_number": 1, "mass_number": 3}
]
```

Keep this record alongside the exact returned PDBQT block. It records labels
from the original graph, with zero-based source indices and one-based output
serials. Unlabelled input returns an empty list. The PDBQT file itself does not
contain these labels; docking programs do not consume this sidecar. It is not a
verified restoration format for modified or reordered output files.

Rust exports `project_pdbqt_isotopes(&Molecule, &[(f64, f64, f64)], &[f64], &str)`
returning `Result<PdbqtIsotopeProjection, PdbqtWriteError>`. The result contains
`pdbqt: String` and `isotopes_omitted: Vec<PdbqtIsotopeLoss>` and is serializable.
WASM exports `smiles_to_pdbqt_isotope_projection(smiles, coords_json,
charges_json, name)`, returning a JSON string with the same two fields.

## Output contract

- Every explicit graph atom, including ordinary H, D and T, is retained in
  source order. Implicit hydrogens are not expanded.
- AutoDock types use the existing graph-based heuristic: isotope labels do not
  affect them. In the example above, O-bound D is `HD` and C-bound T is `H`.
- Source coordinates in angstroms and caller-supplied partial charges in e are
  rounded to three decimals. No charges are recalculated, merged or folded.
- All [checked rigid output](pdbqt-checked-output.md) limits still apply:
  exact array lengths, finite values, fixed widths, printable name, supported
  default types and at most 99999 atoms (10000 via WASM).
- The source molecule is unchanged. The legacy `write_pdbqt`/`Mol.to_pdbqt` APIs
  keep their previous bytes and still omit isotope annotations silently.

This profile uses synthetic LIG/A/1 residue identity and has no torsion tree or
`TORSDOF`. Bonds, formal charges, atom maps, stereochemistry and implicit H
counts are outside the output profile. The sidecar only inventories isotope
omission; it is not a complete conversion-loss report or evidence that AutoDock
typing/protonation is chemically correct.

## Verification

The original corpus in
`validation/openbabel_issues/pdbqt-isotopes/cases.json` contains 23 cases and 48
isotope annotations: explicit H/D/T, labelled heavy atoms and ions, aromatic
labels, disconnected components, atom maps, acetophenone, acetophenone-d8 and
acetophenone with eight carbon-13 labels. Expected atomic numbers and isotope
indices were independently checked with RDKit 2026.09.1 using `removeHs=False`.
The corpus is authored for this repository, not copied from an issue attachment.

With independently generated 3D acetophenone inputs containing 17 explicit atoms
and identical coordinates, Open Babel **3.2.1** on this host reproduced the
reported difference with default `obabel -isdf input.sdf -opdbqt`:

| Source labels | Output atoms | `HD` atoms | `TORSDOF` |
| --- | ---: | ---: | ---: |
| None | 9 | 0 | 1 |
| Eight D | 17 | 8 | 2 |
| Eight carbon-13 | 9 | 0 | 1 |

All three commands exited successfully with only `1 molecule converted` on
stderr. The exact SDF inputs, PDBQT outputs and reproduction metadata are retained
in the dated evidence directory. chematic's separate rigid projection retains
17 atoms in both explicit-H acetophenone and d8, assigns all eight C-bound H/D
atoms type `H`, and emits no `TORSDOF`. This is a narrower workflow than Open
Babel's default flexible-ligand preparation.

Native, Python and Node/WASM tests check label refusal, complete sidecars,
source order, fixed-column coordinates/charges and unchanged ordinary output.
Separate fixed-column decoding and the pinned RDKit observer validate the
cross-binding evidence; see `validation/openbabel-3005-2026-10-24/`.

```sh
cargo test -p chematic-mol --test pdbqt_isotope_contract --test pdbqt_checked_contract
python -m pytest crates/chematic-py/tests/test_pdbqt_isotopes.py crates/chematic-py/tests/test_pdbqt_checked.py
node crates/chematic-wasm/tests/pdbqt_isotopes.test.mjs /path/to/node/chematic_wasm.js
```
