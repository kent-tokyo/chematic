# Isotope-aware total formula

[Open Babel #1464](https://github.com/openbabel/openbabel/issues/1464) requests
non-hydrogen isotope labels and total charge in molecular formulas. The
isotope-aware total formula is available as Rust
`chematic_chem::calc_mol_formula_with_isotopes(&mol)`, Python's read-only
`Mol.formula_with_isotopes` property and WASM `MolHandle.formula_with_isotopes()`.

```python
mol = chematic.from_smiles("[13CH3]C(=O)[O-]")
assert mol.formula == "C2H3O2-"
assert mol.formula_with_isotopes == "C[13C]H3O2-"
assert chematic.from_smiles("[2H]O[2H]").formula_with_isotopes == "[2H]2O"
```

The operation includes explicit atom counts, bracket/inferred hydrogens and
net formal charge using the existing formula descriptor's hydrogen/model
convention, including its original-graph fallback. It does not add or remove
atoms from the caller's graph. The existing isotope-blind `formula` output is
unchanged. Total charge already worked in that API and in the measured Open
Babel 3.2.1; the added capability distinguishes isotope-specific counts.

Carbon and hydrogen families come first, then other symbols alphabetically.
Within an element, unspecified atoms precede numeric isotope labels in ascending
mass order. Isotopes use `[13C]`, `[15N]`, `[2H]` and `[3H]`, with a following
count when greater than one. Explicit `[12C]` remains distinct from unspecified
carbon. Mass zero means unspecified; declared nonzero mass numbers are retained
without validating physical isotope existence or abundance. Net charge uses
`+`, `-`, `+2`, `-3`. Wildcards retain `*` or a bracketed isotope-qualified `*`.

This is a total-composition string. It does not encode isotope locations,
connectivity, stereo, atom maps or radical state and is not a molecular identity
key. Existing `parse_formula` does not preserve isotope labels; it is not a
round-trip reader for this output. The core `Molecule::formula_with_isotopes`
helper remains a separate explicit-atom-only view; use the chemistry API above
for hydrogens and net charge. Existing graph/parser budgets still apply.

Authored controls and the open issue snapshot are retained under
`validation/openbabel-1464-2026-10-24/`; the report and raw outputs are in
`validation/results/openbabel-1464-2026-10-24*`. The packet covers 47 compositions,
including isotope order/multiplicity, charge cancellation, explicit/inferred
hydrogens, aromatic heteroatoms and separate wildcard/mass-label boundaries.
Native tests also verify formula invariance when adding ordinary H graph atoms.
Python/WASM tests verify that formula calls leave the graph and legacy output
unchanged.

```sh
python3 scripts/validation/validate_isotope_formula.py \
  --python-root /path/to/candidate-python \
  --rdkit-python /path/to/rdkit-python \
  --wasm-module /path/to/candidate/chematic_wasm.js
```

The Open Babel lane uses the same SMILES with `-ismi -omolreport`, normalized by
an independent formula observer. D/T aliases count as equivalent to `[2H]`/`[3H]`.
Repeated charge signs (`++`, `--`) also normalize to `+2`, `-2`.
Five wildcard/mass-label boundary controls are candidate-only. Expected fields
are authored before execution; RDKit's isotope-separated output provides an
additional cross-check. Source/runtime/fixture hashes are retained. This is
local development evidence, not a published release or a toolkit-wide ranking.

The candidate Python/WASM APIs match all 47 authored compositions. The old
isotope-blind APIs in the same runtime match 12; their complete composition
implementation is verified equal to main `38b4c92c`. This is not a separately
rebuilt main/release binary. In the 42-input common lane, Open Babel 3.2.1
matches 16 compositions: 23 differ only in isotope-specific counts, and three
also undercount ordinary hydrogen. Those hydrogen regressions are separate
from the new isotope capability. No charge differences are introduced as a
new win. RDKit 2026.09.1 matches composition on 46 accepted controls, with one
equivalent ordering difference, and refuses the chemically invalid H–H chain
control. Workers have 30-second/1 MiB combined pipe budgets, with no process
memory cap.
