# Hydrogen-count accounting

[Open Babel #3008](https://github.com/openbabel/openbabel/issues/3008) reports
incorrect hydrogen counts in formulas for explicit-H and disconnected
SMILES. [#2624](https://github.com/openbabel/openbabel/issues/2624) reports
bracket hydrogens disappearing from formulas for disconnected components.
Both issues remain open as inspected on 2026-10-16; local Open Babel 3.2.1
reproduces these failure classes.

chematic's molecular formulas already retain these hydrogens. Investigation
found a separate defect in `num_hydrogens`: `[NH4+]` returned 8, while
`[H][O][H]` returned 0. The implementation added the stored bracket-H count
twice and did not count H atom nodes. This change fixes that descriptor and
the WASM atom-information field without changing the formula implementation,
SMILES serialization, valence model or hydrogen transformation algorithms.

## Molecular total

Rust: `chematic_chem::num_hydrogens(&mol)`.
Python: `mol.num_hydrogens`.
WASM: the new `num_hydrogens(mol)` export.

```python
import chematic

assert chematic.from_smiles("[NH4+]").num_hydrogens == 4
assert chematic.from_smiles("[H][O][H]").num_hydrogens == 2
assert chematic.from_smiles("[OH2].[OH2]").num_hydrogens == 4
```

The count includes:

- Each real H atom node once, including isotopic, charged, isolated,
  bridging and H-H-bonded hydrogen. Bond degree does not multiply this count.
- Attached bracket H once, or valence-inferred H when no bracket count is
  stored. `implicit_hcount` already returns the stored bracket count.
- A wildcard's explicitly stored attached H. Its backing element field
  does not make the wildcard itself a real H atom.

All components contribute independently; atom order and component order do
not affect the total. An empty graph has zero hydrogens. Isotopes count as
hydrogen nuclei, not isotope-specific molecular mass.

## Per-atom attached H in WASM

`get_atom_info(mol, index).totalHydrogens` now counts attached bracket/inferred
H plus explicit H neighbors, including isotopes. It excludes the queried
atom itself when that atom is H. The JSON shape and out-of-range `null`
behavior are unchanged.

For `[H][O][H]`, the per-atom counts are `[0, 2, 0]`, while the molecular
total is 2. For `[Al][H-][Al]`, the per-atom counts are `[1, 0, 1]` but the
molecular total is 1: both Al atoms share the same hydride. Do not sum
per-atom attached-H counts to obtain the molecule's hydrogen total.

## Verification boundary

The original input table is
`validation/openbabel_issues/hydrogen-counts/cases.tsv`. Rust, Python and real
Node/WASM regressions include bracket H, explicit H graphs, ions, isotopes,
wildcards, aromatic molecules, component/atom-order changes and selected
ordinary add-H/remove-H representation controls. These do not extend the
supported domains of the hydrogen transformation APIs.

`scripts/validation/validate_hydrogen_count_contract.py` uses RDKit 2026.09.1 as an
independent graph reference. It retains H nodes (`removeHs=False`) and uses
an unsanitized graph plus a non-strict property cache for unusual examples.
The molecular oracle counts real H nodes and each atom's attached bracket/
inferred H; the per-atom oracle is `GetTotalNumHs(includeNeighbors=True)`.
A generated cohort varies compact/bracket/explicit-node representations and
atom order; a separate cohort pairs disconnected components. Full input and
output records, source/runtime hashes and retained pre-fix source counts are
stored under `validation/results/openbabel-hydrogen-counts-2026-10-16*`.

Open Babel formula-H diagnostics and SMILES-output-H diagnostics are recorded
separately. A correct formula does not prove that serialization retains H
nodes or charges. The chematic formula is an unchanged control, not a newly
implemented fix. These checks cover hydrogen accounting on a source
candidate, not chemical validity, mass accuracy, released wheels or general
Open Babel/RDKit compatibility.

The 2026-10-16 run covers 280 cases (202 distinct SMILES). Python and real
Node/WASM total counts match all 280, and WASM attached-H fields match all
1,132 atom entries. The retained pre-fix source runtime miscounts 224 cases;
its formula-H control already matches all 280. The source candidate also
matches all 280 after native SMILES round trips.

On the 36 original inputs, installed Open Babel 3.2.1 supplies 34 formulas:
eight have a different aggregate H count. `*` and `[*H2]` return empty formula
fields, recorded separately. All 36 serialized outputs are checked for H
count; two differ. Those serialization differences are diagnostics, not
new writer fixes in this change.

Local checks pass 11 Rust tests (four new regressions, four existing shared
contract tests and three existing H-counter unit tests), 46 Python tests,
the Node/WASM regression, Clippy with warnings denied, formatting, binding/
manifest checks, refactor budget and strict MkDocs.
