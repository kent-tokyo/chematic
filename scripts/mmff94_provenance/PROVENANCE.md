# MMFF94 source provenance (Phase 1B-0, issue #227)

Pinned, not moving-`main`. All RDKit citations below resolve to one fixed
commit; re-fetch commands are given so this can be regenerated if the pin is
ever bumped, but no code in this repository reads RDKit's git history at
build or run time.

## Pinned RDKit revision

- **pip package installed for this audit**: `rdkit==2026.3.3` (`python -c
  "import rdkit; print(rdkit.__version__)"` → `2026.03.3`)
- **Matching RDKit git tag**: `Release_2026_03_3`
- **Commit SHA**: `e74e7b0a5a2fc4e7f77c04ec26a61d4b8edbf22f`
- **Resolved via**: `gh api repos/rdkit/rdkit/tags --paginate` (GitHub's tags
  API, not a local moving clone) — the tag→SHA mapping is immutable once
  published, so this is a fixed reference, not a moving target.
- **License**: BSD-3-Clause (RDKit `license.txt` at this commit) — reuse of
  parameter table data with attribution is permitted; this repo already has
  one precedent (`mmff94_numeric.rs`'s `defaultMMFFPBCI`-derived PBCI table,
  cited the same way before this PR).

## Source files and tables, by MMFF94 subsystem

All paths are relative to the RDKit repo root at commit
`e74e7b0a5a2fc4e7f77c04ec26a61d4b8edbf22f`. Raw URL pattern:
`https://raw.githubusercontent.com/rdkit/rdkit/e74e7b0a5a2fc4e7f77c04ec26a61d4b8edbf22f/<path>`.

| Subsystem | File | Table / function | Notes |
|---|---|---|---|
| Atom typing rules (the actual classification logic — aromatic 5-ring/6-ring alpha/beta detection, N-oxide, imidazolium, etc.) | `Code/GraphMol/ForceFieldHelpers/MMFF/AtomTyper.cpp` | `MMFFMolProperties::setMMFFHeavyAtomType`, `RingMembershipSize`, `isAtomNOxide`, `isRingAromatic` | 3726 lines total; this repo ports only the aromatic C/N/O/S 5-ring/6-ring block (lines ~503–800 at the pinned commit), not the full aliphatic/metal/halogen switch — see the PR body for exact scope. |
| Numeric type definitions + equivalence-class levels | `Code/ForceField/MMFF/Params.cpp` | `defaultMMFFDef` string table (parsed by `MMFFDefCollection`) | Format: `SYMBOL\tTYPE\tEQ_LEVEL2\tEQ_LEVEL3\tEQ_LEVEL4\tEQ_LEVEL5\tDESCRIPTION...`; lines starting with `*` are secondary/alias symbols, skipped by RDKit's own parser (`inLine[0] != '*'`). Frozen copy: `scripts/mmff94_provenance/rdkit_defaultMMFFDef.txt` (extracted verbatim, C-string-escapes decoded, from the pinned commit — not hand-transcribed). |
| Atom-type intrinsic properties (atomic number, coordination, valence, lone-pair/multiple-bond/aromaticity/linear/sbmb flags) | `Code/ForceField/MMFF/Params.cpp` | `defaultMMFFProp` string table (parsed by `MMFFPropCollection`) | Format: `atype\taspec(atomic#)\tcrd\tval\tpilp\tmltb\tarom\tlin\tsbmb`. Frozen copy: `scripts/mmff94_provenance/rdkit_defaultMMFFProp.txt`. This is the field this PR's semantic-compatibility gate is built on — `aspec` is ground truth for "which element is this numeric type allowed to represent." |
| Bond types / bond-stretch parameters | `Code/ForceField/MMFF/Params.cpp` | `defaultMMFFBond` (params), `defaultMMFFBndk` (empirical-rule Badger's-rule-like constants), `defaultMMFFHerschbachLaurie`, `defaultMMFFCovRadPauEle` (covalent radii + Pauling electronegativities feeding the empirical bond rule) | Confirms MMFF94 *does* define a real empirical bond-stretch fallback (Halgren Part V) — chematic has none. Not implemented in this PR (see bond-fallback classification below); cited for the future-PR decision only. |
| Angle types / angle-bend parameters | `Code/ForceField/MMFF/Params.cpp` | `defaultMMFFAngleData` | Matches chematic's `MMFF94_ANGLE_ENERGY` table structurally (angle_type, i, j, k, ka, theta0). |
| Stretch-bend | `Code/ForceField/MMFF/Params.cpp` | `defaultMMFFStbn` (specific rows), `defaultMMFFDfsb` (default/generic stretch-bend constants by periodic-table-row pair — RDKit's own empirical fallback for stretch-bend) | chematic's `mmff94_stbn` has a partial fallback chain already (angle-type→0, then generic `(0,0,type_j,0)`); RDKit's `defaultMMFFDfsb` is periodic-row-keyed, not identical in shape — noted as a residual difference, not implemented here. |
| Out-of-plane (Wilson angle) | `Code/ForceField/MMFF/Params.cpp` | `defaultMMFFOop` (MMFF94), `defaultMMFFsOop` (MMFF94s variant, not used by chematic) | |
| Torsion | `Code/ForceField/MMFF/Params.cpp` | `defaultMMFFTor` (MMFF94), `defaultMMFFsTor` (MMFF94s variant) | |
| van der Waals | `Code/ForceField/MMFF/Params.cpp` | `defaultMMFFVdW` | |
| Charges (partial bond charge increments) | `Code/ForceField/MMFF/Params.cpp` | `defaultMMFFPBCI`, `defaultMMFFChg` | `defaultMMFFPBCI` is already the cited source for chematic's existing `pbci_for` table (pre-dates this PR). |
| Aromaticity perception feeding MMFF typing | `Code/GraphMol/MolOps.cpp` (not vendored here) | `MolOps::setMMFFAromaticity` | RDKit runs its own aromaticity model before MMFF typing, distinct from chematic's `chematic_perception::apply_aromaticity`. Divergences here show up as `aromaticity_perception_divergence` in the atom-typing audit, **not** as numeric-typing bugs — see the audit JSONL's exclusion bucket. |

## Halgren primary literature (secondary/theoretical cross-reference, not the implementation source)

- T. A. Halgren, "Merck Molecular Force Field. I. Basis, Form, Scope,
  Parameterization, and Performance of MMFF94," *J. Comput. Chem.* **17**,
  490–519 (1996).
- T. A. Halgren, "MMFF VI. MMFF94s option for energy minimization studies,"
  *J. Comput. Chem.* **20**, 720–729 (1999) — MMFF94s variant (RDKit's
  `defaultMMFFsOop`/`defaultMMFFsTor`), not used by chematic; noted for
  completeness only.
- T. A. Halgren, "MMFF VII. Characterization of MMFF94, MMFF94s...," *J.
  Comput. Chem.* **20**, 730–748 (1999).

This repo implements against RDKit's transcription of the Halgren tables
(RDKit's own `Params.cpp` comments cite "Copyright (c) Merck and Co., Inc.,
1994, 1995, 1996" directly on the table data), not against a fresh reading
of the papers — RDKit is the practical interoperability target (this whole
program's benchmark is "close the gap to RDKit"), and its numbering is
independently verified against the Halgren type names it carries in the
same table (`CB`, `C5A`, `NPYD`, etc., matching Halgren's own Table I
symbols).

## Regeneration

```bash
SHA=e74e7b0a5a2fc4e7f77c04ec26a61d4b8edbf22f
curl -sL "https://raw.githubusercontent.com/rdkit/rdkit/$SHA/Code/ForceField/MMFF/Params.cpp" -o /tmp/Params.cpp
# then re-run scripts/gen_mmff94_numeric_type_registry.py --extract /tmp/Params.cpp
# to refresh scripts/mmff94_provenance/rdkit_default{MMFFDef,MMFFProp}.txt
python3 scripts/gen_mmff94_numeric_type_registry.py
```
