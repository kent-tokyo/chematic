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
| Stretch-bend | `Code/ForceField/MMFF/Params.cpp` | `defaultMMFFStbn` (specific rows), `defaultMMFFDfsb` (default/generic stretch-bend constants by periodic-table-row pair — RDKit's own empirical fallback for stretch-bend) | chematic's `mmff94_stbn` has a partial fallback chain already (angle-type→0, then generic `(0,0,type_j,0)`); RDKit's `defaultMMFFDfsb` is periodic-row-keyed, not identical in shape. **Priority 2 update (issue #227, 2026-08-05)**: `defaultMMFFDfsb` (29 rows, frozen copy `scripts/mmff94_provenance/rdkit_defaultMMFFDfsb.txt`, extracted programmatically not hand-transcribed) is now ported as a **diagnostic-only** classifier in `mmff94_term_coverage_audit.rs` (`dfsb_default_resolvable` field) — NOT wired into chematic-ff's production energy/gate path. Algorithm verified against `Code/GraphMol/ForceFieldHelpers/MMFF/AtomTyper.cpp` at the pinned commit: `MMFFMolProperties::getMMFFStretchBendParams` (lines ~3566-3612) tries `MMFFStbnCollection::getMMFFStbnParams` (specific/generic by MMFF *type*, single exact lookup after I/K canonicalization — no equivalence-class step) first, and **only on failure** falls back to `MMFFDfsbCollection::getMMFFDfsbParams(periodicRow(atom1), periodicRow(atom2), periodicRow(atom3))`, keyed by periodic-table row via `getPeriodicTableRow` (`AtomTyper.cpp:251-264`: atomic number 1-2→row 0, 3-10→row 1, 11-18→row 2, 19-36→row 3, 37-54→row 4, else 0), canonicalized so `min(row1,row3) <= max(row1,row3)`. Confirmed **no equivalence-class (`eqLevel`) step exists anywhere in the stretch-bend resolution path** — `eqLevel` is used only by RDKit's angle/torsion/OOP fallback functions (`AtomTyper.cpp:527,552,743,768,862`), not stretch-bend — so a `equivalence_fallback_resolvable` bucket does not apply to stretch-bend under RDKit's real algorithm, and building one for it would misrepresent RDKit's own behavior, not just be incomplete. **Priority 2B update (issue #227, 2026-08-06)**: the same table is now ported into PRODUCTION `chematic_ff::mmff94_stbn` (`crates/chematic-ff/src/mmff94_energy/oop_stbn.rs`), unconditional, not gated — the diagnostic-only `dfsb_default_resolvable` field above was replaced by a `dfsb_resolved: bool` field (`mmff94_term_coverage_audit.rs`) that keeps the type-only diagnostic (`present_at_different_classification`, via the newly-split-out `mmff94_stbn_type_only`) and the final production-resolution question separate on purpose: of the 2,107 instances the type-only lookup misses, 1,680 are genuine table gaps that Dfsb closing matches RDKit's real behavior for, but 427 are routing-bug candidates (a real, correctly-typed parameter exists at a *different* classification code) that Dfsb *also* happens to resolve — coverage achieved, but chematic uses RDKit's generic default instead of the correctly-routed specific parameter for those 427, a parameter-selection-parity gap, not fixed by this port (see `validation/results/mmff94_coverage_227_term_audit_summary.json`'s `stretch_bend_dfsb_resolution` for the exact split). **Priority 2C update (issue #227, 2026-08-09)**: root-caused the 427, diagnostic-only, in `mmff94_stbn_equivalence_diagnostic_227.rs` — chematic's production stretch-bend code uses the **angle type** (0-8, `angle_type_for`) directly as the `MMFF94_STBN` key, but RDKit's real `getMMFFStretchBendParams` (`AtomTyper.cpp:3566-3612`) computes a *distinct*, finer-grained **stretch-bend type** (0-11) via `getMMFFStretchBendType(angleType, bondType1, bondType2)` (`AtomTyper.cpp:2480-2508`, ported verbatim) and `MMFFStbnCollection::getMMFFStbnParams`'s I/K canonicalization (`Params.h:601-663`: swap iff `iAtomType > kAtomType`, or tie-break on raw `bondType1 < bondType2` when equal) — not an `eqLevel` gap, confirming the Priority 2 finding above still holds. Also ported and cross-checked: `getMMFFAngleType`'s real ring-offset formula (`AtomTyper.cpp:2412-2447`: `angleType = size; if (bondTypeSum) angleType += bondTypeSum + size - 2`) disagrees with `angle_type_for`'s own ring-offset table for `bt_sum=2` (3-ring) and `bt_sum∈{1,2}` (4-ring) — a real, independent, second bug, measured LATENT on the 265-molecule corpus (0/113 reachable ring-embedded angle triples exercise the diverging branches); and `isAngleInRingOfSize3or4` (`AtomTyper.cpp:357-395`), which is local bond-adjacency, NOT SSSR-based (0/10,107 triples disagreed with chematic's SSSR-based ring check on this corpus, but the two are not the same algorithm). Live-oracle cross-check via `MMFFMolProperties.GetMMFFStretchBendParams`/`GetMMFFBondStretchParams` (`scripts/mmff94_stbn_oracle_validate_227.py`, same pinned RDKit build) confirms the ported formula matches RDKit exactly on 255/427 candidates (100%, zero exceptions) where chematic's own bond-order/aromaticity perception agrees with RDKit's; the remaining 172/427 are confounded by a separate, pre-existing aromaticity-perception gap (chematic trusts certain lowercase-aromatic input SMILES RDKit's sanitizer kekulizes instead — same mechanism also explains all 172 (molecule, triple) instances shared between `Angle`'s 277 and `StretchBend`'s 427 populations, set-identity verified). See the PR for issue #227 for the full breakdown. |
| Out-of-plane (Wilson angle) | `Code/ForceField/MMFF/Params.cpp` | `defaultMMFFOop` (MMFF94), `defaultMMFFsOop` (MMFF94s variant, not used by chematic) | |
| Torsion | `Code/ForceField/MMFF/Params.cpp` | `defaultMMFFTor` (MMFF94), `defaultMMFFsTor` (MMFF94s variant) | |
| van der Waals | `Code/ForceField/MMFF/Params.cpp` | `defaultMMFFVdW` | |
| Charges (partial bond charge increments) | `Code/ForceField/MMFF/Params.cpp` | `defaultMMFFPBCI`, `defaultMMFFChg` | `defaultMMFFPBCI` is already the cited source for chematic's existing `pbci_for` table (pre-dates this PR). |
| Aromaticity perception feeding MMFF typing | `Code/GraphMol/Aromaticity.cpp` | `setMMFFAromaticity` (module-level function, not a `MolOps` member despite earlier notes in this file placing it there) | Priority 1A (issue #227): ported as `compute_mmff94_aromatic_view` in `mmff94_numeric.rs` — a **partial, behaviorally-calibrated** port, not a full one: every rule is a direct, line-cited port (ring-by-ring pi-electron counting at lines ~955-1035, the exocyclic-double-bond/NOS lone-pair-bonus rules, the multi-pass resolution loop) except the hybridization gate at line 1023 (`atom->getHybridization() != Atom::SP2`), approximated as `total_degree(atom) > 3` since chematic has no general hybridization-inference engine to port this faithfully. Measured gap on the 265-molecule Wave 1 corpus (`scripts/mmff94_hybridization_gate_gap_227_report.py`): 4,128/4,172 (98.9%) ring C/N atoms same decision as RDKit, 44 where the approximation under-triggers (misses a real pyramidal-SP3 ring N), 0 where it over-triggers, 0 unclassified — see `validation/results/mmff94_hybridization_gate_gap_227_report.txt`. RDKit's own general aromaticity model (distinct from both this MMFF-specific one and from chematic's `chematic_perception::apply_aromaticity`) is not relevant to MMFF typing and is out of scope here. |

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

`rdkit_defaultMMFFDfsb.txt` was extracted from the same pinned `Params.cpp` (the
`defaultMMFFDfsb` string literal, lines ~4894-4932 at this commit) via a small
one-off Python script that regex-extracts every quoted string segment,
concatenates them, and unescapes `\t`/`\n` — not hand-transcribed. Re-run the
same extraction against a fresh `/tmp/Params.cpp` if the pin is ever bumped;
`getPeriodicTableRow`'s 5-way atomic-number bucketing (ported into
`mmff94_term_coverage_audit.rs` as `rdkit_periodic_table_row`) is from
`Code/GraphMol/ForceFieldHelpers/MMFF/AtomTyper.cpp` at the same pinned
commit, lines 251-264.
