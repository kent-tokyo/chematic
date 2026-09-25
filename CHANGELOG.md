# Changelog

This file keeps concise, user-visible `chematic` release summaries. Detailed release
history through v1.0.25 is preserved in the
[changelog archive](docs/archive/changelog-through-v1.0.25.md). Validation and
benchmark claims remain scoped to their dated records.

The format follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and public releases follow [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

- **Behavior change (Python/WASM SMARTS):** `smarts_match`, `smarts_find`,
  `Mol.has_substructure`, `Mol.find_matches`, `bulk.substructure_search`,
  `bulk.substructure_match`, and the WASM SMARTS highlight functions now match
  against the molecule's perceived (RDKit-parity) aromatic view. A Kekulé
  benzene matches `c`, and no longer matches `[#6]=[#6]`. Returned atom
  indices still refer to the input molecule. When perception fails, matching
  falls back to the molecule as given. The Rust `find_matches` family is
  unchanged; `chematic_smarts::find_matches_perceived` and
  `has_match_perceived` expose the new semantics (#635).
- `Mol.cip_stereo()` E/Z entries now include `bond_idx` and `bond_atoms`.
  `atom_idx` is the double bond's first atom, not a bond index. The RDKit lane
  compares E/Z labels by bond endpoints (#634).
- `CipMode.ACCURATE` assigns E/Z with the hierarchical-digraph substituent
  ranker (`chematic_cip::SubstituentRanker`), for example an aryl above a
  tert-butyl substituent. `assign_ez_bonds` and `assign_ez_bonds_with_mode`
  return bond-keyed E/Z labels. Legacy mode is unchanged (#634).
- MMFF94 now matches RDKit 2026.03.6 per term on shared coordinates (#637):
  - hydrogen types come from the parent atom's MMFF type (HOCC, HOP, HOS,
    HNSO/HNCS, HNR+ and related, HP);
  - no vdW R* asymmetry correction applies to donor pairs;
  - stretch-bend lookup has no stretch-bend-type-0 retry, and no term is
    added at linear centres;
  - OOP and angle lookups use the full equivalence-level step-down;
  - no torsion term is added for i-j-k-i in 3-rings.
  - Result: on the 265-row set, 262/262 comparable rows are within
    1 kcal/mol (was 230), and the maximum delta falls from 9.87 to
    0.32 kcal/mol. MMFF94 energies, gradients and minimized geometries change
    accordingly.
- Evidence: RDKit rebaseline rows, residual classification, and the #632 long
  relabel audit on v1.0.25-based source; the RDKit per-term MMFF94 oracle
  (`scripts/mmff94_same_explicit_h_energy.py`, schema v3) and
  `scripts/mmff94_atom_type_census.py`.

## [1.0.25] - 2026-09-25

- Added atom-output and source-atom provenance for SMILES, fragments, and Rust
  reaction products.
- Corrected the named RDKit-compatible Python HBA profile and plain SMILES
  aromatic/non-aromatic ring-closure spelling.

## [1.0.24] - 2026-09-24

- Optimized ring perception, RDKit-parity aromatic preparation, and
  SMARTS-existence checks without changing the checked source outputs.

## [1.0.23] - 2026-09-24

- Improved RDKit-defined output agreement for compatible fingerprints, MACCS,
  QED, Murcko scaffolds, and selected descriptors.
- Aligned SMARTS implicit-bond semantics and added `rdkit_tpsa`.

## [1.0.22] - 2026-09-24

- Added derived caches and a versioned RDKit operation matrix with explicit
  equivalence and environment boundaries.

## [1.0.21] - 2026-09-23

- Fixed aromatic-stash E/Z canonical SMILES and added strict documentation-site
  checks.

## Earlier releases

See the [v1.0.0–v1.0.20 archive](docs/archive/changelog-through-v1.0.25.md)
for release summaries and historical details.
