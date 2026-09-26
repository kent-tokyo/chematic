# Changelog

This file keeps concise, user-visible `chematic` release summaries. Detailed release
history through v1.0.25 is preserved in the
[changelog archive](docs/archive/changelog-through-v1.0.25.md). Validation and
benchmark claims remain scoped to their dated records.

The format follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and public releases follow [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

## [1.0.27] - 2026-09-26

- Added Rust `Molecule::set_tag` / `atom_tag` for caller-managed atom labels
  preserved by molecule clone, core graph edits, reaction apply, fragments,
  and aromaticity perception. Private, lazily allocated tag storage preserves
  existing `Atom` struct literals and equality. Tags do not affect SMILES or
  canonicalization; write/parse requires explicit atom-order remapping. Labels
  are `1..=u16::MAX`, need not be unique, and `None` / `Some(0)` clear them.
  Tag changes invalidate cached molecular views.
- Restored the WASM torsion-scan demo API and added browser smoke coverage for
  its typed response; this does not expand the supported 3D chemistry domain.
- Strengthened the RDKit 2026.03.6 SMARTS residual checker to account for
  every classified query/target cell. The remaining 200 of 310,000 cells are
  documented compatibility boundaries, not exact parity.
- Published a separately pinned v1.0.26-wheel MMFF94 quality packet with all
  265 input rows retained and independent geometry/stereo scoring. This is
  historical evidence for v1.0.26, not a v1.0.27 remeasurement; broader A6
  typing, convergence, and conformer-quality gates remain open.

## [1.0.26] - 2026-09-25

- Improved RDKit-compatible SMARTS matching, aromatic/ring perception, and
  selected fingerprint hot paths. The dated source differential preserves the
  measured output boundary; shared-VM timing is not a package or universal
  performance claim.
- **Behavior change (Python/WASM SMARTS):** public SMARTS APIs now match a
  perceived RDKit-parity aromatic view. A Kekulé benzene matches `c`, and no
  longer matches `[#6]=[#6]`; returned atom indices still refer to the input.
  Rust's default matcher is unchanged.
- Made accurate E/Z CIP output bond-keyed, including bond endpoints, and
  improved same-coordinate MMFF94 per-term agreement with RDKit 2026.03.6.
  These are scoped source-evidence improvements, not complete parity,
  convergence, conformer-quality, or published-package claims.

Detailed inputs, results, and limits are in [validation](docs/validation.md)
and the [benchmark index](benchmarks/README.md).

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
