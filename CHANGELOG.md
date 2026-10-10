# Changelog

This file keeps concise, user-visible release summaries. Full release details
through v1.0.41 are preserved in the
[changelog archive](docs/archive/changelog-through-v1.0.41.md); validation and
performance claims remain scoped to their dated records.

The format follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and public releases follow [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Changed

- Expanded Rust coverage measurement to exercise native InChI reconstruction,
  canonical SMILES, stereochemical MOL readers, reaction caches, MMFF94 terms,
  and cross-binding contracts. The report now excludes the ignored manual 3D
  differential harness from the production-code denominator.
- Reduced RDKit-compatible Morgan allocation for folded-bit callers, reused
  parser-proven acyclic/single-cycle ring facts, and removed duplicate MMFF
  force-field preparation. Existing detailed fingerprint and chemistry results
  remain covered by the same output gates.
- Reused parser-known SMILES cycle rank for SSSR counts and removed discarded
  S/P typing and per-atom output allocation from whole-molecule RDKit TPSA and
  Labute ASA calculations. ChEMBL/NCI agreement remains unchanged.
- Upgraded the RDKit/COSMolKit benchmark to retain raw alternating-order
  samples, report operation-only and parse-inclusive pipelines separately, and
  calculate paired bootstrap confidence intervals.

### Evidence boundary

- A clean-commit macOS arm64 candidate beats RDKit 2026.03.1 and COSMolKit
  0.5.0rc15 in all 21 measured blocks for the seven recorded parse-inclusive
  pipelines. Five operation-only lanes remain slower, and registry artifacts
  and additional hosts have not been rerun. See the
  [dated record](benchmarks/2026-10-10-rdkit-cosmolkit-performance-candidate.md).
- A second candidate trims direct ring/descriptor work while retaining exact
  ChEMBL/NCI outputs. TPSA and Labute pipelines still win, but their isolated
  calls remain slower than both comparators; see the
  [follow-up record](benchmarks/2026-10-10-direct-descriptor-fastpaths-candidate.md).

### Documentation

- Added a crates.io-first `chematic-mcp` quick start with direct stdio client
  configuration (#785).
- Condensed the roadmap, open-work ledger, Trust Release rules, and public
  documentation entry points. Detailed history remains versioned in Git and
  the documentation archive.

## [1.0.41] - 2026-10-10

### Added

- Added RDKit 2026.03.1-compatible tautomer enumeration and scoring, CXSMILES,
  random and fragment SMILES, distance matrices, sparse fingerprints,
  chemistry-problem detection, reaction writers, and related compatibility
  helpers.
- Expanded the RDKit/COSMolKit comparison harness to 161 operations on ChEMBL
  5k and RDKit.js 10k, plus a 4,072-row unusual-SMILES stress corpus.

### Changed

- RDKit-compatible descriptors, fingerprints, SMARTS, and reactions use the
  RDKit-model molecular view for the bounded structures where chematic's native
  aromaticity or cleanup model differs.
- Reaction compatibility now handles grouped components, RDKit-style cleanup,
  carried atom maps, explicit-hydrogen removal, and rejected products.
- Dummy atoms, isotope ordering, formulas, and elements 113–118 now follow the
  documented RDKit-compatible semantics.

### Fixed

- Normalized zero TPSA to `0.0`, preserved tetrahedral orientation when
  explicit hydrogen is removed, retained exocyclic E/Z through aromatic
  conversion, and refused contradictory RDKit-compatible 2D MOL output.

### Evidence boundary

- On the recorded source comparison, chematic matches RDKit at least as often
  as COSMolKit 0.5.0rc15 on all 161 covered operations over both main corpora.
  The stress corpus retains one weaker native canonical-SMILES round trip.
  These are source-build results, not universal parity or a substitute for
  published-package reruns. See
  [the dated record](benchmarks/2026-10-10-cosmolkit-050.md).

## [1.0.40] - 2026-10-09

- Added the preceding RDKit-compatible writer, reader, stereoisomer,
  alignment, MolHash, force-field, and fixed-seed ETKDG surfaces. Full details
  remain in the archive below.

## Earlier releases

- [Full history through v1.0.41](docs/archive/changelog-through-v1.0.41.md)
- [Earlier compact archive through v1.0.25](docs/archive/changelog-through-v1.0.25.md)
