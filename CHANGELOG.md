# Changelog

This file summarizes user-visible changes. Exact validation and performance
claims belong to dated records because they depend on the artifact, corpus,
comparator, options, and host.

Full notes through v1.0.41 are preserved in the
[changelog archive](docs/archive/changelog-through-v1.0.41.md). Later release
details remain available in the immutable tags and
[GitHub Releases](https://github.com/kent-tokyo/chematic/releases).

The format follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and public releases follow [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Added

- Added opt-in proof-backed removal of one redundant tetrahedral declaration in
  Rust/Python/WASM. An anchored graph witness preserves all other configurations
  and source labels; failed proofs retain the declaration, and unsupported inputs
  or exhausted budgets return errors. Parsing/writing defaults remain unchanged.

### Performance

- Parse complete buffered SDF records directly, reserving owned record storage
  only for the existing line-reader fallback. Preserve diagnostic/lightweight
  results, limits, errors, metadata and logical input positions.
- Append SD data fields directly and skip stereo-loss collection for plain
  V2000 output callers that discard it. Reporting APIs retain every loss field;
  stereo depiction and serialized bytes remain unchanged. See the
  [bounded source comparison](benchmarks/2026-10-16-openbabel-sdf-loss-work.md).

## [1.1.5] - 2026-10-16

### Added

- Added checked molecular and materials interchange for HyperChem HIN,
  wwPDB CCD, ChemDoodle JSON, DMol3, MOPAC, FHI-aims, ORCA-AICCM EXYZ,
  RDFile, TURBOMOLE, GRO, Gaussian-style Z-matrices, mmCIF-to-PDB, and
  PQR-to-PDBQT. Readers retain source metadata where declared and refuse
  unsupported or lossy projections unless the caller explicitly opts in.
- Added indexed and bounded SDF selection, including slices, integer indices,
  strict boolean masks, source-indexed batches, and V3000 CTAB dispatch while
  retaining SD fields and malformed-record positions.
- Added bounded graph utilities for heavy-atom distance matrices, RDKit-style
  path enumeration, atom-radius environments, adjacency matrices, shortest
  paths, and checked atom permutations.
- Added RDKit 2026.09.1-compatible Layered fingerprint options and
  target-density fingerprint folding with explicit unmet-target reporting.
- Added checked stereo-group construction, inversion, alternative expansion,
  and conservative meso reduction for supported tetrahedral and geometric
  members. Unsupported axes continue to fail explicitly.
- Added conservative query-hydrogen folding and chirality-aware compiled
  Python SMARTS matching.
- Added checked single-model PDB selection, general mmCIF residue grouping,
  and molecular-template projection while preserving source identities and
  reporting metadata loss.
- Added an opt-in checked RDKit-model 3D distance matrix to Rust, Python, and
  WASM with allocation, finite-value, weighted-dummy, and extreme-distance
  safeguards.
- Corrected native stereoisomer candidate selection so identical terminal
  ligands and CH2 sites do not consume the bounded enumeration budget.

### Performance

- Reduced repeated ring and aromaticity work by reusing cyclic-component
  labels, normalized descriptor views, and preallocated traversal storage.
- Reduced TPSA and E/Z setup for molecules that cannot contribute to the
  requested result, without changing the declared chemistry contract.
- Reduced V3000 parsing and writing overhead through token reuse, declared
  capacity reservations, fast common numeric forms, and direct zero-map
  serialization. Output-retention tests protect exact bytes where promised.

### Fixed

- Restored the v1.1.3 exhaustive-match surface of public `CmlError` and
  `CifError` enums while preserving explicit failures for invalid inputs.
- Added a targeted `cargo-semver-checks` lane and a downstream compile fixture
  for the established v1.1.3 `chematic-mol` contract.
- Pinned known RDKit 2026.03.6/2026.09.1 `CalcNumRings` disagreements as a
  versioned contract instead of treating one comparator as universal truth.
- Updated transitive `rustls` to 0.23.46 for RUSTSEC-2026-0347,
  RUSTSEC-2026-0348, and RUSTSEC-2026-0349.

### Documentation

- Synchronized the three README files, documentation home, roadmap,
  validation summary, open-work ledger, and Trust Release rules with the
  v1.1.5 release boundary.
- Replaced the obsolete v1.0.8 verification diary with a short historical
  pointer to the current validation and benchmark records.

## [1.1.4] - 2026-10-14

### Added

- Added bounded cation-resonance enumeration, screening-enrichment metrics,
  checked CIF output and protein projections, and pre-smoothing internal
  constraints for the bounded 3D pipeline.
- Expanded RDKit-compatible name, typed property, conformer RMSD,
  fingerprint, atom-masking, and dense NumPy fingerprint APIs.
- Added local protein-projection and selected interchange improvements while
  retaining explicit loss reports and input limits.

### Changed

- Improved molecular-weight, Morgan, Labute ASA, MOL/V3000, MOL2, CML, CIF,
  and selected 3D compatibility paths under their recorded test boundaries.
- Restored exhaustive public-enum compatibility for `Chirality`; later CML
  and CIF error compatibility work remains in Unreleased source.

### Verification boundary

- The published release is v1.1.4. Its complete cross-platform artifact packet
  was still open at documentation cutoff; source results and older package
  packets are not relabelled as v1.1.4 package evidence.
- See the [validation report](docs/validation.md) and
  [benchmark index](benchmarks/README.md) for exact denominators and artifacts.

## [1.1.3] - 2026-10-13

- Unified QED component evaluation on one normalized molecular view, improved
  selected V3000 stereo handling, and added persistent MMFF/UFF evaluation,
  pickle support, and lazy SDF reading.
- Full details: [v1.1.3 release](https://github.com/kent-tokyo/chematic/releases/tag/v1.1.3).

## [1.1.2] - 2026-10-13

- Continued the bounded reaction, standardization, R-group, descriptor, and
  3D compatibility program without broadening unsupported chemistry claims.
- Full details: [v1.1.2 release](https://github.com/kent-tokyo/chematic/releases/tag/v1.1.2).

## [1.1.1] - 2026-10-12

- Improved meaning-preserving file I/O, WASM reaction and labelled R-group
  diagnostics, and force-field recovery behavior.
- Full details: [v1.1.1 release](https://github.com/kent-tokyo/chematic/releases/tag/v1.1.1).

## [1.1.0] - 2026-10-11

- Established the 1.1 release line around typed failures, cross-binding
  contracts, broader interchange, and reproducible validation records.
- Full details: [v1.1.0 release](https://github.com/kent-tokyo/chematic/releases/tag/v1.1.0).

## [1.0.42] - 2026-10-10

- Improved bounded Morgan and MMFF performance while retaining the versioned
  chemistry and output gates used by the published comparison packet.
- Full details: [v1.0.42 release](https://github.com/kent-tokyo/chematic/releases/tag/v1.0.42).

## Older releases

Release notes through v1.0.41 are in the
[versioned archive](docs/archive/changelog-through-v1.0.41.md). The archive is
historical evidence and does not establish current behavior.
