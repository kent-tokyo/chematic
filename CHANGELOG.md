# Changelog

This file records public releases and the current unreleased changes to
`chematic`. Detailed development notes are retained in
[`docs/archive/detailed-development-history.md`](docs/archive/detailed-development-history.md).

The format follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and public releases follow [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

No unreleased changes.

## [1.0.7] - 2026-09-05

### Performance

- Optimized canonical rank normalization, SMILES ring-label handling, and
  V2000/SDF parsing and serialization while preserving output bytes.
- Added a resumable alternating A/B benchmark with exact-output checks. On
  the recorded v1.0.6-to-v1.0.7 source comparison, canonical SMILES improved
  1.176x, SDF graph/property read 1.180x, and reusable-buffer SDF write
  1.419x at the median of seven paired ratios. SMILES parse improved 1.034x
  and remains below the 1.10x target.

### Validation

- Added bounded ring-label and stereo-partner regressions, direct-append SDF
  streaming boundary tests, and benchmark-runner interruption/resume tests.
- Recorded exact output parity on two 5,000-input SMILES corpora and the
  365-record SDF fixture. The measurements remain local-source evidence and
  are not universal cross-platform guarantees.
- Updated the release, security, benchmark, and roadmap documents to the
  v1.0.7 publication boundary.

## [1.0.6] - 2026-09-05

### Fixed

- Added descriptor provenance metadata and a shared descriptor contract fixture
  consumed by Rust, Python, and Node/WASM tests.
- Added a reproducible 4,999-molecule core descriptor parity lane and a
  2,000-pass SDF/MOL/XYZ streaming benchmark record. The report keeps the
  file-backed streaming and Python block-parser boundaries explicit.

- Prevented canonical SMILES E/Z carrier selection from choosing the
  suppressed (ring-closure close-side) occurrence. Ring-adjacent stereobonds
  now retain their directional marker when the canonical output is reparsed.
  The coupled-E/Z fail-closed boundary remains unchanged.
- Added a bounded default-Hückel fused/non-alternant fallback for all-carbon
  odd/odd envelopes (including azulene), while retaining the separate
  `RdkitLike` model and its holdout gate for broader parity work.
- Made the accurate CIP held-out suite assert zero wrong confident labels and
  zero regressions for the non-phosphorus scope: 140/140 cases currently
  match the modern RDKit oracle. The 15 phosphorus cases remain separately
  reported and now fail closed as `OracleUnstable`; no confident phosphorus
  label is emitted until a representation-independent oracle exists.

- Strengthened fused/non-alternant aromaticity handling and added explicit
  residual manifests for canonical E/Z and MMFF94 cycle diagnostics.

### Added

- Added descriptor provenance and a shared Rust/Python/Node/WASM descriptor
  contract fixture with a reproducible 4,999-molecule parity lane.
- Added versioned machine-readable release metadata, schema validation, and a
  tag-driven GitHub Release asset for downstream integrations.

## [1.0.5] - 2026-09-05

### Fixed

- Closed the legacy-coordinate UFF stereo-rescue residual tracked by #210.
  The five named ibuprofen, naproxen, testosterone, cholesterol, and
  atorvastatin cases now return finite, stereo-satisfied results.

### Performance

- Reduced canonical E/Z setup work; the recorded 5,000-molecule source-level
  median improved from 120.74 ms to 108.28 ms (1.115x).
- Reduced file-backed SDF parsing overhead; the recorded 365-record median
  improved 1.130x.
- Added reusable V2000/SDF output buffers and allocation-free fixed-width
  integer emission. Returned-`String` serialization improved 3.042x and the
  reusable-buffer path improved 3.139x over the pre-change v1.0.4 median.

### Validation

- Added deterministic canonical benchmark digests and regressions for SDF
  byte identity, fixed-width formatting, buffer reuse, and invalid UTF-8.
  Scope and reproduction details are in
  [`benchmarks/2026-09-05-hot-path-follow-up.md`](benchmarks/2026-09-05-hot-path-follow-up.md).

## [1.0.4] - 2026-09-04

### Added

- Added immutable `PreparedReaction` templates for repeated reaction matching
  and application, with reusable ring perception.
- Added typed, loss-aware reaction documents and RXN V2000 adapters that
  preserve agents, coefficients, conditions, provenance, mappings, and steps.
- Added bounded multi-page CDXML edits and bounded Markush/polymer expansion
  with source-to-expanded mappings across Rust, Python, WASM, and Node.
- Added deterministic occupancy-weighted crystal composition summaries and a
  file-backed SDF/MOL/XYZ benchmark lane.

### Fixed and changed

- Made UFF rescue constraint- and stereo-aware after minimization.
- Synchronized intra-workspace dependency pins used by crates.io publication.
- Removed Spectrophores from the public Rust and Python APIs pending
  independent patent/FTO review; Issue #464's proposed replacement is not
  shipped.

### Performance and release engineering

- Reduced canonical automorphism bookkeeping and SDF property-serialization
  allocations without changing the output contract.
- Made npm publication idempotent and refreshed compatibility, provenance,
  licensing, and benchmark records.

## [1.0.3] - 2026-09-04

- Added reusable MMFF94 topology preparation with cached bonded terms and
  non-bonded exclusions.
- Added reproducible local MMFF94 energy, minimization, ETKDG, and 3D pipeline
  benchmarks. Results remain scoped to the recorded macOS arm64 environment.

## [1.0.2] - 2026-09-04

- Rejected impossible explicit-hydrogen valence states and corrected the
  N-alkylation retrosynthesis template (#455).
- Reduced canonical ranking/search and SDF graph/property I/O overhead.
  Recorded performance claims remain limited to the corpora, operation
  boundaries, hardware, and configurations in [`benchmarks/`](benchmarks/).
- Added a version-pinned, resumable competitive benchmark protocol with corpus
  hashes and explicit unsupported/failure outcomes.
- Normalized the public product name to `chematic` and consolidated local-only
  milestone versions into the v0.89.0 history.

## [1.0.1] - 2026-09-03

- Added up to three deterministic distance-geometry starts to bounded UFF
  stereo rescue. Every accepted candidate must retain finite coordinates,
  sane bond lengths, sound minimization, and declared stereo.
- Kept the v1.0 compatibility and Experimental 3D/MMFF94 boundaries unchanged.

## [1.0.0] - 2026-09-03

- Established the stable bounded API boundary for CDXML/polymer editing,
  partial Python `RWMol`, fail-closed canonical identity, explicit
  aromaticity/CIP modes, and Experimental 3D/MMFF94.
- Added repository-local fuzz, dependency, cross-binding, focused-Miri, and
  Linux sanitizer gates plus three-browser CI evidence.
- Published the GitHub Release, crates.io packages, and PyPI wheels. Later
  releases also established the npm publication path.

## [0.89.0] - 2026-09-01

This release consolidated the local hardening milestones developed between
v0.49.0 and v0.89.0:

- Added finite defaults and typed resource-limit errors across the public
  SMILES, molecule-format, reaction, CLI, Python, WASM, MCP, and 3D boundaries.
- Added shared Rust/Python/Node/WASM adversarial fixtures, fuzz targets,
  focused Miri, sanitizers, dependency/license checks, SBOM/provenance
  generation, immutable workflow pins, and release-key verification.
- Added fail-closed canonical identity, explicit aromaticity/CIP selection,
  bounded standardization/parent reports, and the documented v1.0 scope.
- Kept unsupported, failed, and unmeasured comparison outcomes distinct.

The intermediate labels `0.50.0` through `0.88.0`, and the former local-only
`0.90.0` through `2.32.0` labels, were development milestones rather than
independently published releases. Their changes are represented by v0.89.0 and
the archived detailed history, not by synthetic release entries.

## Earlier public releases

- **v0.24.0-v0.49.0 (2026-08-30 to 2026-08-31):** expanded parent identity,
  comparison contracts, RDKit-oriented APIs, format conversion, and CLI
  workflows.
- **v0.16.0-v0.23.0 (2026-08-15 to 2026-08-30):** added crystal/materials
  formats, Python/WASM format bindings, 3D stereo safety, tautomer work, and
  canonical/standardization correctness fixes.
- **v0.8.0-v0.15.0 (2026-07-29 to 2026-08-14):** developed the bounded 3D
  pipeline, MMFF94 coverage, reaction/stereo handling, and coordination
  chemistry support.
- **v0.1.x-v0.7.0 (2026-05-26 to 2026-07-26):** established the Rust molecular
  graph, SMILES/SMARTS, descriptors, fingerprints, depiction, file formats,
  Python/WASM bindings, and initial validation corpora.

The authoritative list of published tags and release artifacts is the
[GitHub Releases page](https://github.com/kent-tokyo/chematic/releases). Exact
historical implementation notes remain available in the archived detailed
history and Git history.

[Unreleased]: https://github.com/kent-tokyo/chematic/compare/v1.0.7...HEAD
[1.0.7]: https://github.com/kent-tokyo/chematic/compare/v1.0.6...v1.0.7
[1.0.6]: https://github.com/kent-tokyo/chematic/compare/v1.0.5...v1.0.6
[1.0.5]: https://github.com/kent-tokyo/chematic/compare/v1.0.4...v1.0.5
[1.0.4]: https://github.com/kent-tokyo/chematic/compare/v1.0.3...v1.0.4
[1.0.3]: https://github.com/kent-tokyo/chematic/compare/v1.0.2...v1.0.3
[1.0.2]: https://github.com/kent-tokyo/chematic/compare/v1.0.1...v1.0.2
[1.0.1]: https://github.com/kent-tokyo/chematic/compare/v1.0.0...v1.0.1
[1.0.0]: https://github.com/kent-tokyo/chematic/releases/tag/v1.0.0
[0.89.0]: https://github.com/kent-tokyo/chematic/compare/v0.49.0...v0.89.0
