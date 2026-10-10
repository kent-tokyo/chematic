# Changelog

This file keeps concise, user-visible release summaries. Full release details
through v1.0.41 are preserved in the
[changelog archive](docs/archive/changelog-through-v1.0.41.md); validation and
performance claims remain scoped to their dated records.

The format follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and public releases follow [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Added

- Added bounded, label-stable R-group decomposition for Rust and Python.
  Terminal mapped wildcards such as `[*:1]` and mapped core atoms now produce
  stable `R1`/`R2` rows and RDKit-style column dictionaries. Input atom order
  is normalized before matching, output wildcards retain their labels, and
  duplicate labels, non-terminal placeholders, ambiguous multiple
  attachments, and match-limit exhaustion fail explicitly. Multi-core/MCS
  alignment, tautomer/enumeration expansion, and RDKit's GA scoring remain out
  of scope for this bounded API.
- Added a versioned Open Babel 3.2.1 file-I/O comparison contract for eight
  production formats, a static CI checker, and a unified record-accounting
  runner. The initial 20-repetition record proves fixture provenance and
  record/failure accounting only; semantic round-trip and equivalent-work
  speed claims remain open.
- Added a shared-observer semantic round-trip gate and a paired 21-block CLI
  round-trip benchmark for V3000, MOL2, CML, and CDXML. Evidence remains scoped
  to the checked-in fixtures and source candidate. A same-process harness now
  separates parse, write, and round-trip hot loops; all 12 bounded lanes pass
  the paired speed gate against Open Babel 3.2.1. Broad corpora, large-file
  throughput, memory, and published artifacts remain open.
- Added loss-aware `parse_mol2_record` and `write_mol2_record` Rust APIs. They
  retain Tripos atom types, partial charges, residue data, status bits,
  `UNITY_ATOM_ATTR` formal charges, and opaque extension sections.
- Exposed the v1.0.40+ RDKit interoperability surface to WASM/Node: SMARTS
  and structure writers, PDB/XYZ/MOL2 readers, Murcko/stereo/hash/Morgan
  helpers, bounded alignment/RMSD, and separately named seeded ETKDG and
  bounds-matrix operations. JSON atom indices and coordinate rows use the
  zero-based `MolHandle` atom order; malformed, unsupported, and oversized
  inputs throw stable JS errors instead of panicking (#784).
- Added checked Rust and Python entry points for Balaban J, IPC, and the legacy
  MMFF94 BCI charge model. They return typed errors when the calculation is
  outside its documented size or element table instead of making a valid
  `0.0` result indistinguishable from an unavailable calculation. WASM also
  exposes a typed JSON result for checked legacy MMFF94 charges.
- Added `num_rings()` to the WASM molecule handle and regression fixtures for
  coordination and disconnected-component ring-closure notation.
- Added opt-in RDKit-profile reaction rejection diagnostics across Rust,
  Python, and WASM. Rejected product sets now report stable valence,
  aromaticity, kekulization, or unknown reasons plus bounded atom indices,
  template maps, elements, explicit H counts, and reliable valence limits.

### Changed

- Reduced temporary allocation in loss-aware MOL2 serialization while
  retaining the common-observer semantic result.
- Kept rejected RDKit-profile reaction products in a new opt-in detailed
  report. The existing Rust `TracedReactionTransformReport` again has its
  original two-field struct-literal shape, avoiding a source-compatibility
  break for downstream Rust users.
- Updated published-wheel chemistry and A6 workflow defaults to v1.0.42 and
  added the hash-pinned v1.0.42 chemistry baseline. Historical v1.0.38 results
  remain available under their own expected-value file.
- Rebased the active performance gate to RDKit 2026.09.1 and COSMolKit
  0.5.0rc22. The benchmark adapter follows the rc22 Morgan and SMILES APIs,
  treats any per-row engine error as an invalid timing, and supports one
  interpreter override for all three engines without creating a spurious
  fourth engine.
- Cached the shared RDKit-model disagreement predicate during validated SMILES
  parsing, skipped a duplicate cleanup scan for ordinary molecules, and
  removed Labute ASA bond-order and per-atom heap allocations. Descriptor,
  fingerprint, and ring-count outputs remain unchanged on the pinned ChEMBL
  and NCI agreement corpora.
- Avoided per-bond closure-map probes when RDKit-order traversal has no SMILES
  ring closures. The current arm64 source candidate keeps TPSA and chiral
  Morgan faster than RDKit 2026.09.1 in all 21 measured blocks; direct Labute
  ASA and first-use ring count remain explicit deficits.

### Fixed

- Matched RDKit 2026.09.1's RingDecomposerLib ring counts for large symmetric
  macrocycles while retaining the legacy-equivalent fast path for compact
  ring systems. The pinned 10k corpus improves from 9,994 to 10,000 exact ring
  counts; this is source-candidate evidence, not yet a published-package claim.
- Stopped interpreting MOL2 partial charges as rounded formal charges and
  retained `.ar` atom aromaticity independently of aromatic bond rows.

- Fixed RDKit-compatible ring counts for SMILES closure notation that creates
  a coordination or disconnected-component bond without a graph cycle. Four
  v1.0.42 regressions now return zero rings instead of one.
- Preserved caller-supplied explicit hydrogen atoms in the RDKit 2026.03.6
  reaction profile. This fixes a deletion/edit case that emitted a carbon
  radical after an `AddHs`-style input; the fuzz report now records overlapping
  deletion, bond-break, radical, and re-sanitization facets separately.
## [1.0.42] - 2026-10-10

### Changed

- Expanded Rust coverage measurement to exercise native InChI reconstruction,
  canonical SMILES, fixed-column PDB parsing, reaction SMARTS grammar,
  stereochemical MOL readers, reaction caches, MMFF94 terms, and cross-binding
  contracts. The report now excludes the ignored manual 3D differential
  harness from the production-code denominator.
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
