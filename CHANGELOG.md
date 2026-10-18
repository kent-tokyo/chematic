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

- Python `rdInfoTheory.InfoBitRanker` supports class-biased entropy and Pearson
  ranking with typed, bounded, transactional `SetBiasList` inputs. Preferences
  can be replaced without resetting votes. See
  [class-biased ranking](docs/biased-bit-ranking.md) for semantics and scope.

## [1.1.7] - 2026-10-18

### Added

- Rust/Python `rdkit_pick_tautomer` selects from caller-supplied candidates
  using the existing RDKit-port score and canonical-SMILES ties. Empty input,
  oversized input and unsupported candidates return explicit errors; Python
  generators own their yielded molecules and preserve selected metadata.
  See [candidate selection](docs/tautomer-candidate-selection.md).
- Rust/Python `gro_to_lammps_atomic` converts GRO coordinates, periodic cells
  and optional velocities with explicit real/metal units and caller-supplied
  type/mass tables. Orthogonal and restricted triclinic cells retain their
  geometry; exact source and row identity are returned as sidecars.
  See [GRO to LAMMPS units](docs/gro-lammps-units.md) for the atomic-only scope.
- Rust/Python `write_cml_named_checked` exports caller-supplied unique atom
  IDs and optional Cartesian 3D coordinates, including matching bond
  references. Invalid IDs and unsupported graph semantics are explicit
  errors. See [named CML atoms](docs/cml-named-atoms.md) for the MOL2 recipe.
- Rust/Python `sasa_with_radii_checked` accepts explicit coordinates and one
  positive radius per atom, with checked values, work limits and typed errors.
  It returns per-atom and total native Shrake-Rupley areas without generating
  coordinates or choosing default radii. See [checked SASA](docs/checked-sasa-radii.md).
- Rust/Python prepared MMFF94/MMFF94s can use caller-supplied partial charges
  with checked counts, finite values and coefficient range. Explicit charge
  products drive energy, gradient, dielectric updates and minimization.
  See [custom MMFF charges](docs/mmff-partial-charges.md).
- Rust/Python `euclidean_distance_matrix_checked` accepts rectangular finite
  descriptor vectors and returns condensed lower-triangle distances. Bounded
  inputs/work, fallible buffers and scaled norms retain tiny/huge distances.
  See [descriptor distances](docs/descriptor-euclidean-distances.md).
- Added Rust/Python checked sparse binary Tanimoto condensed matrices from
  on-bit indices or Python `GetOnBits()` sources, with bounded extraction,
  resource limits and explicit empty-set policy. See
  [sparse Tanimoto matrices](docs/sparse-tanimoto-matrix.md).
- Added Rust/Python checked single-step reaction CDXML output with explicit
  reactant/agent/product roles, unique object IDs and a caller-supplied arrow
  and coordinates. Unsupported graph annotations fail explicitly; see
  [reaction CDXML output](docs/reaction-cdxml-output.md) for the bounded profile
  and external-reader differences.
- Rust/Python `parse_mopac_archive` reads one MOPAC archive restart deck,
  preserving Cartesian/internal values, references, flags, announced charges
  and raw summary/ancillary sections. See [archive input](docs/mopac-archive-input.md).
- Python compatibility `rdAlignment.GetAlignmentTransform` accepts typed
  point rows, numeric arrays or x/y/z points, and returns weighted SSD plus
  a checked proper/improper transform. See [point alignment](docs/point-alignment.md)
  for input limits and the default-only iteration option.
- Rust/Python `write_gamess_zmat` emits a checked GAMESS-US geometry
  template with source-frame Cartesian `$DATA` and caller-defined `$ZMAT`
  distance/angle/torsion definitions. See [GAMESS output](docs/gamess-zmat-output.md)
  for precision, card limits and calculation-setup scope.
- Python compatibility `rdInfoTheory.InfoBitRanker` accepts typed explicit/
  sparse binary vectors and returns bounded entropy-based feature rankings
  with per-class counts, atomic vote validation and deterministic tie ordering.
  See [entropy bit ranking](docs/info-bit-ranker.md) for the supported profile.
- Python `split_pdb_models` produces validated coordinate-only PDB documents
  with names padded to the total model count. Atom/ANISOU/TER and global CONECT
  lines retain their identity; metadata loss requires explicit consent, and
  aggregate output limits refuse partial results. See [PDB model splitting](docs/pdb-model-splitting.md).
- Python `rdkit_compat.LayeredFingerprint` exposes checked native Layered
  options and accumulates caller-owned list/uint32 NumPy atom counters after
  successful calculation. Strided views, seeded counts and bit masks are
  supported with typed inputs and bounded native enumeration. See
  [Layered caller counters](docs/layered-caller-counts.md).
- Rust/Python `write_xyz_checked` supports caller-selected fixed/scientific
  XYZ coordinates, fractional digits and minimum field widths, with finite
  geometry checks, bounded output and unchanged legacy writer defaults.
  See [XYZ output precision](docs/xyz-output-precision.md).
- Python `plan_sdf_split` creates an exact-count, balanced SDF partition plan
  and bounded chunk iterators. Source bytes and record order are preserved;
  whole-source/part hashes and file identity detect changes without loading
  the complete file. See [SDF file partitions](docs/sdf-file-partitions.md).
- Python `rdkit_compat.rdMIF.CoulombDielectric` adds a typed local-solvation
  model with copied inputs, explicit solvent/solute permittivities and a
  1,024-site construction limit. See [the model and validation](docs/coulomb-dielectric-field.md).

- Python `rdkit_compat.rdMIF.Coulomb` adds typed, copied explicit-charge vacuum
  fields with bounded real inputs and scalar kJ/mol evaluation. See the
  [Coulomb field profile](docs/coulomb-interaction-field.md) for squared-distance
  thresholds, narrowed settings and measured RDKit rounding-boundary differences.
- Python `write_xbs`, `XbsSpecies` and `XbsBondRule` emit checked static XBS
  scene files with explicit gray styles and distance drawing rules. XBS is
  distinct from Open Babel's `bs` dialect. See [XBS output](docs/xbs-output.md).
- Python `parse_xyz_atom_data` / `XyzAtomDataFrame` preserve headerless XYZ
  scalar, vector or flat nine-component payloads with explicit widths, source
  retention and named extxyz conversion. See [XYZ per-atom data](docs/xyz-atom-data.md).
- Python `write_xsf_checked` exports ordered periodic geometry from native
  structures/CIF to XSF, refusing partial occupancy and known unexpanded
  symmetry and requiring explicit metadata-loss consent. See [CIF-to-XSF output](docs/cif-xsf-output.md).
- Rust/Python `FeaturePatternSet` owns up to 32 caller-supplied SMARTS feature
  definitions and produces full-width atom masks or native circular feature
  fingerprints. Matching budgets and size limits refuse partial output;
  empty definitions use no defaults. See [custom feature patterns](docs/custom-feature-patterns.md)
  for the native hashing contract and RDKit comparison scope.
- Rust `SdfInput` / `SdfCompression` and Python `iter_sdf_checked(...,
  compression="plain" | "gzip" | "auto")` support independent per-input
  gzip detection and concatenated members. Plain input remains the default;
  native SDF limits count decoded data and compression failures are reported.
  See [SDF gzip input](docs/sdf-gzip-input.md).
- Python `ForwardSDMolSupplier` reads binary file objects, including BytesIO
  and gzip streams, through the native SDF parser. Both native and compatibility
  APIs reject missing/non-callable `read` at construction and preserve read
  exceptions, titles and SD fields. See [stream input contracts](docs/python-forward-sdf.md).
- Added Rust `Molecule::clear()` and Python `RWMol.Clear()` for complete
  molecule reuse, including metadata/cache reset and Python title/typed
  properties. Existing snapshots stay independent. The compatibility wrapper
  exposes the same chematic extension. See [molecule reset](docs/molecule-clear.md)
  for the Open Babel #3032 reproduction and bounded validation.

### Fixed

- `rdkit_compat.DataStructs.ConvertToNumpyArray` now replaces zero bits as well
  as on-bits, validates source/destination types and dimensions, and reads the
  complete source before modifying an existing array. Bounded inputs and numeric
  NumPy views are covered in [NumPy fingerprint output](docs/numpy-fingerprint-output.md).

- Preserve SDF data fields whose headers contain optional trailing information,
  including RDKit SDWriter's `> <Field> (record-number)` headers, across native
  stream, file and in-memory readers.

## [1.1.6] - 2026-10-16

### Added

- Python file readers, suppliers, writers and HTML report outputs now accept
  `bytes` and `os.PathLike` alongside `str`, using native filesystem encoding
  and preserving custom path-protocol errors. The 15 entry points and local
  platform limits are recorded in [Python filesystem paths](docs/python-filesystem-paths.md).
- Add Rust/Python `write_cif_with_bonds_checked` for explicit molecular CIF
  bond geometry (S/D/T/A), with checked Cartesian distances and no periodic
  bond perception. Related to Open Babel #2726; that upstream report was a
  CLI option mix-up, not missing Open Babel support.
- Add opt-in stereo-checked rooted SMILES in Rust, Python and WASM, with an
  original-atom correspondence, validated tetrahedral frames and a 256-atom
  bound. Addresses the implicit-H bridgehead serialization failure class in
  RDKit #9646 without changing the RDKit-compatibility writer.
- Added opt-in proof-backed removal of one redundant tetrahedral declaration in
  Rust/Python/WASM. An anchored graph witness preserves all other configurations
  and source labels; failed proofs retain the declaration, and unsupported inputs
  or exhausted budgets return errors. Parsing/writing defaults remain unchanged.
- Added opt-in MOL2 Kekulé bond encoding in Rust/Python/WASM. A bounded neutral
  aromatic profile changes only `ar` bond tokens to `1`/`2`, retaining all other
  source bytes and refusing failed matching, valence/H changes and aromaticity
  loss. Existing source round-trips and writer defaults are unchanged.
- Thirteen Python serialized-text readers/validators now accept immutable UTF-8
  `bytes` alongside `str`, retaining native parsing, return shapes and keyword
  names. Invalid UTF-8 raises `UnicodeDecodeError`; unrelated object types raise
  `TypeError` without implicit protocol conversion. Paths and identifiers retain
  their existing contracts.
- Added bounded source-preserving multi-record MOL2 inspection and explicit
  all-reference heavy-atom RMSD in Rust/Python. Atom correspondence and indexed
  bonds are checked; every reference is returned or an error is reported, with
  optional proper rigid alignment and no implicit remapping.
- Added a Rust/Python four-atom alkene stereo query with endpoint-aware
  traversal normalization, explicit cis/trans/unspecified results and invalid
  path errors. Existing RDKit bounds matrices and embedding defaults retain
  their compatibility behavior.
- Added an opt-in MOL2 ATOM → BOND → UNITY section order to Rust and the
  existing Python/WASM graph writers for readers requiring adjacent core
  sections (Open Babel #2350). Default bytes and immutable source retention
  remain unchanged. Open Babel 3.2.1 ignores trailing UNITY charges in six
  original graph controls; see the [compatibility record](docs/mol2-section-order.md).

### Performance

- Skip unused stereo-loss collection for plain V3000 output. Reporting APIs
  preserve all loss fields and serialized bytes remain identical. The fixed
  source comparison measures a 1.5% prepared-write gain on both diverse slices;
  the disclosed parsing retreat and complete scope are retained in the
  [V3000 source record](benchmarks/2026-10-16-openbabel-v3000-loss-work.md).
- Parse complete buffered SDF records directly, reserving owned record storage
  only for the existing line-reader fallback. Preserve diagnostic/lightweight
  results, limits, errors, metadata and logical input positions.
- Append SD data fields directly and skip stereo-loss collection for plain
  V2000 output callers that discard it. Reporting APIs retain every loss field;
  stereo depiction and serialized bytes remain unchanged. See the
  [bounded source comparison](benchmarks/2026-10-16-openbabel-sdf-loss-work.md).
- Fix native XYZ/extxyz readers to recognize uppercase two-letter species
  aliases (CL/BR/CA/FE and the other supported elements), addressing the
  element-identity failure class in Open Babel #2333. Core symbol lookup,
  SMILES and PDB atom-name rules are unchanged.
- Expose opt-in CDXML Bold/Hash/Dash stereo inference to Python and WASM
  (RDKit #9359). Correct screen-coordinate parity, solid/hashed sense and
  directional narrow-end ownership in the native reader. Default
  non-directional inference remains disabled.
- Fix total hydrogen counts to include H atom nodes once and avoid double
  counting bracket H. Correct WASM atom-info attached H and expose the total
  counter there. Found while reproducing Open Babel #3008/#2624; molecular
  formulas and hydrogen transformation algorithms are unchanged.
- Added native explicit PDB CONECT graph reading in Rust, Python and WASM.
  The new API preserves atom/coordinate/serial order and declared bond
  multiplicities without distance-based bonding, bounds input and graph
  resources, and requires opt-in for extended target slots. Existing PDB
  readers keep their behavior. See [the accepted profile and losses](docs/pdb-conect-graph.md).
- Added checked constant-dielectric updates to prepared Rust/Python MMFF94 and
  MMFF94s fields. Energy, gradient and minimization immediately use the new
  screening while retaining prepared atom types, charges and pairs. Invalid
  updates leave the field unchanged; the default remains 1.0. See
  [the API and numerical contract](docs/mmff-dielectric-updates.md).
- Added explicit AUTO/CDX/CDXML byte-block reading in Rust/Python/WASM with a
  bounded concrete-molecule CDX profile, exact source sidecars and reported
  presentation-loss opt-in. Unsupported CDX chemistry and marked stereo refuse.
  The new projection also retains supported geometry-derived E/Z in graph
  output. Existing CDXML entry points are unchanged. See
  [the accepted profile and limits](docs/chemdraw-block.md).
- Added checked rigid PDBQT output in Rust/Python/WASM, rejecting missing or
  excess values, non-finite coordinates/charges, fixed-column overflow after
  rounding, unsafe names and unsupported default atom types. Existing writers
  retain their behavior.
- Added a Python type-stub gate for source and built wheels, covering Python 3.9
  grammar, compiler syntax, typed markers and exact source/package stub bytes.
  CI checks source on Python 3.9/3.13 and source wheels on Linux/macOS/Windows;
  the PyPI workflow checks all wheel stubs before artifact upload and publication.

### Fixed

- Kept degree-two non-aromatic sp centres linear in 2D depictions and improved
  recovery of tetrahedral and alkene stereochemistry from unmarked 3D SDF,
  PDB, and MOL2 inputs under the checked format contracts.
- Corrected aromatic-ligand CIP/E/Z priorities, explicit-hydrogen accounting,
  and dative-carbonyl/ring E/Z rendering without broadening unsupported
  stereochemistry into guessed results.
- Reconciled the integrated Python and WASM binding surfaces and added stable
  MMFF dielectric error-contract tests, restoring the Rust coverage gate after
  the feature integration.

### Verification boundary

- The merged source candidate passed the required hosted CI and the Rust
  coverage gate. Registry publication and exact PyPI/npm/crates.io artifact
  chemistry remain post-publication evidence and are not inferred from source
  tests.

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
