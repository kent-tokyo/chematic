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

- Tinker XYZ (`txyz`) single-frame read/write in Rust, Python, WASM/Node and
  the CLI. Editable records retain opaque atom names, force-field type numbers,
  declared connectivity, Cartesian coordinates and optional periodic cells.
  Invalid topology, extra frames and unrepresentable JSON fields are refused.
- DL-POLY CONFIG/REVCON single-record read/write in Rust, Python, WASM/Node
  and the CLI, retaining source-unit positions, velocities, forces, cell vectors,
  opaque names/suffixes and optional IDs/counts. Unsupported boundary modes,
  inconsistent levels, unknown JSON fields and precision-losing output are refused.
- Molden geometry/vibration record read/write in Rust, Python, WASM/Node and
  the CLI. Explicit atom units, Bohr vibration coordinates/displacements,
  frequency signs, mode/atom IDs and optional IR/Raman data survive edits.
  Electronic, convergence and trajectory sections are refused rather than
  discarded. This is a bounded profile, not a general orbital-file reader.

### Fixed

- XYZ and ASE-style ExtXYZ real reads now share finite binary64 validation
  across Rust, Python and WASM. Nonzero decimal coordinates, cell entries and
  real atom data that round to zero are refused; literal signed zero and finite
  subnormals remain readable. The 3D XYZ path also refuses NaN, infinities and
  overflow. String columns/raw info and existing writer precision are unchanged.
- Python Extended XYZ dictionaries now retain declared R/I/S/L atom-column
  kinds, widths and original comments for typed read/write round trips,
  including empty frames and signed 64-bit identifiers. The writer accepts an
  explicit property schema and refuses lossy integer-to-real coercion;
  duplicate reader columns are refused before dictionary overwrite.
- The CML Cartesian/CJSON projection now resolves namespace prefixes and keeps
  foreign annotations out of the molecular graph. Forward bond references are
  resolved after collecting source-order atoms. Malformed XML, missing element
  or bond-order declarations, and ambiguous/self bond references are refused.
- Checked rigid PDBQT output now refuses isotope annotations that the format
  cannot store (Open Babel #3005). Explicit isotope projection in Rust, Python
  and WASM returns PDBQT plus a per-atom isotope-loss sidecar, retaining every
  explicit atom and caller-supplied charge without ligand torsion preparation.
- V2000 MOL/SDF to ChemicalJSON conversion now retains source x/y/z values
  instead of dropping the coordinate array. Rust and Python expose a bounded
  source-triplet reader, refusing missing/nonfinite fields and retaining flat
  or tiny z values independently of stereo perception. Python also supplies
  this geometry to PDB, XYZ and PDBQT output within each writer's precision.
- MCS stops symmetric search after a valid mapping reaches the absolute atom
  and bond upper bounds, restoring the focused benzene/toluene case from about
  315 µs to 3.93 µs without changing its result. Mixed chiral three-input
  results are input-order invariant; RDKit's order-dependent oracle result is
  recorded rather than copied.
- RDKit-compatible Cleanup serialization now follows transformed molecules
  through Python and pickle, so the public Python and dedicated Rust paths use
  the same writer contract.

### Validation

- Advance Rust, Python/stub, WASM declaration and shared JSON compatibility
  checks to v1.2.0.
- Freeze and score an independent v2 holdout selected after the candidate
  commit. SDF graph/properties/stereo passed 64/64; broad MCS and exact 3D
  coordinate parity remain explicitly incomplete.

## [1.2.1] - 2026-10-23

### Documentation

- Consolidated the active roadmap, open-work ledger, compatibility and error
  contracts; shortened entry documents, the cookbooks, RDKit migration table,
  language-binding guide, MCP reference and GJF integrity contract; repaired
  moved links; retained referenced benchmark evidence and retired 88
  unreferenced intermediate records to the immutable `v1.0.41` history.

### Fixed

- SDF named properties now accept registry information before the field name
  and retain header-shaped lines inside values. Only an empty line ends a data
  item; headers without that separator are no longer guessed as new fields.
- Checked CML → CML conversion now retains the exact XML source, including
  both coordinate sets, properties and unknown extensions, with an explicit
  opaque-semantics report shared by Rust, CLI, Python and WASM.
- CML 2D coordinate output now retains supplied finite binary64 values instead
  of rounding to four decimal places, preserving positions through a rewrite.
- CML readers now decode decimal/hex XML character references and reject
  malformed charge/isotope/2D-coordinate values, conflicting isotope aliases
  and duplicate decoded atom IDs instead of silently changing the graph.
  Shared XML attribute decoding also avoids double-decoding literal entity
  spellings in CDXML/MRV attributes.
- Checked CJSON conversion now preserves standard root molecule names as well
  as the legacy `properties.name` alias. Conflicting names remain in exact
  source records but refuse MOL2 conversion instead of selecting a title.
  MOL2 → CJSON emits and verifies both matching name fields.

- Checked CJSON → MOL2 now retains signed formal charges in explicit UNITY
  attributes, independently of partial charges. Candidate charge verification
  requires exact graph correspondence. UNITY precedes BOND for Open Babel
  3.2.1 readback; isotope/stereo limitations and strict refusal remain explicit.
- Accurate CIP assignment now reports `oracle_unstable` when replacing a
  wildcard with carbon or fluorine reverses a distant tetrahedral R/S label,
  rather than returning a confident label for an unspecified atom.
- Fingerprint similarity search now applies inclusive Tanimoto thresholds
  exactly; a candidate just below the requested threshold is excluded.
- Cleanup no longer adds a spurious positive charge to a metal when a
  disconnected bracket carbon radical remains neutral, including mixed
  organozinc compounds.
- MCS now constructs a common non-induced subgraph, so extra closure bonds in
  one input no longer make the maximum depend on input order (`C1CC1` versus
  `CCC` now returns three atoms and two bonds in either order). Python
  `find_mcs_checked` raises a typed `mcs_resource_limit` error when a search
  limit is reached instead of returning `was_timed_out=False`; WASM reports the
  distinct `wasResourceLimited` state.

## [1.2.0] - 2026-10-21

### Added

- Python `read_xbs` / `XbsScene` read bounded static gray-style XBS scenes,
  retaining exact source, original record order and physical line numbers in
  immutable data. Unsupported commands and unresolved references are errors;
  drawing rules are retained without chemical bond inference.
- Rust `check_geometry` and Python `Mol.check_geometry` inspect supplied 3D
  coordinates without mutation, reporting bond-length outliers, nonbonded
  clashes, skipped references and declared-stereo outcomes independently.
- Typed `rdkit_compat.PatternFingerprint` facade for the native Pattern
  calculator, with exact folding to power-of-two widths through 2048 and checked
  structural mask/list/NumPy counter inputs. Counter and mask arguments remain
  compatibility no-ops, matching the pinned RDKit implementation.

- Checked Python TURBOMOLE ⇄ POSCAR geometry conversion retains 3D lattice
  vectors, unwrapped positions and fully fixed atoms, with source text and an
  explicit atom correspondence map. Unsupported partial constraints and restart
  data are rejected instead of silently discarded.

### Changed

- Route macOS symbol stripping through Apple's linker so source-built Python
  extensions remain loadable under macOS 27's stricter Mach-O `LINKEDIT`
  alignment checks.
- Update the MCP HTTP client's certificate-path validator and bundled root
  store to the releases that remove the distrusted emSign chain.
- Split oversized chemistry, reaction-stereo, standardization, CLI, CJSON,
  depiction, cache-state and evidence helpers into focused modules. The
  refactor is intended to preserve the public 1.x API and is guarded by the
  Rust, Python and WebAssembly compatibility lanes.
- Native RDKit-compatible cleanup uses sparse inline candidate indices and
  borrowed adjacency, preserving correction order and internal molecule state
  while avoiding a molecule-sized flag allocation and redundant scan.
- Python SMILES parsing no longer computes descriptor-model, aromatic-view or
  ring caches eagerly. The first consumer populates the same invalidated
  caches. On the fixed 1,000-row macOS source lane this reduces parse time
  from 5.715 ms to 3.185 ms and improves parse-plus-operation lanes while
  moving preparation cost into first use; prepared reuse is essentially
  unchanged. See the [dated record](benchmarks/2026-10-19-python-lazy-preparation.md).

### Fixed

- `rdkit_smiles_checked` now performs the bounded atom-order checks documented
  for v1.1.9 in addition to its existing one-step fixed-point check. Failures
  expose the first differing atom order and a typed reason. The bound is
  configurable from 2 through 256 and defaults to 32.
- Previous-release compatibility jobs now compare all changed publishable Rust
  crates, Python declarations/dependencies and generated WebAssembly TypeScript
  declarations against v1.1.9. Rust semver checks use the 1.2.0 minor-release
  policy while still rejecting breaking public API changes.
- RDKit-compatible Cleanup now handles three bounded metal representations
  and seven existing-anion acid/base transfers. Exact agreement on the pinned
  stress corpus rises from 3,923 to 3,933/4,069 with zero regressions. Generic
  Reionize remains separate because applying it unconditionally regresses
  neutral Cleanup rows. See the
  [dated record](benchmarks/2026-10-19-rdkit-cleanup-acid-metal-followup.md).
- Added a one-shot scorer and immutable row-level result for the frozen MCS,
  SDF-attribute and stereo/3D holdouts. The score keeps success, mismatch,
  refusal, invalid input and internal error separate and is retained as
  historical evidence rather than tuning data.

## [1.1.9] - 2026-10-19

### Added

- `inchi_keys_checked` preserves every accepted batch input position with an
  immutable success/failure record. Standard IUPAC keys are returned only after
  exact input-identity readback; rejected rows retain their source and reason.

- `fragment_on_bonds_checked` accepts writable one-dimensional NumPy `uint32`
  cut-count destinations, including reversed and strided views. Invalid layout,
  dtype or native fragmentation failures leave the destination unchanged.

- Python `rdkit_smiles_checked` adds a one-step textual fixed-point guard for
  RDKit-compatible canonical SMILES. The v1.1.9 implementation did not yet run
  the bounded atom-order permutations stated in the original release note;
  that contract is completed under `[Unreleased]` above.

### Changed

- Benchmark guidance now includes Open Babel as a separate file-interchange
  comparator, with semantic-equality gates and explicit boundaries between
  historical v1.0.42 evidence, newer one-host source records, and current
  published-package claims. The comparison is now directly available from the
  documentation navigation.
- Benchmark and validation overviews now distinguish the v1.1.8 one-operation
  public-wheel result, the v1.1.9 hosted chemistry workflow, the v1.1.3 durable
  cross-platform packet, current source experiments, and historical results.
- ROADMAP and the open-work ledger now separate strategy from operational
  denominators. The English, Japanese, and Chinese READMEs and current install
  examples distinguish the v1.1.9 release boundary from dated v1.1.8 evidence.
- Added a toolkit-selection guide and clearer use-case entry points for local
  browser chemistry, Rust embedding, Python batch work, MCP, and bounded RDKit
  migration. Removed unpinned throughput claims and examples that discarded
  failed SDF rows or implied experimental property conclusions.

- Morgan hashing, numeric descriptor reads, strict-rotor predicates, Layered
  fingerprints, and CIP assignment avoid redundant preparation while retaining
  their recorded output contracts.

- Native stereo-aware reaction application now rejects alkene templates whose
  product geometry is unspecified instead of silently selecting one outcome.

## [1.1.8] - 2026-10-19

### Added

- Added post-candidate, outcome-blind MCS, SDF-attribute and stereo/3D holdouts
  from hash-pinned NCI and PubChem sources. CI verifies acquisition ordering,
  row hashes and selection policy so the evaluation set cannot drift silently.
- Added previous-release compatibility gates for every changed publishable
  Rust crate, Python declarations and dependency contracts, and generated
  WebAssembly TypeScript declarations.
- Added one source-built RDKit Normalize/Cleanup/Reionize classifier and a
  separate CTfile readback classifier, retaining every residual, refusal and
  oracle failure with the source revision that produced it.
- Python `read_mmcif_chains` returns immutable model/label-chain atom records
  with full residue/atom aliases, source indices and original text. Explicit
  input/work limits and native residue checks protect ambiguous reads. See
  [mmCIF chain records](docs/mmcif-chain-records.md) for the first-loop scope.
- Python `read_pdb_atom_records` exposes immutable single-model ATOM/HETATM
  attributes, including altloc, occupancy and B factor, with source indices
  and original text. Blank optional numbers remain distinct from explicit
  zero. See [PDB atom attributes](docs/pdb-atom-records.md) for limits and scope.
- Python `merge_pdb_coordinate_blocks` assembles validated PDB coordinate blocks
  while preserving ATOM/HETATM identity and numeric fields. It explicitly maps
  chains, remaps serial/ANISOU/CONECT references and retains source/loss sidecars.
  See [PDB coordinate merging](docs/pdb-coordinate-merge.md).
- Python `read_cube_datasets` reads bounded multi-field CUBE inputs, retaining
  source text and dataset IDs. Explicit selection exports one ordinary scalar
  CUBE with an omission sidecar and native geometry validation. See
  [CUBE dataset selection](docs/cube-dataset-selection.md).
- Native MOPAC Cartesian readers/writers accept bounded standalone `+`
  keyword continuations, preserving physical keyword text and the following
  title/comment lines. Signed values and quoted markers do not continue lines.
  See [continuation semantics](docs/mopac-keyword-continuation.md).
- Python `rdkit_compat.GetAdjacencyMatrix` and `GetDistanceMatrix` expose the
  existing bounded native graph profiles as independent NumPy arrays, with
  strict options and refusal of atom-index changes during preparation. See
  [NumPy graph matrices](docs/compat-graph-matrices.md) for supported weights
  and the uncached, 1024-input-atom scope.
- Python `rdkit_compat.rdMetricMatrixCalc` provides typed Euclidean and binary
  Tanimoto condensed matrices as owned float64 NumPy arrays. It reuses checked
  native calculators with staged binary inputs and RDKit's empty-pair policy.
  See [metric matrix facade](docs/metric-matrix-facade.md).
- Python `rdkit_compat.Scoring.CalcEnrichment` adapts bounded ranked records,
  generators and NumPy arrays to the existing native enrichment calculator,
  with explicit binary labels, checked column selection and corrected endpoint
  behavior. See [ranked-record enrichment](docs/scoring-enrichment-facade.md).
- Python `rdInfoTheory.InfoBitRanker` supports class-biased entropy and Pearson
  ranking with typed, bounded, transactional `SetBiasList` inputs. Preferences
  can be replaced without resetting votes. See
  [class-biased ranking](docs/biased-bit-ranking.md) for semantics and scope.
- Python `rdkit_compat.rdInfoTheory.BitCorrMatGenerator` counts selected binary
  bit-pair co-occurrences with typed explicit/sparse inputs, staged state updates,
  bounded work and independent float64 results. See
  [binary bit-pair co-occurrence](docs/bit-cooccurrence.md).
- Python `InfoBitRanker.SetMaskBits` accepts typed integer iterables, validates
  masks before replacing state and limits both future on-votes and ranked
  candidates. Existing counts remain intact; skipped votes are not restored.
  See [feature mask semantics](docs/info-ranker-masks.md).
- Python `rdkit_pick_tautomer_scored` selects supplied candidates with typed
  integer scoring callbacks. Input snapshots and canonical keys are prepared
  before callbacks, scoring views are isolated, and graph/state/count limits
  refuse oversized input. See [scored selection](docs/scored-tautomer-selection.md).

### Changed

- Labute ASA now reuses bounded RDKit-compatible scalar preparation from the
  Python parse boundary and avoids redundant bond/radius work. On the fixed
  10,000-row macOS source lane, values match RDKit and COSMolKit while
  first-use, repeated-use, and parse-inclusive calls win all 20 paired blocks.
  This is not yet public-wheel or two-host evidence; see the
  [dated record](benchmarks/2026-10-19-python-labute-asa-cosmolkit.md).
- Python SMILES parsing now prepares the RDKit-compatible aromatic-view
  decision and ring count for bounded ordinary molecules, so first-use Labute
  ASA and ring-count calls reuse work that RDKit also performs during parsing.
  Larger or highly cyclic inputs remain lazy, graph edits still invalidate
  derived state, and the Rust parser is unchanged. Exact one-host evidence and
  the parse/first-use tradeoff are recorded in the
  [dated benchmark](benchmarks/2026-10-18-python-first-use-ring-asa.md).
- Repeated QED reads on an unchanged molecule now reuse the computed scalar.
  The cache is cleared by every graph edit and clones continue to start with
  empty derived caches; first-use values are unchanged.
- The Python speed-frontier protocol now measures QED and accepts explicit
  comparator module paths while continuing to separate parse-inclusive,
  first-use, repeated-use and peak-process-RSS evidence.

### Fixed

- Release validation now sends only publishable library crates to
  `cargo-semver-checks` and refreshes the current-source standardization
  evidence after SMILES parser changes, keeping patch-release CI fail-closed
  without treating the binary-only CLI as a Rust API surface.
- Native MOL/SDF V2000 and V3000 writers bound titles to 80 UTF-8 bytes,
  keep complete characters and replace controls/line separators with spaces.
  Original metadata stays unchanged; serialized titles can lose content.
  See [title policy](docs/mol-title-width.md) for preservation guidance.

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
- See the [validation report](docs/reference/validation.md) and
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
