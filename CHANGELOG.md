# Changelog

This file summarizes user-visible changes. Exact validation and performance claims belong to dated records tied to the artifact, corpus, comparator, options and host.

Older notes are preserved in the [v1.0.42–v1.1.9 archive](docs/archive/changelog-v1.0.42-through-v1.1.9.md) and [earlier archive](docs/archive/changelog-through-v1.0.41.md). Tagged releases and [GitHub Releases](https://github.com/kent-tokyo/chematic/releases) retain publication history.

The format follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/), and public releases follow [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Added

- `stereo_depiction::with_stereo_depiction_on` draws a molecule's stereo on
  caller-supplied coordinates (for example a file's own drawing): one wedge or
  hash per tetrahedral centre, chosen and checked on those coordinates, as the
  MOL writers do. A centre whose wedge was on a hydrogen that
  `remove_hydrogens` dropped is drawn with a wedge again.

### Fixed

- The CDXML reader finds atoms and bonds whose element name is followed by a
  line break (`<n` then `id="..."` on the next line), as ChemDraw writes them;
  such files used to parse as an empty molecule without an error.

## [1.2.5] - 2026-11-26

### Fixed

- Reaction SVGs position components by their bounding boxes, with consistent
  spacing and theme-aware arrows and plus signs (#1294).
- Neutral three-coordinate boron uses trigonal-planar geometry in rule-based
  placement, ETKDG constraints and the legacy minimizer (#1296).
- Selected main-group centers with five or six electron domains use VSEPR
  placement and suitable angle handling, including ring cases (#1297).

## [1.2.4] - 2026-11-23

### Added

- Added opt-in checked Python conformer ensembles, including a mixed
  compatibility/corrected torsion-gradient profile, declared-stereo checks
  and bounded near-linear nitrogen repair. See the
  [checked ensemble guide](docs/reference/checked-conformer-ensembles.md).

### Fixed

- Corrected aromatic substituent placement, ETKDG ring-bond handling and
  DREIDING X–H bond parameters for either atom order (#1287).
- Corrected hybridization of two-coordinate cationic hydrogen for UFF
  bridging-hydrogen cases (#1282).
- Isolated generation and scoring errors by engine in the timed conformer
  evaluator. Previously unscored outputs require new measurements (#1282).

### Validation

- Recorded Linux and macOS source-wheel development and holdout evaluations
  for conformer recovery and common-coordinate force-field comparisons.
  The [evaluation index](benchmarks/README.md) links each protocol and
  its limitations. These measurements are tied to their recorded candidate
  wheels and hosts.


## [1.2.3] - 2026-11-22

### Compatibility

- Rust `cleanup(&Molecule)` retains its best-effort return type. The new
  `cleanup_checked` reports metal transformation errors. Python and WASM
  use the checked path. See [Cleanup migration](docs/reference/cleanup-migration.md).

### Fixed

- Cleanup preserves metal and ligand hydrogen counts in supported cases, retains selected aromatic copper coordination across input spellings, and preserves explicit zero hydrogen on a disconnected imine. It reports unresolved metal hydrogen changes, unstable detached metals and unsupported five-ligand tetrahedral centers.
- Bonded atom labels render charge, hydrogen count and isotope consistently in the shared SVG depiction path (#1281). EPS renders formal charges on atom labels (#1283).
- Corrected diagnostic format strings and widened traversal arithmetic in the bundled native InChI backend (#1273–#1277, #1279).

### Release checks

- Accept either order of compatible manylinux wheel tags and run the published Python 3.9 chemistry gate on Ubuntu 24.04.
- Keep the [v1.2.2 artifact verification](validation/results/release-channel-verification-v1.2.2.json) and [published comparison](benchmarks/2026-11-21-v122-published-comparison.md) as the baseline for release review.

## [1.2.2] - 2026-11-21

### Added

- Added bounded, editable readers and writers for DL-POLY HISTORY and CONFIG/REVCON, NWChem GEOMETRY, Tinker XYZ, Molden geometry/vibrations, and nonperiodic BIOGRF 200 across documented Rust, Python, WASM/Node and CLI surfaces. Unsupported variants and lossy edits return explicit errors.
- Extended opt-in RDKit-compatible fingerprints for Morgan, atom-pair, torsion, path, Layered and MACCS/count-simulation workflows.
- Added a bond-first MCS objective and bounded SMARTS output for observed atom/bond alternatives.

### Fixed

- Strengthened CIP and stereo handling for isotope priorities, auxiliary rules, lone-pair centers, cumulenes and allene ligand order. Unresolved configurations continue to abstain.
- Improved focused aromatic MCS search while retaining complete-ring, multi-input stereo and resource-limit checks.
- Kept RDKit-compatible Cleanup's transformed writer contract through Python and pickle; corrected bounded acid/base, charge, aromaticity and explicit-haptic-iron cases.
- Preserved typed source data in XYZ/ExtXYZ, V2000 MOL/SDF, CML, MolJSON and checked PDBQT conversion. Invalid or lossy values are refused or reported explicitly.

### Validation

- Moved Rust, Python/stub, WASM declaration and shared JSON compatibility checks to the v1.2.1 baseline.
- A pre-publication macOS arm64 wheel gate with RDKit 2026.03.6 recorded one CIP difference in the exposed 10,000 corpus: a native-parsed bridgehead amine is R where RDKit reports S; the RDKit-parsed chematic model reports S. Published-wheel gates retain this in the denominator.
- The v1.2.2 source Cleanup gate against RDKit 2026.09.1 matches 4,022/4,069 comparable stress rows and 9,996/10,000 broad rows. Relative to the 2026-10-24 source record, 47 and 2 previously exact rows differ. The record names every row and keeps the older result separate. Reionize matches 4,072/4,072 stress rows. These are fixed-corpus string comparisons, not chemical truth or published-wheel results.
- Earlier frozen holdout evidence remains historical: SDF graph/properties/stereo matched on 64/64 rows; broad MCS and exact 3D coordinate parity remained incomplete. It does not validate all v1.2.2 changes.

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
