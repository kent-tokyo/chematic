# Changelog

This file records public releases and the current unreleased changes to
`chematic`. Historical roadmap and audit notes are retained in
[`docs/archive/README.md`](docs/archive/README.md).

The format follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and public releases follow [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

## [1.0.23] - 2026-09-24

- Raised RDKit 2026.03.6 output agreement for RDKit-defined operations
  (row-level, versioned record in `benchmarks/`): RDKit-compatible atom-pair
  and torsion fingerprints (hypervalent S/P/Se pi-electron gate and
  RDKit-perceived aromatic input), MACCS (verbatim key 26, RDKit ring info for
  key 125, key 1 unset like RDKit), QED (every property per
  `rdkit.Chem.QED.properties`), Bemis–Murcko scaffolds (linkers of any length,
  exocyclic double bonds, preserved stereo), and Kekulé-input Crippen hydrogen
  and TPSA nitrogen typing. ChEMBL 5k: atom pair 4,249 -> 5,000, MACCS
  4,714 -> 4,995, QED 3,201 -> 4,962, Murcko 1,591 -> 4,998.
- SMARTS follows Daylight/RDKit semantics more closely: an unspecified bond
  matches single or aromatic bonds (not any bond), full bond-expression
  precedence (`=@`, `=;@`, `-&@`, `=@,-`), `.`-separated components, bare `r`
  and multi-digit `rN`. Queries that relied on an implicit bond matching
  double/triple bonds must now write `~` or the explicit bond.
- Added `rdkit_tpsa` (Rust and Python): RDKit-default TPSA (N and O only);
  `tpsa` keeps including S and P.
- SSSR candidate selection deduplicates with a bond bitmask and computes
  tie-break keys per ring length; the selected rings are unchanged.

## [1.0.22] - 2026-09-24

- Optimized hot paths without changing the checked source outputs: molecule-local
  derived caches, lazy SSSR materialization, cached RDKit-parity aromaticity and
  canonical ranks, a smaller SMARTS/VF2 search state, and reusable MACCS/Crippen
  preparation. The accompanying differential digest record covers 1,074,928
  operation rows across 46,736 molecules with zero main-to-optimized differences.
- Added a versioned Python operation matrix against RDKit 2026.03.6, including
  raw repetitions, corpus/wheel hashes, and streamed line-by-line output digests.
  Of 43 output-comparable operations, 21 have exact output agreement and a faster
  median; 18 are faster in every recorded repetition. This is a source-level,
  shared 2-vCPU x86 record, not an Apple-silicon, WASM, or published-package claim.
- Corrected the earlier aggregate speed wording: non-equivalent operations and
  the non-standard `Mol.inchi` representation are not counted as RDKit wins.

## [1.0.21] - 2026-09-23

- Fixed canonical SMILES for aromatic-stash E/Z systems so two directional
  carriers at one alkene end cannot encode the same side. The RDKit 2026.03.6
  exposed 10,000-row semantic lane improves from 18 stereo-identity losses to
  zero, while graph identity, CIP, Morgan, and SMARTS counts remain unchanged.
- Condensed the active roadmap, Trust Release plan, performance plan, and
  open-work ledger around current priorities and acceptance gates; refreshed
  the English, Japanese, and Chinese READMEs plus validation and benchmark
  indexes without relabeling historical measurements as current results.
- Added pinned MkDocs build dependencies, a contributor-facing strict-build
  command, and a strict Pages build. Repository artifacts referenced from the
  documentation now use explicit GitHub links so the published site validates
  without suppressing link warnings.

## [1.0.20] - 2026-09-23

- Completed the process-level Criterion gate calibration from issue #70. A
  hosted +10% regression was blocked, and a +5% regression was detected in
  8/10 independent hosted runs without weakening the predeclared sign test.
- Made the independent-build null control actionable at the configured 95%
  Wilson/0.2 boundary by increasing it from seven to ten blocks. A 1.04
  practical-effect floor prevents small build/codegen differences from
  contaminating a run, while a hosted one-sided-bias calibration now downgrades
  an otherwise clear Stage 2 failure to environment-inconclusive.
- Kept tetrahedral chirality and 2D wedge/hash depiction synchronized when
  `invert_stereocenter` is called, including the existing WASM-facing API.
- Canonical-SMILES validation now includes the exact RENKIN issue #128 target
  behind issue #372. The retained exhaustive and current orbit-pruned engines
  produce identical output; the exact witness drops from 24 search leaves to
  1 and measured 4.86x faster by the five-run paired median on the recorded
  local release build.
- Closed the measured canonical E/Z shared-carrier residual from issue #149.
  A clean-main audit now passes all 28 coupled components across 1,024 seeded
  atom relabelings each with zero divergent outputs or correspondence failures;
  all four historical residuals are permanent equivalent-spelling regressions.
  Unmeasured coupled shapes keep the conservative stable-key refusal boundary.

- Preserved real carbon, generic `*`, unnumbered `R`, and numbered
  `R1..R9999` atom identities across MOL V2000, MOL V3000, SDF, and CDXML
  round trips instead of serializing wildcard-backed atoms as carbon.
- Added a private molecule sidecar plus stable Rust/WASM inspection and edit
  APIs for R-group labels. Unsupported CDXML labels and out-of-range R-group
  numbers fail closed; real Ru, Rh, Re, and Rn atoms remain elements.
- Completed the RDKit 2025.09.3 A0 core-eight descriptor gate: a frozen
  candidate passed 2,000/2,000 development rows and a separately sourced,
  overlap-audited one-time sealed holdout at 8,000/8,000 for molecular weight,
  HBA, HBD, TPSA, LogP, molar refractivity, Fsp3, and aromatic-ring count.
- Replaced broad descriptor aromaticity promotion with the bounded RDKit-parity
  model, preserved authoritative zero-aromatic results, and removed a
  shape-based cyclic-ether Crippen shortcut that misclassified one fused system.
- Fixed mixed aromatic/Kekulé perception so RDKit-parity descriptors promote
  morphine/codeine's aromatic bridge oxygen without changing the established
  aromatic partition of large fused cages or unrelated cyclic ethers.
- Added post-freeze acquisition attestation, strict eight-field evaluation,
  commit-safe provenance summaries, and a fail-closed checked-in A0 validator.
- Resolved MMFF94 issue #337's six pyridinium/macrocycle residuals by retaining
  the complete symmetrized ring family, restoring RDKit-compatible MMFF ring
  processing order, and treating `/`/`\` markers as directional single bonds.
  The 265-molecule gate now has zero aromatic atom/bond mismatches; MMFF types
  and charges are exact on every comparable row, with only the separately
  declared unsupported probe remaining outside parity.
- Fixed Python similarity-search result indices so filtered inputs retain their
  original positions, and aligned scalar and indexed Tanimoto handling for
  unequal fingerprint lengths.
- Added reproducible RDKit 2026.03.6 Python/npm comparison packets with exact
  artifact provenance and complete 10,000-row accounting. Residual SMILES,
  CIP, SMARTS, and unsupported Morgan cases remain explicitly classified.
- Added immutable prepared RDKit-compatible Morgan fingerprints and retained a
  separate native ECFP profile, avoiding repeated preparation without changing
  configured-bit output on the supported 9,999-row comparison domain.
- Promoted the gated analytic MMFF94 gradient to the production 3D pipeline,
  increased L-BFGS history and iteration budget, and retained typed convergence
  diagnostics.
- Added a production stereo-safe MMFF94 benchmark arm and constrained line
  search for the bounded expanded-H/heavy-only stereo mismatch. The fixed
  265-molecule source gate now has 265/265 independently sound, stereo-clean
  outputs with zero gross clashes, closing the prior 12 typed failures and four
  clash rows. This quality lane is slower than RDKit in the recorded source
  packet; no quality-equivalent MMFF94 speed win is claimed.

## [1.0.19] - 2026-09-22

- Replaced ad hoc TPSA functional-group exceptions with the complete RDKit
  2025.09.3 Ertl N/O/S/P environment table and one shared total/per-atom
  contribution path.
- Added rare atom-type regressions covering 63/63 public cases and
  10,000/10,000 rows across two exposed public corpora.
- Froze the replacement candidate before source acquisition, audited a new
  14,764-row ChEMBL source against 20,000 exposed rows, and recorded TPSA
  agreement on 2,000/2,000 development rows plus a one-time 8,000/8,000 sealed
  holdout at `1e-6` tolerance with maximum absolute error `0.0`.
- Made sealed-cohort acquisition resumable with bounded retries and response
  caching, and added a provenance-bound frozen-candidate evaluator.
- This closes the known TPSA residual only. The complete multi-field A0 gate
  remains open, and the exposed sealed rows cannot validate later candidates.

## [1.0.18] - 2026-09-21

- Added a shared Rust/Python/WASM/Node contract for unknown-length SMILES
  canonicalization streams. Normal EOF processes every observed row; an
  interrupted source reports the processed prefix, observed-but-unprocessed
  rows, and an explicitly unknown unread suffix.
- Added typed stream terminal reasons (`cancelled`, `time_limit`,
  `resource_limit`, `producer_error`, and `consumer_closed`). Unknown reasons
  are rejected by Python/WASM bindings instead of causing a panic or inventing
  successful/skipped rows.
- Extended the shared binding fixture, manifest validation, and generated
  compatibility dashboard to cover the 58-operation contract inventory.
- Documented the remaining T3.7 work explicitly: retry/export semantics,
  Worker/MCP adapters, clean-install coverage, and measured 10k resource
  limits remain open.

## [1.0.17] - 2026-09-20

- Added bounded SMILES+ extended ring-closure support: `%nn` remains the
  legacy two-digit form and `%(n)` is accepted and emitted for labels of 100
  or greater. Invalid, overflow, and ambiguous legacy `%100` forms reject
  rather than changing the molecular graph.
- Preserved the existing `u8` SMILES error variants for legacy labels and
  added distinct extended-label errors only above that range, avoiding a
  breaking public error-field change.
- Added a 102-cycle canonical-SMILES round-trip regression and a 121-cycle
  non-canonical writer/parser regression.
- Added explicit CXSMILES `_AP<n>` parsing and degree-one wildcard attachment
  identity inspection. This does not claim MDL attachment-point collapse or
  V3000 `ENDPTS`/`ATTACH` semantic support.

## [1.0.16] - 2026-09-16

- Updated the transitive `rustls` dependency to 0.23.45, addressing
  RUSTSEC-2026-0285. The patched lockfile has passed local `cargo audit`; the
  hosted Security Audit rerun remains required before a release claim.
- Preserved the typed `UnsupportedCoordinationSanitization` outcome in both
  RDKit Morgan/ECFP4 evidence-dump examples instead of failing their complete
  corpus runs on that documented unsupported boundary.
- Made native `chematic-wasm` V3000 tests exercise their fallible internal
  helpers rather than aborting through `wasm-bindgen` host conversions, and
  corrected the fixtures' SGROUP counts.
- Added a published-npm, 10,000-row, 20-repetition Chromium cold local
  download-to-ready record against official RDKit.js. It records the complete
  JS/WASM asset path separately from workload timing and does not claim
  internet or cross-host startup performance.
- Added a dispatchable Linux/macOS/Windows PyPI-wheel runtime smoke that saves
  version and basic SMILES-operation evidence as per-platform artifacts.
- Added a separate public-npm Chromium process-tree RSS record, explicitly
  preserving its shared-page double-counting limitation.
- Added version-pinned RDKit and Indigo V3000 semantic round-trip gates,
  including SGROUP ordering/edit preservation and explicit unsupported
  coordination/ENDPTS boundaries.
- Added a 300-structure Stereo Torture Suite and a 5,115-input SMILES-spelling
  invariance gate for the development corpus.
- Froze an annotated Trust evaluation candidate and created a post-freeze,
  attested 2,000-row development / 8,000-row sealed split. No sealed-holdout
  accuracy result has been calculated yet.
- Updated maintenance dependencies, including `ureq` 3.4.1, `jsonschema`
  0.56.0, and provenance action 4.2.2.
- Consolidated duplicated Node/WASM cross-binding dump scripts into one
  mode-driven entry point and simplified the active roadmap and validation
  documents.

## [1.0.15] - 2026-09-13

- Added a release-candidate safety boundary: prerelease tags no longer create
  GitHub Releases or publish crates.io, PyPI, or npm artifacts.
- Added a Linux-enforced parser-security corpus gate with bounded execution,
  memory accounting, and isolated network namespace coverage for malformed
  SMILES, SMARTS, V2000, V3000, and SDF inputs.
- Moved Local Compound Explorer parsing and descriptor analysis into a Worker,
  with a 10,000-record input cap, bounded rendering, deterministic cancellation
  feedback, and native scalar/batch/Worker parity checks.
- Added cross-browser Explorer validation for Chromium, Firefox, and WebKit,
  including the 10,000-record Worker workflow and a fresh-runner native CLI
  parity build.
- Added provenance-sealed cohort preparation and attestation checks for the
  Trust Release evaluation workflow; preparation evidence is kept distinct
  from a completed external-oracle measurement.

## [1.0.14] - 2026-09-13

- Reduced false `RingModelAmbiguous` refusals in opt-in RDKit-parity SMARTS
  matching by keeping `[R0]` ring-membership negation outside the positive
  ring-count model; the 5,021-molecule × 31-pattern gate now records 18
  explicit refusals and 155,618/155,651 parity matches (99.9788%), with zero
  regressions or worsened cells.
- Added a diagnostic potential-stereocenter parity gate against RDKit's
  `CalcNumAtomStereoCenters`; it compares 5,000/5,000 rows and records five
  explicit residuals after narrow aromatic-sulfur and deeper-branch fixes,
  without promoting the profile to compatibility status.
- Fixed a panic in `depict_data_with_coords` when callers provide fewer
  coordinates than atoms; missing positions now follow the documented
  `(0, 0)` fallback and have regression coverage.
- Stopped silently replacing MMFF94 charge-calculation failures with zero
  charges. Energy, torsion-scan, and minimization APIs now return a typed
  `MinimizerError::ChargeCalculation` instead.
- Refactored Python bulk descriptor output into a dedicated row type and NumPy
  materialization helper without changing the public result schema.
- Centralized the `pipeline_v2` ring-torsion applicability predicate so
  diagnostic evidence and fail-closed policy checks cannot diverge.
- Removed redundant descriptor-selection conditions and revalidated the
  affected Rust, Python, depiction, and 3D paths.

## [1.0.13] - 2026-09-11

- Rejected multiple top-level elements in the opt-in strict CML parser while
  preserving the historical lenient parser behavior.
- Added the multiple-root rejection to the shared Rust/Python/Node/WASM CML
  contract so all supported bindings exercise the same fail-closed boundary.
- Added typed CDXML text/caption style access for font, size, and alignment,
  with non-finite numeric values rejected and unknown presentation attributes
  preserved.
- Added typed CDXML page `BoundingBox` access with finite-number and arity
  validation.
- Added typed CDXML object classification for atoms, bonds, groups, arrows,
  captions, text, and drawing objects, with an explicit unsupported variant.
- Fixed loss-preserving CDXML parsing for valid self-closing empty pages.
- Included typed CDXML object kinds in the document JSON summary for binding
  consumers.
- Added provenance-backed `contracted_smiles` to semantic expansion JSON so
  Python and WASM callers can restore the exact expansion base graph.
- Added typed semantic S-group and polymer-linkage definitions with stable IDs,
  source-atom validation, repeat-unit references, and explicit unsupported-kind
  handling.
- Added stable-ID S-group kind editing and explicit rejection of unimplemented
  biomolecule payloads.
- Added S-group and polymer-linkage source-to-expanded provenance mappings to
  semantic expansion results.
- Added stable-ID editing for R-group allowed substituent sets, with selection
  reset and validation of the replacement set.
- Added the fallible `rdkit_ecfp4` profile to `PreparedFingerprintIndex`,
  including original-index failure reporting and deterministic top-k ties.
- Added the three-lane native/native, RDKit-compatible/RDKit, and cross-profile
  top-10 comparison gate; the v1.0.12 valid-input run reached 100% mean recall
  in the compatible lane and reports preprocessing failures separately.
- Added versioned WASM/JS document-binding APIs for rich reaction and CDXML
  JSON, bounded stable-ID/page edits, exact CDXML source reserialization, and
  structured error categories with paths.
- Added a deterministic 155-row CIP oracle-corpus atom-order gate with eight
  generated permutations per row (1,240 checks); resolved labels and explicit
  fail-closed skip reasons must remain invariant after target remapping.
- Added a bounded typed V3000 SGROUP syntax view for ID, parent, kind, atom
  references, and attributes while retaining the original logical line for
  lossless writing; polymer expansion and cross-engine semantics remain
  explicitly outside this contract.
- Exposed the typed V3000 SGROUP syntax view through a bounded WASM JSON API,
  including explicit handling for unknown kind tokens and output-size limits.
- Added fail-closed validation for V3000 SGROUP duplicate IDs, missing parent
  groups, and atom references outside the parsed molecule.

## [1.0.12] - 2026-09-11

- Added bounded reaction requirements and stoichiometry analysis for
  `PreparedReaction`, including conservative prefiltering and explicit
  diagnostics for unsupported or oversized queries (#510, #513).
- Added reaction application support for the documented atomic-number query
  subset across Rust, Python, and Node/WASM bindings, with fail-closed errors
  for unsupported SMARTS and regression coverage for product enumeration.
- Preserved opaque V3000 `SGROUP` metadata through MOL parsing and writing,
  including the WASM round-trip boundary, and added generated V3000 gate
  coverage for binding behavior.
- Added shared semantic expansion fixtures and APIs for Markush selections
  and polymer repeat commands, with exact source-to-expanded atom mappings and
  explicit limits for bounded expansion.
- Added CIP regression coverage for negative-charge resonance systems and
  atom-order shuffling, ensuring that resolved labels remain invariant to
  input atom order and malformed internal state fails closed.
- Added reproducible browser/WASM comparison gates against RDKit.js, covering
  artifact size, initialization, SMILES parsing/writing, ECFP4, and memory
  accounting; recorded the generated and official RDKit.js comparison
  artifacts separately from native parity results.
- Updated the browser compatibility workflow and release evidence manifests
  so WASM validation, benchmark provenance, and publishable package versions
  are checked together.
- Synchronized current README, roadmap, benchmark, compatibility, and binding
  documentation with the implementation and corrected stale WASM size/date
  references while preserving historical records and their pinned versions.

## [1.0.11] - 2026-09-10

- Split V3000 MOL atom and bond record parsing from the state-machine driver,
  preserving the existing format, stereo, coordinate, and error contracts.
- Split pipeline configuration, torsion configuration, and authoritative final
  stereo verification helpers from the 3D embedding orchestration path.
- Hardened the experimental CIP resolver so malformed internal rank/position
  state fails closed instead of panicking, with a regression test.
- Added a repository-wide code audit record and revalidated the standard local
  release gates.
- Published the v1.0.11 Rust crates to crates.io and the WASM package to npm;
  the PyPI publication workflow remains in progress.
- Added bounded #510 reaction-template application for supported atomic-number
  primitives, with explicit unsupported-query failures and cross-binding tests.
- Added the #513 `ReactionRequirements` API derived from `PreparedReaction`,
  including conservative lower bounds and a fail-open `could_match()` prefilter.
- Added opaque V3000 `SGROUP` metadata round-tripping and the explicit WASM
  `roundtrip_mol_v3000_block` boundary.

## [1.0.10] - 2026-09-09

- Added a versioned `MoleculeExtension`/`ExtensionRegistry` contract in
  `chematic-core` and documented contributor corpus provenance, licensing,
  minimization, and oracle-version requirements.
- Added a reproducible 4,500-library/500-query exact similarity-search
  comparison runner and v1.0.9 RDKit benchmark artifact.
- Optimized prepared similarity search by caching database popcounts and using
  partial top-k selection while preserving the existing result contract.

- Added a bounded reaction-application compatibility path for simple SMARTS
  atomic-number primitives (\`[#N]\` and \`[#N:map]\`). It deterministically
  preserves map text, distinguishes valid aliphatic/aromatic alternatives,
  applies them through reusable \`PreparedReaction\` and ring-aware paths, and
  exposes additive per-variant diagnostics, and fails closed for compound
  primitives or expansion blow-ups.
- The per-variant diagnostics also cover caller-provided ring perception,
  without changing product enumeration.
- Added a 5,000-row native ECFP4 byte-parity report across Rust, Python, and
  Node/WASM, with the exact 2,048-bit representation and corpus hash recorded.
- Added the corresponding native MACCS 166-bit binding-parity report across
  Rust, Python, and Node/WASM, while keeping RDKit key mapping separate.
- Added a 5,000-row native `topo_path` binding-parity report across Rust,
  Python, and Node/WASM, separate from RDKit-compatible path accuracy.
- Added a 5,000-row native torsion binding-parity report across Rust, Python,
  and Node/WASM, separate from RDKit-compatible torsion accuracy.
- Added a 5,000-row explicit standardization-profile parity report across
  Rust, Python, and Node/WASM, separate from expected-identity accuracy.
- Added a 5,000-row RDKit-compatible ECFP4 binding-parity report across Rust,
  Python, and Node/WASM, including typed preprocessing-failure accounting.
- Added a 5,000-row RDKit-compatible hashed torsion binding-parity report
  across Rust, Python, and Node/WASM.
- Added a 5,000-row RDKit-compatible ECFP4 sparse identifier/count parity
  report across Rust, Python, and Node/WASM.
- Added a 5,000-row RDKit-compatible ECFP4 folded bitInfo provenance parity
  report across Rust, Python, and Node/WASM.
- Added a 5,000-row RDKit-compatible ECFP4 raw bitInfo provenance parity
  report across Rust, Python, and Node/WASM.
- Added a 5,000-row RDKit-compatible path fingerprint parity report and a
  corresponding WASM `rdkit_path_bitvec` API, separate from native `topo_path`.
- Added a reproducible 5,000-row Rust/Python/Node-WASM descriptor binding
  parity report. All bindings parsed the same corpus and agreed on MW, TPSA,
  HBD, HBA, and heavy-atom count; the report explicitly remains separate from
  the external RDKit accuracy ledger.
- Added a shared semantic-expansion fixture contract for Markush selection and
  polymer repeat commands, consumed by Rust, Python, and Node/WASM tests with
  exact source-to-expanded atom mappings.
- Expanded the shared standardization contract to the ten Phase 1 holdouts,
  with Rust, Python, and Node/WASM checking the same largest-fragment outputs.
- Added a same-input 20-repetition SDF record-accounting run for chematic,
  RDKit, and the installed Open Babel CLI, with process boundaries and zero
  failures recorded separately; the report now pins the Open Babel executable
  version in machine-readable metadata.
- Extended the same-input Open Babel accounting runner to V3000 and MOL2,
  including singular-record conversion output and version-pinned reports.
- The same-input runner now records the installed RDKit version alongside the
  Open Babel version in every machine-readable comparison report.
- Preserved the existing implicit Open Babel SDF comparison when no executable
  flag is supplied, while keeping new V3000/MOL2 comparisons opt-in.
- Extended the Explorer browser smoke with the stable empty-paste error and
  recovery path before the cancellation/display-cap workload, and made the
  Cancel assertion wait for its asynchronous hidden transition.
- Added an Explorer browser smoke contract for malformed pasted records,
  cancellation, and the 2,000-record display cap, with a Chromium-local run
  and a CI matrix hook for the existing browser lanes.
- Added a checked `ReactionDocument::from_json_str` boundary and routed the
  Python/WASM RXN-document writers through it, preventing invalid typed
  documents from bypassing component and metadata validation.
- CDXML document parsing now rejects unmatched and nested page elements with
  typed invalid-document errors instead of silently changing page structure.
- Semantic polymer JSON parsing now rejects repeat endpoint indices outside
  the `u32` contract instead of truncating them during deserialization.
- Rich reaction components now expose optional, validated atom-map identities
  (map number plus component-local atom index), populated by RXN/SMILES-derived
  documents without breaking older authored JSON.
- V2000 MOL/RXN atom-map fields now round-trip for three-digit values and fail
  closed with a typed loss for values the fixed-width RXN dialect cannot hold.
- Added an experimental, bounded periodic-structure geometry fingerprint with
  deterministic stored-representation semantics, provenance, JSON Schema, and
  WASM validation/fingerprint support.
- Added opt-in deterministic SVG publication preflight diagnostics for Rust,
  JSON, and WASM consumers. The conservative contract reports stable paths for
  invalid geometry, clipping, overlaps, crossings, degenerate bonds, and
  resource limits, and includes a reproducible input/style fingerprint.
- Normalized parser-side aromatic E/Z direction stashes during canonical orbit
  coloring without changing the physical aromatic bond order. The remaining
  three held-out representation-dependent outputs continue to fail closed.
- Added a versioned, bounded NMR spectrum interchange contract with finite peak
  validation, explicit normalization, opaque vendor metadata, stable
  diagnostics, JSON Schema, and a WASM validation entry point.
- Expanded the streaming malformed-input safety corpus to ten cases for each
  of the ten runner formats, including Extended XYZ (100/100 negative cases),
  preserving the separate oversized-input and 20-case gzip controls.
- Added a bounded, explainable genotoxicity-reactivity slice: source-referenced
  PubChem structure fixtures, epoxide/aziridine/Michael-acceptor checks, and a
  deterministic two-site electrophile heuristic with topological spacer data.
  This remains structural triage only, not a biological predictor or score.
- Promoted the connectivity-ordered 3D coordinate engine to the default
  `generate_coords` path, covering fused-ring seam and chain-bridged ring-island
  layouts with deterministic ring-entry placement.
- Extended the Python file-backed streaming batch contract from SDF to plain
  XYZ and Extended XYZ trajectories. The new iterators preserve input order
  and bounded batch boundaries, support cancellation, and expose deterministic
  progress manifests with rejected-frame counts.
- Extended the common Rust streaming benchmark and safety gate to Extended
  XYZ, including plain/gzip controls, ten malformed cases, and an oversized
  post-decompression limit check.
- Added the same-input Extended XYZ record/failure contract report for
  chematic and RDKit, with source-byte accounting and explicit parser-boundary
  notes; it is not a same-process parity or speed claim.
- Added an Extended XYZ parse fixture to the shared Rust/Python/Node/WASM
  contract manifest, checking coordinates, lattice, typed per-atom properties,
  and frame metadata in each binding.
- Added the typed RXN document round-trip to the shared binding contract,
  covering authored component roles and SMILES across Rust, Python, and
  Node/WASM while retaining the legacy V2000 loss boundary.
- Extended the deterministic gzip streaming record/failure contract to
  Extended XYZ, keeping compressed Rust input bytes separate from decompressed
  RDKit frame-parser bytes.
- Strengthened the shared ECFP4/MACCS binding contract with exact outputs for
  four fixtures across Rust, Python, and Node/WASM; the full held-out parity
  corpus remains separate.
- Added a WASM `topo_path_bitvec` entry point and exact native topological-path
  bit fixtures across Rust, Python, and Node/WASM.
- Extended the shared native fingerprint contract with exact topological-torsion
  bit positions across Rust, Python, and Node/WASM.
- Added a WASM `rdkit_torsion_bitvec` entry point and exact four-fixture bits
  for the RDKit-compatible hashed torsion operation.

## [1.0.9] - 2026-09-07

- Expanded the streaming format safety gate to 27 malformed cases across nine
  supported formats, with explicit oversized-input and gzip limit coverage.
- Added arithmetic and version consistency validation for held-out parity
  reports, keeping measured, not-measured, and mismatch states distinct.
- Hardened Node/WASM batch and screening JSON contracts for malformed records,
  continuation, delimiter errors, and resource limits.
- Synchronized the browser demo's web-target WASM artifact, cache-buster, and
  displayed version with the 1.0.9 release candidate.

## [1.0.8] - 2026-09-06

- Added a versioned cross-binding fingerprint contract for core ECFP4 and
  MACCS outputs, including bit width, packed-byte shape, bit order,
  configuration, sparse/count semantics, and implementation provenance.
- Added Rust, Python, and Node/WASM contract tests for the shared fingerprint
  shape and non-empty output boundary. Held-out value parity and explanation
  contracts remain separate follow-up gates.
- Added `PeriodicStructure::identity_bytes()` with a versioned deterministic
  binary representation for exact cache keys and content-addressed storage.
  The representation preserves validated lattice/site/species/occupancy and
  label data without adding a digest dependency.

- Updated the release documentation and package metadata to the v1.0.8
  publication boundary.

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

[Unreleased]: https://github.com/kent-tokyo/chematic/compare/v1.0.15...HEAD
[1.0.15]: https://github.com/kent-tokyo/chematic/compare/v1.0.14...v1.0.15
[1.0.14]: https://github.com/kent-tokyo/chematic/compare/v1.0.13...v1.0.14
[1.0.13]: https://github.com/kent-tokyo/chematic/compare/v1.0.12...v1.0.13
[1.0.12]: https://github.com/kent-tokyo/chematic/compare/v1.0.11...v1.0.12
[1.0.11]: https://github.com/kent-tokyo/chematic/compare/v1.0.10...v1.0.11
[1.0.10]: https://github.com/kent-tokyo/chematic/compare/v1.0.9...v1.0.10
[1.0.9]: https://github.com/kent-tokyo/chematic/compare/v1.0.8...v1.0.9
[1.0.8]: https://github.com/kent-tokyo/chematic/compare/v1.0.7...v1.0.8
[1.0.7]: https://github.com/kent-tokyo/chematic/compare/v1.0.6...v1.0.7
[1.0.6]: https://github.com/kent-tokyo/chematic/compare/v1.0.5...v1.0.6
[1.0.5]: https://github.com/kent-tokyo/chematic/compare/v1.0.4...v1.0.5
[1.0.4]: https://github.com/kent-tokyo/chematic/compare/v1.0.3...v1.0.4
[1.0.3]: https://github.com/kent-tokyo/chematic/compare/v1.0.2...v1.0.3
[1.0.2]: https://github.com/kent-tokyo/chematic/compare/v1.0.1...v1.0.2
[1.0.1]: https://github.com/kent-tokyo/chematic/compare/v1.0.0...v1.0.1
[1.0.0]: https://github.com/kent-tokyo/chematic/releases/tag/v1.0.0
[0.89.0]: https://github.com/kent-tokyo/chematic/compare/v0.49.0...v0.89.0
