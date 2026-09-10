# chematic roadmap

> Revised 2026-09-11. The current release is v1.0.12. The workspace version is
> 1.0.12; checked-in benchmark artifacts retain their measured version.

The detailed roadmap and completed gate-by-gate evidence through 2026-09-05 is
retained in
[`docs/archive/roadmap-through-2026-09-05.md`](docs/archive/roadmap-through-2026-09-05.md).
This file is the active plan.

## North star

Make chematic a dependable, embeddable cheminformatics runtime for Rust,
Python, JavaScript/WASM, and AI agents. Compatibility claims must be named and
measured; unsupported work must fail explicitly; untrusted input must remain
bounded.

RDKit remains the broadest and most mature reference, Open Babel the format
conversion reference, CDK a reaction/SMARTS/QSAR reference, sdfrust a Rust
dataset reference, and kekule a polymer/modeling reference. COSMolKit is not
part of the active comparison program.

## Competitive response: official RDKit.js

RDKit now maintains an official JavaScript/WASM distribution path through
`@rdkit/rdkit`, with the MinimalLib WASM artifact and TypeScript-facing package
surface released from the RDKit source tree. This removes much of the former
browser-installation friction. The official distribution is therefore the
primary browser comparator for the next candidate; the existing v1.0.10
comparison artifacts remain historical measurements, not a claim about every
future RDKit.js release. See the [RDKit MinimalLib README](https://github.com/rdkit/rdkit/blob/master/Code/MinimalLib/README.md)
and the [rdkit-js repository](https://github.com/rdkit/rdkit-js).

The competitive objective is not feature-count parity. It is to make a named,
reproducible browser workload measurably better while retaining explicit
failure boundaries:

1. **P0 — official comparison gate.** Pin the first stable `@rdkit/rdkit`
   release, record package/WASM/type-definition digests, and rerun the same
   corpus and harness for raw/gzip size, cold initialization, SMILES parse and
   write, ECFP4/Morgan, peak RSS, and browser-engine timing. Keep cold-start,
   steady-state throughput, and memory as separate dimensions.
2. **P1 — correctness moat.** Add negative-charge resonance CIP cases,
   including `[CH-]1C=CC=C1`, with atom-order shuffles and equivalent-operation
   checks. A result is accepted only when labels are invariant or the API
   returns a documented fail-closed outcome; RDKit agreement alone is not the
   invariant.
3. **P1 — browser developer experience.** Compare the smallest useful API
   surface: initialization, typed errors, cancellation, limits, JSON
   serialization, and generated TypeScript declarations. Publish a single
   install-to-first-result example for both packages and measure its startup
   path.
4. **P2 — interoperability wedge.** Add V3000 `SGROUP`/`COLLECTION` ordering
   round trips and preserve unsupported richness. Keep strict parsing separate
   from a future relaxed preview mode; incomplete editor input must never be
   silently accepted by the normal API.
5. **P2 — positioning gate.** Update public comparison claims only after the
   preceding measurements are checked in. Claims may say “smaller”, “faster”,
   or “more compatible” only for a pinned version, corpus, runtime, hardware,
   and failure policy.

The next candidate should ship only the local portions of P0–P2 that have
reproducible evidence. Missing stable packages, browser engines, or
   cross-platform runners remain environment-dependent evidence gaps rather than
   successful competitor claims.

## Current candidate

Competitive-response gate status (2026-09-10):

- [x] Official `@rdkit/rdkit` Node and Playwright Chromium comparison, pinned
  package/artifact digests, 1,000-row ECFP4/Morgan parity, and separate size,
  startup, throughput, and memory evidence.
- [x] Negative-charge resonance CIP corpus with four checked-in cases and
  deterministic atom-order permutations; labels remain invariant in the
  declared scope.
- [x] Browser API/package contract covering malformed input, stable limits,
  resumable batches, serialization, TypeScript declaration digest, and a
  generated-artifact install-to-first-result path.
- [x] V3000 `SGROUP`/`COLLECTION` local interoperability slice: opaque SGROUP
  preservation, either input order, SGROUP-before-COLLECTION output, strict
  topology-handle rejection, and generated Web-WASM execution evidence.
- [ ] Broader CIP generated permutations and independent oracle corpus.
- [ ] Exact canonical-SMILES parity, typed SGROUP semantics, and Indigo/RDKit
  cross-engine V3000 fixtures.
- [x] Add typed CDXML text/caption style access for font, size, and alignment,
  with finite-number validation and regression coverage. Full ChemDraw style
  inheritance and rendering semantics remain outside the bounded contract.
- [x] Add typed CDXML page `BoundingBox` access with finite-number and arity
  validation; full page-layout/rendering semantics remain outside the bounded
  contract.
- [x] Add typed CDXML object classification with an explicit unsupported
  variant, preserving unknown raw XML and diagnostics.
- [x] Accept self-closing empty CDXML pages while preserving page order and
  source representation.

### Next priority: compatible similarity search

- [ ] Add a fallible `rdkit_ecfp4` prepared-search profile alongside the fast
  native `ecfp4` default. It must use the existing RDKit-bit-exact Morgan
  implementation, preserve original input indices, and report preprocessing
  failures instead of silently falling back.
- [ ] Make top-k tie ordering deterministic (`score` descending, original
  index ascending) across native and compatible search lanes.
- [ ] Re-run the 4,500-library/500-query top-10 gate with native/native,
  compatible/RDKit, and cross-profile lanes kept separate. Treat the existing
  72.6% native-vs-RDKit overlap as a disclosed cross-profile diagnostic, not
  an RDKit-compatibility score; target at least 99% compatible top-10 recall
  on the valid-input scope and publish failure counts separately.

- [x] Expose the exact provenance-backed contraction graph in semantic
  expansion JSON (`contracted_smiles`) across the existing Python and
  WASM/Node entry points; inference from edited expanded graphs remains
  rejected by design.
- [x] Add typed S-group membership and polymer-linkage definitions to the
  semantic JSON contract, including stable IDs and source-reference checks.
- [x] Add stable-ID S-group kind editing and fail-closed biomolecule boundary;
  first-class biomolecule expansion remains a separate future schema.
- [x] Include S-group and polymer-linkage IDs in expansion provenance mappings
  so downstream editors can address their source atoms deterministically.

Completed through the v1.0.11 release tree (with historical v1.0.10 and v1.0.9
milestones retained below):

- [x] Add issue #509's experimental bounded geometry fingerprint for
  `chematic-crystal`, including deterministic lattice/site/species/occupancy/
  fractional-coordinate hashing, provenance/schema fields, missing-structure
  and resource-limit errors, JSON serialization, and a WASM entry point. It is
  explicitly not symmetry, polymorph, property, or retrieval equivalence.
- [x] Resolve issue #508 with an opt-in, typed SVG publication preflight in
  `chematic-depict`. The report has stable object paths, deterministic
  diagnostics for invalid geometry, clipping, label/atom overlap, crossings,
  degenerate bonds, and resource limits, plus a deterministic input/style
  fingerprint. Rust and JSON callers share the same implementation, and the
  WASM binding exposes `preflight_smiles_json` with the existing 1 MiB and
  10,000-atom input bounds. Font metrics are explicitly conservative and
  renderer-authoritative; no PDF dependency is added.
- [x] Normalize issue #503's aromatic direction stash in canonical orbit edge
  coloring: a parser-side carrier keeps the physical bond's Aromatic edge class
  instead of becoming a directional edge class. The exhaustive oracle,
  existing stereo checks, and fail-closed stable-key boundary remain green.
  Canonical partition sensitivity now follows the adjacent structural
  exocyclic double-bond endpoint rather than the input-selected aromatic edge;
  equivalent carrier spellings therefore expose the same stereo-sensitive
  neighborhood. The three held-out families still produce two valid spellings
  after this bounded step, so full E/Z convergence remains explicitly open.
  Issue #503 is not marked complete.
- [x] Add issue #495's versioned, bounded NMR interchange contract. Rust types,
  deterministic validation, opaque vendor metadata preservation, the checked-in
  JSON Schema, and WASM `validate_nmr_spectrum_json` are available. This is an
  interchange/visualization model only; vendor parsing, peak picking,
  assignment, prediction, and spectral-accuracy claims remain out of scope.
- [x] Promote issue #486's streaming safety cases into the checked-in
  `validation/streaming_format_safety_cases.json` corpus. The dependency-free
  gate now validates the corpus schema and exact nine-format/four-case shape
  before running the 36 negative cases; lenient CML/CDXML/PDB boundaries remain
  explicit rather than being reported as chemistry-malformed parity.
- [x] Add issue #487's same-input SDF record/failure contract across chematic,
  RDKit, and Open Babel. The 2026-09-08 report records 40/40 valid records,
  zero failures, fixture digest, and each parser/process boundary; it does not
  promote those differently scoped runs into a speed or semantic-parity claim.
- [x] Add issue #488's same-input XYZ record/failure contract across chematic
  and RDKit. The 2026-09-08 report records 40/40 valid records, zero failures,
  identical 3,220 input bytes, fixture digest, and the file-backed Rust versus
  Python block-parser boundary; it does not claim same-process parity or speed
  equivalence.
- [x] Add issue #490's same-input V2000 MOL record/failure contract across
  chematic and RDKit. The 2026-09-08 report records 40/40 valid records, zero
  failures, identical 12,660 source bytes, fixture digest, and the file-backed
  Rust versus Python block-parser boundary; it does not claim same-process
  parity or speed equivalence.
- [x] Add issue #491's same-input V3000 record/failure contract across chematic
  and RDKit. The 2026-09-08 report records 20/20 valid records, zero failures,
  identical 5,960 source bytes, fixture digest, and the materialized Rust versus
  Python block-parser boundary; it does not claim same-process parity or speed
  equivalence.
- [x] Add issue #492's same-input MOL2 record/failure contract across chematic
  and RDKit. The 2026-09-08 report records 20/20 valid records, zero failures,
  identical 7,220 source bytes, fixture digest, and the materialized Rust versus
  Python block-parser boundary; it does not claim same-process parity or speed
  equivalence.
- [x] Add issue #493's same-input CML record/failure contract across chematic
  and Open Babel. The 2026-09-08 report records 20/20 valid records, zero
  failures, identical 7,660 source bytes, fixture digest, and the materialized
  Rust versus CLI boundary; it does not claim same-process parity or speed
  equivalence.
- [x] Add issue #494's same-input CDXML record/failure contract across chematic
  and Open Babel. The 2026-09-08 report records 20/20 valid records, zero
  failures, identical 5,960 source bytes, fixture digest, and the materialized
  Rust versus CLI boundary; it does not claim same-process parity or speed
  equivalence.
- [x] Add issue #496's same-input mmCIF record/failure contract across chematic
  and Open Babel. The 2026-09-08 report records 20/20 valid records, zero
  failures, identical 13,960 source bytes, fixture digest, and the materialized
  Rust versus CLI boundary; it does not claim same-process parity or speed
  equivalence.
- [x] Add issue #497's same-input PDB record/failure contract across chematic
  and Open Babel. The 2026-09-08 report records 20/20 valid records, zero
  failures, identical 5,720 source bytes, fixture digest, and the materialized
  Rust versus CLI boundary; it does not claim same-process parity or speed
  equivalence.
- [x] Add issue #498's deterministic gzip SDF record/failure contract. The
  2026-09-08 report records 40/40 valid records and zero failures for chematic
  gzip input and RDKit over the decompressed fixture, while keeping compressed
  and decompressed byte stages distinct; it does not claim same-process parity
  or speed equivalence.
- [x] Extend issue #499's deterministic gzip contract to XYZ. The 2026-09-08
  report records 40/40 valid records and zero failures for chematic gzip input
  and RDKit over decompressed frames, while keeping compressed and decompressed
  byte stages distinct; it does not claim same-process parity or speed
  equivalence.
- [x] Expand issue #489's checked-in streaming safety corpus to an initial four
  cases per format (36/36 negative cases), retaining explicit line-limit
  handling for lenient CML/CDXML/PDB readers and the initial 9 oversized plus 4
  gzip gates.
- [x] Extend issue #500's gzip safety gate to XYZ. This established the SDF
  and XYZ gzip controls plus post-decompression input limits, preserving the
  explicit decompressed-limit boundary; issue #504 later generalized it to all
  nine runner formats.
- [x] Add issue #501's same-process SDF parity contract: `SdfFileReader` and
  `parse_sdf_with_limits` now compare molecule atoms, bonds, metadata, and
  property shape on the shared two-record fixture. This closes the bounded
  Rust streaming/materialized slice while cross-language parity remains open.
- [x] Add issue #502's reproducible all-format benchmark matrix:
  `scripts/benchmark_streaming_matrix.py` now measures all ten runner formats
  (including Extended XYZ) in plain and optional gzip modes, records fixture
  hashes and parser boundary metadata, and fails on unexpected record/failure
  counts. A local 20-row run (10 formats × 2 compression modes, one repeat)
  is green; this is measurement coverage, not a cross-engine speed claim.
- [x] Add issue #504's gzip safety coverage for every runner format. The
  dependency-free gate now checks one valid gzip control and one
  post-decompression `max-input-bytes` rejection for each of the nine formats
  (18 gzip cases total), while malformed-corpus breadth remains separate.
- [x] Add issue #505's second malformed case pair for every runner format. The
  corpus and dependency-free gate now cover six cases per format (54/54
  negative cases), while the typed line-limit boundary remains explicit for
  lenient CML/CDXML/PDB readers.
- [x] Add issue #506's third malformed case pair for every runner format. The
  corpus and dependency-free gate now cover eight cases per format (72/72
  negative cases), retaining the explicit typed line-limit boundary for
  lenient CML/CDXML/PDB readers.
- [x] Complete issue #507's fourth malformed case pair for every runner
  format. The corpus and dependency-free gate now cover ten cases per format
  (90/90 negative cases), while the nine oversized and 18 gzip controls remain
  green; lenient CML/CDXML/PDB inputs continue to use the explicit line-limit
  safety boundary.
- [x] Complete the bounded first slice of issue #303: retain explainable
  epoxide, aziridine, and Michael-acceptor findings; add source-referenced
  PubChem structure fixtures; and report at most 64 deterministic pairs of
  independent electrophilic sites with a topological spacer estimate. No
  biological label, genotoxicity score, 3D claim, or external validation is
  implied; broader category coverage remains open.

- [x] Close #210's five named legacy-coordinate UFF stereo-rescue residuals.
  Every returned geometry is finite, bond-sane, and independently checked
  against declared stereo.
- [x] Exceed the additional local 1.10x hot-path target: canonical SMILES
  1.115x, file-backed SDF read 1.130x, and V2000 SDF serialization 3.042x at
  the recorded medians.
- [x] Preserve output-byte, malformed-input, workspace-test, Python-binding,
  rustfmt, and warnings-as-errors clippy gates for those changes.
- [x] Reorganize public documentation: concise README files, release-focused
  CHANGELOG, active ROADMAP, corrected publication status, and archived
  detailed development evidence.
- [x] Close the bounded #246 bridged-bicyclic 2D layout residual: anchored
  ring placement now evaluates already-placed shared atoms and keeps every
  named regression bond near `BOND_LEN`; the full depiction test suite passes.
- [x] Add the bounded `PeriodicStructure::identity_digest()` API over the
  versioned exact-identity bytes, with a dependency-free SHA-256 implementation
  and deterministic/change-sensitive regression tests; this remains an exact
  stored-representation key, not symmetry canonicalization.
- [x] Implement #478's deterministic `chematic-smiles` batch canonicalization
  API: lazy input-order results, reusable parser limits, retained source text,
  and per-record accepted/rejected diagnostics without batch aborts.
- [x] Add the #478 newline-delimited `BufRead` adapter with the same result
  contract; parse failures remain ordered records and underlying I/O failures
  remain explicit `io::Error` items.
- [x] Add the #478 compiled exact-identity index over
  `canonical_smiles_stable_key()`, preserving duplicate input positions and
  fail-closed rejected records.
- [x] Add a shared #478 batch-canonicalization fixture consumed by Rust,
  Python, and Node/WASM bindings, including deterministic accepted/rejected
  ordering and canonical output expectations.
- [x] Wrap the #478 Python/WASM batch JSON APIs in the common versioned
  partial-result envelope (`schema_version`, `operation`, `status`, and
  `record_count`) while retaining per-record diagnostics.
- [x] Extend the SDF/MOL/XYZ streaming benchmark runner with explicit input,
  record, line, frame, and atom limits, and record effective limits in JSON;
  V3000/MOL2/CML/CDXML/mmCIF parser coverage is now included as explicitly
  materialized one-shot rows; true same-process streaming parity remains open.

The performance figures are source-level measurements recorded before the
release artifact was built; they remain scoped to their named corpus and
configuration.

## Priority phases

| Priority | Goal | Exit evidence | Cost |
|---|---|---|---|
| P0 | Trust and measurement | Version-pinned schema, corpus hashes, typed statuses, offline validation | Critical / light |
| P1 | Interchange throughput | Same-input SDF/MOL/XYZ measurements, memory/error gates, loss-preserving streaming | Critical / heavy |
| P2 | Identity and ML primitives | Stable canonical key, held-out descriptor/fingerprint reports, explainable output | Critical / heavy |
| P3 | Portable production surface | One Rust/Python/Node/WASM fixture contract, clean-install and browser evidence | High / heavy |
| P4 | Chemistry workflows | Reaction/SMARTS/standardization provenance and ambiguity contracts | High / heavy |
| P5 | 3D and materials | Soundness, unit/frame round trips, per-class quality and failure rates | Medium / heavy |
| P6 | Ecosystem durability | Reproducible dashboards, migration paths, contributor corpus policy | Medium / external |

Complete the highest-priority reproducible local slice before starting lower-
priority breadth. A feature is not complete until implementation, tests,
documentation, and required measurement agree.

## P0 — Trust and measurement

- [x] Separate capability claims from measured results and preserve
  `unsupported`, `failure`, and `not_measured` as distinct states.
- [x] Pin benchmark protocols, corpus hashes, environment metadata, and
  historical result records.
- [x] Maintain a capability matrix for RDKit, Open Babel, CDK, sdfrust,
  kekule, and chematic.
- [x] Exclude COSMolKit from the active comparison scope.
- [x] Add a two-build Criterion null-control: independently compile `main` and
  `main-null` and compare them through the same process-level pipeline, so
  build/codegen variance is measured before a real regression can block.
- [x] Check in the #70 synthetic calibration contract for +5%, +10%,
  sub-threshold build noise, and contamination routing; hosted real-run
  evidence remains a separate gate.
- [x] Add versioned machine-readable release metadata and a tag-driven GitHub
  Release attachment. The schema, versioned raw JSON, validator, and historical
  benchmark separation are checked in under `docs/`, `release-metadata/`, and
  `scripts/`.
- [x] Add a scorecard validator that rejects stale release versions, missing
  corpus/configuration metadata, and claims derived from unsupported rows;
  `scripts/validate_scorecard.py` is dependency-free and fail-closed.

## P1 — Interchange throughput and safety

- [x] Add file-backed SDF/MOL/XYZ benchmark fixtures and a resumable runner.
- [x] Keep the RDKit block-parser comparator explicitly separate from
  chematic's file-backed streaming boundary.
- [x] Record a 2,000-pass SDF/MOL/XYZ streaming lane with zero malformed
  fixture failures and explicit cross-engine boundary notes.
- [x] Keep the dated benchmark record index complete and link-checked with a
  dependency-free CI gate; historical records remain versioned snapshots.
- [x] Make streaming benchmark generators derive `target_version` from the
  workspace manifest, preventing future candidate reports from inheriting the
  v1.0.9 metadata while preserving existing historical outputs.
- [x] Extend the common benchmark to V2000/V3000, XYZ/Extended XYZ, MOL2, CML,
  CDXML, PDB/mmCIF, and gzip, including bounded malformed and oversized inputs.
  The dependency-free gate was revalidated on 2026-09-10 and passes 800
  malformed-input attempts (800 unique payloads), one oversized case for every
  runner format (10/10), plus ten gzip controls and ten post-decompression
  limit cases. The
  safety runner enforces the per-format minimum so aggregate totals cannot
  hide a format-specific gap. Its twelve checked-in base cases per format now
  also carry a versioned category manifest, and the gate requires every
  format-specific malformed category to be represented before supplemental
  parser-path cases run. The generated parser-entry wave now also requires
  eight declared six-case entry families per format, including format-shaped
  truncation and numeric-or-vocabulary corruption; the result is recorded
  in `validation/results/streaming-parser-entry-categories-v1.0.10.json`.
  The generated wave also records a non-empty typed failure variant for all
  80 format/category combinations in
  `validation/results/streaming-parser-entry-failure-kinds-v1.0.10.json`.
  Exhaustive malformed-corpus coverage remains a separate open item.
- [x] Record the current v1.0.10 bounded streaming safety gate with 800
  malformed, 10 oversized, and 20 gzip cases passing across the ten runner
  formats. The 120-case base corpus also has reproducible error-variant
  taxonomy output for all ten formats; exhaustive malformed-corpus and parser-
  state coverage remain separate open items.
- [x] Add bounded Rust streaming batch APIs with cancellation, pull-based
  backpressure, deterministic ordering, and an explicit partial-result
  manifest (`SdfBatchReader`); cross-language streaming parity remains open.
- [x] Add a versioned partial-result manifest envelope to every CLI batch
  operation. It preserves input-order records and exposes the operation,
  status, record count, and effective input limits; cancellation and
  backpressure remain separate open gates for the streaming API.
- [x] Harden the Python file-backed SDF batch iterator with bounded batch
  sizes, explicit cancellation, deterministic progress manifests, and lazy
  input-order emission; full cross-language streaming parity remains open.
- [x] Extend the Python file-backed batch contract to plain XYZ and Extended
  XYZ frames. Both bindings preserve input order and batch boundaries, expose
  cancellation and deterministic progress manifests, count rejected frames,
  and recover after malformed content at an unambiguous count-line boundary;
  WASM now has bounded resumable SDF/XYZ/Extended XYZ manifest slices plus the
  same bounded recovery rule. Browser-artifact parity and broader
  error-recovery semantics remain open.
- [ ] Measure only equivalent operations against installed RDKit and Open
  Babel versions on identical inputs; report sdfrust separately. A same-
  process Python SDF semantic contract now confirms record/failure counts and
  canonical SMILES parity for 20 repetitions. Two malformed controls produce
  zero accepted records and one failure in both engines. V2000 MOL block
  parsing now also has a same-process 20-repetition contract with zero
  record/failure or canonical-SMILES mismatches, and V3000 MOL now has the
  same contract with zero record/failure or structural-signature mismatches;
  MOL2 now has the same-process contract with zero record/failure or
  structural-signature mismatches (charge/SMILES parity intentionally not
  claimed); XYZ now has the same-process 20-repetition frame/coordinate
  contract with zero mismatches; Extended XYZ now has the same-process
  coordinate/frame contract with zero mismatches while property parity remains
  schematic-only; PDB now has the same-process coordinate/signature contract
  with zero valid-input mismatches. The Rust core now also provides an
  opt-in fixed-column `parse_pdb_atoms_strict` validator while the Python/WASM
  compatibility APIs remain intentionally lenient; Python now exposes the same
  opt-in boundary as `from_pdb_strict`, and the WASM source now exposes
  `mol_from_pdb_strict`. Generated-artifact/runtime verification, broader
  malformed-corpus coverage, and full binding-level strict-mode coverage
  remain open. Broader formats,
  file-backed reports now cover SDF, V2000/V3000 MOL, MOL2, XYZ, Extended XYZ,
  CML, CDXML, mmCIF, and PDB with explicit parser/process boundaries. Open
  Babel remains recorded separately because only its CLI boundary is
  installed; the 20-repetition SDF, V3000, and MOL2 lanes now record matching
  counts and zero failures for all three available engines with those
  boundaries explicit.
  The Rust-only common matrix now records 20 plain/gzip rows across all ten
  runner formats with normalized fixture-relative paths. A 20-repetition
  same-input matrix now records matching record counts and zero failures for
  all ten formats across chematic and the installed RDKit/Open Babel lanes
  where available; the reproducible result is
  `validation/results/cross-engine-matrix-v1.0.10.json`. Open Babel is CLI-only
  here, so its rows are subprocess boundary evidence rather than same-process
  speed evidence. This is contract evidence only because the parser and
  process boundaries are not equivalent. The same-process semantic contracts
  currently cover SDF, V2000/V3000 MOL, MOL2, XYZ, Extended XYZ, PDB, and CDXML;
  `scripts/check_same_process_contracts.py` validates the complete eight-format evidence bundle, current-version metadata, malformed-case reports, and explicit non-ranking timing boundary;
  Same-process timing fields are now required for all eight semantic-contract
  reports by `scripts/check_same_process_contracts.py`; the dedicated SDF
  context run on 2026-09-10 used a locally built v1.0.10 Python extension and
  RDKit 2025.09.3, with 40 alternating rounds over the same 633-byte fixture,
  2/2 records, zero failures, and P50 values of 14,176.7 versus 374,594.2
  ns/parse. These are bounded context results, not a general speed ranking.
  Equivalent timing across all ten formats and full malformed/signature parity
  remain separate open claims.
  Deterministic gzip record/failure contracts now cover SDF, XYZ,
  Extended XYZ, V3000, MOL2, CML, CDXML, mmCIF, and PDB with
  compressed/decompressed byte stages explicit; same-
  process cross-engine equivalence remains open. The 2026-09-10 v1.0.10
  refresh passes all available lanes with zero failures and records current
  RDKit/Open Babel versions; its throughput fields remain non-ranking context
  because parser and process boundaries are not equivalent.

## P2 — Identity and ML primitives

- [x] Make `canonical_smiles_stable_key()` the only recommended dedup/cache
  path; it fails closed when stability is not proven, including coupled E/Z
  systems using aromatic direction stashes.
- [x] Keep native and RDKit-compatible fingerprint modes separate in API names
  and documentation.
- [x] Exclude Spectrophores and Issue #464's proposed replacement pending
  independent patent/FTO review.
- [ ] Finish canonical atom-order and E/Z invariance for the supported domain.
  The 2026-09-09 offline standardized NCI-5k, descriptor-census, and ChEMBL
  idempotency lanes all pass (3/3, 637.30s); this strengthens round-trip
  evidence but does not close atom-order invariance or the three held-out
  aromatic-stash representation residuals.
- [x] Add a charged-conjugated-ring atom-order regression probe for
  `[CH-]1C=CC=C1` and two equivalent spellings; all three canonicalize to one
  output. This is a targeted invariant case, not completion of the full
  canonical corpus gate.
- [x] Prevent ring-closure close-side carrier selection from erasing an E/Z
  marker; the regression fixture and focused canonical suite are green. This
  is a bounded residual fix, not completion of the full corpus gate.
- [x] Resolve coupled E/Z carrier choices jointly across shared physical
  bonds, with 19 permutation-invariance fixtures; retain fail-closed
  abstention when no conflict-free assignment is proven.
- [x] Re-audit the committed 5,000-line carrier corpus with 64 seeded
  relabelings per coupled molecule: 28 coupled components, all size 2, and 0
  correspondence failures. The audit now reproducibly exposes 3 residual
  molecules with two canonical outputs; those three inputs and their observed
  output pairs now have held-out Rust regression coverage. The Rust gate also
  replays 256 deterministic atom relabelings per residual, verifies that the
  relabeling-only axis preserves one output, then combines it with the two
  measured aromatic-stash spellings to retain exactly two outputs;
  `canonical_smiles_stable_key()` remains fail-closed for each. Keep #149 open
  until their aromatic carrier/stereo traversal is normalized.
- [x] Extend the #149 residual oracle gate across the corpus spelling and both
  observed aromatic-stash output variants for all three held-out molecules.
  The orbit-pruned search agrees with the unpruned exhaustive oracle for all
  nine spellings; the two representation-dependent canonical outputs remain
  intentionally distinct and fail-closed rather than receiving an
  index-based winner.
- [x] Add a bounded default-Hückel fused/non-alternant fallback for
  all-carbon odd/odd envelopes, with azulene regression coverage; keep the
  broader `RdkitLike` model separately gated.
- [x] Add a held-out CIP parity gate requiring zero wrong confident labels and
  zero regressions in the non-phosphorus scope (140/140 current cases), and
  fail closed for all 15 representation-unstable phosphorus rows.
- [x] Freeze descriptor/fingerprint shape, sparse/count, configuration,
  provenance, and explanation contracts. The shared contract now covers the
  RDKit-exact ECFP4 detail operation (folded bytes, sparse counts, and raw /
  folded atom-radius provenance) across Rust, Python, Node, and WASM; held-out
  parity reports remain a separate open item.
- [x] Freeze the core bit-packed ECFP4/MACCS shape, configuration, bit order,
  and implementation provenance in `validation/cross_binding_contract.json`,
  with Rust/Python/Node/WASM shape tests. Sparse/count and explanation
  contracts remain open for the held-out parity phase.
- [x] Record a local #372 canonical-symmetry proxy with per-fixture nodes,
  orbit tests, leaves, pruning, timing, and old/new correctness checks; keep
  the exact downstream RENKIN witness as an external held-out gate.
- [x] Add descriptor field provenance and a shared Rust/Python/Node/WASM core
  descriptor fixture; run the 5,000-row MW/TPSA/HBD/HBA/heavy-atom lane.
- [x] Add a reproducible 5,000-row same-source descriptor binding-parity lane
  across Rust, Python, and Node/WASM. The committed report records the corpus
  SHA-256, parse-status counts, and exact agreement for MW/TPSA/HBD/HBA and
  heavy atoms; this is a binding contract, not an RDKit accuracy claim.
- [x] Add a reproducible 5,000-row native ECFP4 binding-parity lane across
  Rust, Python, and Node/WASM. The committed report checks the complete
  2,048-bit LSB-first byte representation and keeps RDKit-compatible and
  sparse/explanation parity as separate gates.
- [x] Add the matching 5,000-row native MACCS 166-bit binding-parity lane;
  keep the RDKit key mapping as a separate accuracy gate.
- [x] Add a 5,000-row native `topo_path` binding-parity lane across Rust,
  Python, and Node/WASM; keep RDKit-compatible path parity separate.
- [x] Add a 5,000-row native torsion binding-parity lane across Rust, Python,
  and Node/WASM; keep RDKit-compatible torsion separate.
- [x] Replay the explicit standardization profile on 5,000 rows across Rust,
  Python, and Node/WASM; keep the ten expected-identity holdouts as the
  separate external accuracy boundary.
- [x] Add a 5,000-row RDKit-compatible ECFP4 binding-parity report across
  Rust, Python, and Node/WASM, comparing typed preprocessing status and exact
  2,048-bit values; keep sparse/count/bitInfo parity separate.
- [x] Add the corresponding 5,000-row RDKit-compatible hashed torsion
  binding-parity report across Rust, Python, and Node/WASM.
- [x] Add a 5,000-row RDKit-compatible ECFP4 sparse identifier/count parity
  report across Rust, Python, and Node/WASM; keep raw/folded bitInfo separate.
- [x] Add a 5,000-row RDKit-compatible ECFP4 folded bitInfo parity report
  across Rust, Python, and Node/WASM.
- [x] Add a 5,000-row RDKit-compatible ECFP4 raw bitInfo parity report across
  Rust, Python, and Node/WASM.
- [x] Add a 5,000-row RDKit-compatible path fingerprint binding-parity report
  across Rust, Python, and Node/WASM; keep native `topo_path` separate.
- [x] Add held-out parity reports for Morgan/ECFP, MACCS, topological,
  torsion, descriptors, and standardization across Rust/Python/WASM. The
  v1.0.8 coverage ledger remains validated by
  `scripts/validate_held_out_parity_manifest.py`, including report-internal
  row/count arithmetic: Morgan/ECFP and the five core descriptors have named
  reports, and MACCS now has a Rust-CLI report
  (4,714/5,000 exact bit-set matches, 94.28%). The RDKit-compatible
  topological `rdkit_rdk` lane is also measured by the Rust CLI (5,000/5,000
  exact bit-set matches). The RDKit-compatible hashed topological torsion lane
  is measured by the Rust CLI as well (4,248/5,000 exact bit-set matches,
  84.96%); native `topo_path` remains explicitly `not_measured`, and
  exact ECFP4 bit positions, MACCS bytes, native `topo_path` bit positions, and
  native torsion and RDKit-compatible torsion bit positions are now checked
  across Rust, Python, and Node/WASM for four shared fixtures; cross-binding
  parity for the held-out measured operations is recorded in the v1.0.9
  `validation/held_out_parity_manifest-v1.0.9.json` and passes its validator.
  External RDKit accuracy, native-vs-compatible semantic equivalence, and the
  four-fixture bit-position contract remain separate gates. The explicit
  Rust CLI standardization profile is measured against the 10 committed Phase
  1 holdouts (8/10 expected fragment identities, 80%); Rust, Python, and
  Node/WASM now share and reproduce the ten returned largest-fragment outputs.
  Two acetate cases still differ from the external expected identities after
  charge neutralization and remain disclosed mismatches.

### Performance acceleration track

Prioritize measured, semantics-preserving speedups before invasive force-field
work. Expected multipliers are hypotheses until a pinned release-mode lane
records them.

- [x] Precompute MAP4 circular environment hashes once per atom/radius instead
  of rescanning and sorting the whole molecule inside every atom-pair loop;
  preserve the shingle set and MinHash output byte-for-byte.
- [x] Make `bulk.descriptors_array(smiles, columns)` execute only the selected
  descriptor dependency groups, including shared ring/logP/pKa values and
  precomputed ADMET/filter formulas; retain the full `bulk.descriptors()` API
  as the explicit all-fields path.
- [x] Add release-mode scaling benchmarks for MAP4 and selected descriptor
  columns (3, 8, and all fields), with Python-visible allocation and
  deterministic output-digest checks; native allocation accounting remains
  explicitly outside `tracemalloc` scope.
- [x] Add a prepared fingerprint index for repeated database queries; exact
  top-k search now reuses database fingerprints while preserving original
  indices through the Python binding.
- [x] Add a parallel row-wise Tanimoto matrix path with serial output parity;
  tiled/threshold-aware search remains open.
- [x] Share descriptor topology context so heavy-atom extraction and topology
  membership are computed once across related Wiener/Kappa/Chi groups; single
  selected columns retain lazy scalar computation.
- [x] Add a shared distance descriptor bundle for AutoCorr2D, Moran, and
  Geary; retain lazy single-family APIs and exact scalar parity.
- [x] Prepare MMFF94 nonbonded parameter combinations and electrostatic charge
  products once per topology; preserve finite-difference output parity and
  record a release-mode energy/minimization lane. Remaining analytic terms,
  electrostatic neighbor handling, and full coordinate-dependent neighbor-list
  integration remain separate experimental gates.
- [x] Parallelize independent prepared MMFF94 finite-difference atom probes
  for molecules with at least 16 atoms while retaining the allocation-light
  sequential path for small molecules; retain exact scalar parity.
- [x] Apply the same bounded finite-difference probe parallelism to the public
  steepest-descent MMFF94 path; keep its existing step size and convergence
  semantics unchanged.
- [x] Optimize canonical rank normalization and V2000/SDF hot paths with
  frozen-binary output parity and resumable alternating A/B evidence. On the
  v1.0.6 local-source lane, median paired speedups are canonical 1.176x,
  SDF graph/property read 1.180x and reused-buffer serialization 1.419x;
  [protocol and load caveats](benchmarks/2026-09-05-hotpath-110.md).
- [x] Replace SMILES ring hashing with bounded direct indexing and inline
  partner storage; validate all labels, reuse and stereo-buffer spill.
- [ ] Confirm at least 1.10x additional SMILES parse throughput under a
  low-noise controlled run; current paired median is 1.034x, not a pass. The
  maintainer explicitly abandoned this stretch target on 2026-09-08; retain
  it as historical context, not as a v1.0.10 release gate.
- [ ] Replace MD/UFF/MMFF94 finite-difference force paths with prepared
  topology, analytic gradients, and bounded neighbor lists. This is a heavy
  experimental-3D task and requires energy/gradient/stereo soundness gates.
- [x] Add the bounded UFF prepared-energy topology slice: bond, angle, and
  graph-based vdW exclusion terms are prepared once per minimization and
  reused by every energy/gradient evaluation. The UFF minimizer now uses the
  analytic bond/angle/vdW gradient, verified against a finite-difference
  reference; MMFF94 vdW now has a bounded coordinate-dependent 10 Å
  neighbor-list slice with exact cutoff parity, while broader MMFF94
  electrostatic and full minimizer integration remain open.
- [x] Add the bounded MMFF94 prepared gradient slices: the reusable energy
  model now exposes separately verified analytic bond/angle, stretch-bend,
  vdW/electrostatic, non-singular torsion, and non-singular Wilson
  out-of-plane gradients while the existing L-BFGS finite-difference path
  remains unchanged. All five slices and their bounded aggregate are verified
  against central differences; OOP branch-point handling is fail-closed; an
  opt-in analytic L-BFGS smoke is green; coordinate-dependent neighbor-list
  work remains open. MMFF94 vdW now uses a coordinate-dependent 10 Å cell-list
  cutoff with exact key/energy parity against the all-pair reference; broader
  electrostatic neighbor-list, triclinic periodic, and performance gates remain
  open. The public non-periodic direct Coulomb cutoff path now also uses a
  finite-coordinate cell list with strict-boundary parity; see
  `validation/results/ewald-coulomb-cutoff-neighbor-list-v1.0.10.json`.
- [x] Add the bounded opt-in MMFF94 electrostatic 10 Å cell-list path. Its
  coordinate-dependent candidate set, strict cutoff boundary, negative-coordinate
  handling, energy parity, and central-difference gradient parity are covered
  by package tests and recorded in
  `validation/results/mmff94-electrostatic-neighbor-list-v1.0.10.json`.
  The opt-in analytic L-BFGS path now switches cutoff energy and gradient
  together; the compatibility-default all-pair path and public minimizer remain
  unchanged. Broader production integration and performance gates remain open.
- [ ] Replace periodic neighbor all-pairs enumeration and Ewald real-space
  all-pairs evaluation with validated cell-list/cutoff paths.
- [x] Add and verify the bounded non-periodic Ewald real-space cell-list path;
  the PME real-space term now applies the configured cutoff and matches the
  direct reference on negative-coordinate and cell-boundary fixtures. The
  periodic crystal neighbor path and broader numerical/performance gate remain
  open.
- [x] Add and verify the bounded orthorhombic periodic-neighbor cell-list
  path. It preserves the complete `(center, neighbor, image)` result set and
  deterministic ordering against the existing exact search; triclinic cells
  now use the same Cartesian replicated-image path with fractional-delta
  cutoff arithmetic, while oversized replicated indexes retain the validated
  fallback. The 4 hand-selected plus 12 fixed-seed triclinic fixtures and
  their 10-cutoff sweep are recorded in
  `validation/results/periodic-neighbor-triclinic-cell-list-v1.0.10.json`.
- [ ] Optimize symmetry-heavy canonical search and SMARTS VF2 state only after
  exact-output and budget-exhaustion parity gates are expanded. The focused
  2026-09-09 local gate now passes 263/263 (`chematic-smarts` 184/184 and the
  `chematic-smiles` canonical module 79/79), including typed budget exhaustion
  and fail-closed paths; broader canonical invariance and optimization remain
  open.

- [x] Make VF2 query expansion prefer mapped-neighbor and higher-degree atoms,
  preserving exhaustive semantics while resolving the known symmetric PAINS
  negative within the production visit budget; keep the typed budget outcome.

## P3 — Portable production surface

- [x] Provide typed, bounded parser and binding failures with shared
  adversarial fixtures.
- [x] Run Chromium, Firefox, and WebKit smoke/adversarial lanes for the
  published v1.0 boundary.
- [ ] Make Rust, Python, Node, and WASM consume one fixture schema and one
  versioned expected-result manifest for every shared stable operation. The
  manifest now has an explicit, CI-validated inventory of 50 currently shared
  operations with four declared binding surfaces and checked-in test anchors;
  this inventory is deliberately scoped to the covered surface and does not
  yet claim completeness for every public stable API.
  The single-frame XYZ parser is now included as a four-binding contract with
  explicit coordinates and the hydrogen/``heavy_atoms`` counting boundary.
  Its serializer now also has a four-binding semantic round-trip contract;
  coordinate text is compared within the writer's six-decimal precision rather
  than byte-for-byte.
  CML now also has a four-binding semantic round-trip contract; generated
  layout coordinates are intentionally checked only for successful structural
  reparse.
  V2000 MOL now also has a four-binding semantic round-trip contract; layout
  coordinates are intentionally outside the cross-binding equality boundary.
  V3000 MOL has the same four-binding structural round-trip contract, with
  generated layout coordinates likewise outside the equality boundary.
  MOL2 now has a bounded four-binding round-trip contract: Rust/Python verify
  atom/bond counts, while Node verifies the documented SMILES conversion
  boundary rather than claiming a Node molecule-handle API.
  MolJSON now also has an explicit four-binding parse/serialize/parse round-trip
  contract; serialization output is checked semantically rather than byte-for-byte.
  The current source audit covers 251 WASM and 132 Python exported names and
  pins the set in `validation/results/binding_surface_inventory-v1.0.10.json`;
  `scripts/check_wasm_artifact_surface.py` reports exactly four current source
  exports absent from the checked-in Node artifact (`mol_from_cjson`,
  `mol_from_pdb_strict`, `mol_from_pdbqt`, and `reaction_smarts_match`), with
  no stale artifact exports.
  The opt-in strict PDB path now has a shared source-only fixture and passing
  Rust/Python/WASM tests; its evidence is recorded in
  `validation/results/pdb-strict-source-only-contract-v1.0.10.json` while the
  Node generated-artifact lane remains explicitly open.
  SMARTS validity and reaction SMARTS matching now have separate source-only
  Rust/Python/WASM contracts
  contract with four accepted and four rejected cases; evidence is recorded in
  `validation/results/smarts-source-only-contract-v1.0.10.json`. This is not a
  four-binding claim because the checked-in Node artifact has no matching
  `is_valid_smarts` export and is still version-drifted.
  standardization profile is now included in the shared contract and the Node
  contract test passes against the checked-in generated Node artifact; an
  explicit version probe reports 1.0.9 while the workspace is 1.0.10, so this
  is not current-candidate artifact evidence and the rebuild remains open.
  The shared manifest now also covers Extended XYZ parsing, the typed RXN
  document round-trip, and Markush/polymer semantic expansion across Rust,
  Python, and Node/WASM. It now also carries the bounded XYZ/Extended XYZ
  batch-recovery contract, verified by Rust/WASM and a clean Python wheel;
  descriptor fixtures now additionally freeze exact mass and aromatic ring
  count across the same four source surfaces;
  the same descriptor contract now also freezes Hill-notation molecular
  formula output and rotatable-bond counts;
  the checked-in Node artifact exports the same batch functions and passes the
  Node contract smoke; its explicit version probe reports 1.0.9, and the
  broader all-stable-operation manifest remains open. Standardization now also
  has four-binding contracts for InChI generation, fragment/charge/isotope/
  stereo parents, canonical tautomer selection, charge neutralization, and
  largest-fragment selection, budgeted tautomer/super-parent computation, and
  explicit hydrogen addition/removal, and composed super-parent stage reports;
  the current 50-operation result, six source-only contracts, and Python 106/106 evidence are recorded in
  `validation/results/cross-binding-manifest-v1.0.10.json`.
- [x] Extend the shared versioned manifest to cover the core ECFP4/MACCS
  fingerprint shapes and configurations across Rust, Python, and Node/WASM.
- [x] Add versioned deterministic exact-identity serialization for
  `PeriodicStructure`, preserving validated lattice/site/species/occupancy and
  label data without introducing a digest dependency.
- [x] Publish current clean-install, cold-start, throughput, peak-memory, and
  WASM-size evidence with explicit platform/configuration metadata. The
  v1.0.10 macOS arm64 Python candidate is recorded in
  `benchmarks/2026-09-09-clean-install-cold-start-v1.0.10.json` and the
  matching raw/gzip WASM artifact snapshot in
  `benchmarks/2026-09-09-wasm-size-v1.0.10.json`; cross-platform and browser
  memory evidence remain separate open lanes.
- [x] Publish current v1.0.10 clean-install, cold-start, Python SMILES
  throughput, and single-process peak-RSS evidence with platform and build
  metadata. The same clean CPython 3.13 environment also runs all 945 Python
  binding tests with zero failures; cross-platform and browser-memory evidence
  remain open.
- [x] Publish the current v1.0.10 WASM-size snapshot with raw/gzip byte counts,
  SHA-256, toolchain, target, and reproduction metadata; clean-install,
  cold-start, throughput, and browser-memory evidence remain open.
- [x] Add a reproducible Node/WASM comparison gate against an installed
  official RDKit.js package. The 2026-09-09 v1.0.10 artifact run records
  initialization, parse, SMILES write, RDKit-compatible ECFP4/Morgan timing,
  per-process peak RSS, raw/gzip artifact sizes, and 1,000/1,000 exact
  fingerprint matches; browser-engine and future-package reruns remain open.
- [x] Add a same-process paired Node/WASM timing lane for the official RDKit.js
  package. The 2026-09-09 follow-up removes separate-process startup from the
  timed operations and retains explicit fingerprint parity boundaries; shared
  process memory and browser-engine evidence remain separate.
- [x] Publish the v1.0.9 WASM-size snapshot with toolchain, target, digest, and
  reproduction metadata; clean-install, cold-start, throughput, and
  peak-memory evidence remain separate open lanes.
- [ ] Extend browser and agent adversarial cases for cancellation, malformed
  records, limits, and stable JSON errors. The Node/WASM shared contract now
  covers malformed inline records, later-record continuation, empty-delimiter
  rejection, batch-count limits, input-size limits, and screening error
  envelopes. The historical `1.0.9` web-target artifact was regenerated into
  `demo/pkg` and passed the real-WASM `initSync` smoke (success and typed
  timeout); this is not current `1.0.10` artifact evidence. The Explorer
  browser smoke now covers malformed pasted records and cancellation locally
  and is wired into the Chromium/Firefox/WebKit CI
  matrix. The current Node-hosted agent-side pipeline/MCS adversarial run is
  recorded in `validation/results/wasm-agent-adversarial-v1.0.10.json`; the
  broader agent matrix and full browser evidence remain open.
- [x] Prepare a standalone Chromium Node/WASM comparison harness for the
  official RDKit.js package; browser timing/parity evidence remains open until
  a stable headless browser run can be captured.

## P4 — Chemistry workflow depth

- [x] Add typed reaction documents, RXN adapters, loss-preserving CDXML
  commands, and bounded semantic Markush/polymer expansion.
- [x] Preserve agents, coefficients, conditions, atom maps, source mappings,
  and unsupported-richness errors instead of flattening them.
- [x] Publish the issue #473 downstream capability matrix, including the
  explicit first-class nucleic-acid/biopolymer non-support boundary and the
  Rust/Python/WASM/Node entry points.
- [x] Add the bounded #510 reaction-template normalization bridge for
  \`[#N]\`, \`[#N:map]\`, and the common explicit-hydrogen forms \`[#N;H1]\`
  / \`[#N;H1:map]\`, with deterministic aromatic/aliphatic expansion,
  map-preserving normalization, explicit unsupported-primitive errors, and
  regression tests. Full query-aware application semantics remain open.
  The shared reaction-application fixture now exercises fifteen cases in Rust,
  Python, and Node/WASM, including product size and fail-closed
  compound-primitive handling; the reproducible breadth result is recorded in
  `validation/results/reaction-application-breadth-v1.0.10.json`.
- [x] Add a bounded 38-case reaction SMARTS presence-matching contract for
  Rust, Python, and WASM source. It covers elemental queries, negative matches, multiple
  patterns, mapped atoms, matching and mismatching target map labels, later
  embedding selection, a transition-metal element query, disconnected
  reactants with no-target-reuse rejection, explicit single-bond and
  hydrogen-count matching, and an agent section;
  the contract result is recorded in
  `validation/results/reaction-smarts-bounded-contract-v1.0.10.json`.
- [ ] Expand reaction/SMARTS/medicinal-chemistry coverage only after P0-P3
  gates have current evidence.
- [ ] Add curated reaction/query precision, recall, invalid-product, timeout,
  and ambiguity reports.

## P5 — 3D and materials

- [x] Keep 3D generation and MMFF94 Experimental, with typed failure and
  explicit force-field/fallback provenance.
- [x] Expose UFF `worst_bond_length` alongside its independent `sound` result
  in Rust, Python, and WASM; keep `converged` separate from geometry validity.
- [x] Complete the connectivity-ordered placement slice for issues #255 and
  #256: route `generate_coords` through the validated engine, repair fused-ring
  seams and chain-bridged ring islands, and retain deterministic direct-bond
  new-island anchoring. The 33-molecule evaluation is raw sound 33/33,
  deterministic 33/33, and UFF-only success 33/33; this does not close the
  separate MMFF94/UFF force-field residuals below.
- [x] Separate long-running 3D and corpus-scale canonical tests into explicit
  ignored lanes and retain their execution manifest.
- [x] Add crystal composition and materials-format foundations without
  conflating periodic structures with molecular bond graphs.
- [ ] Close MMFF94/UFF typing, charge, parameter, convergence, and stereo gaps
  with independent oracle and soundness gates. The six-molecule/32-atom
  pyridinium/macrocycle residual boundary is pinned in
  `validation/manifests/mmff94_issue337_pyridinium_sssr_residual.json`; do not
  replace this with a local atom-type heuristic.
- [x] Re-run the bounded independent RDKit MMFF94 availability oracle over the
  pinned 265-molecule A/B corpus. The current RDKit 2025.09.3 run parses,
  embeds, constructs MMFF94, and produces finite energies for 265/265; the
  minimize return-code split is retained as an explicit oracle outcome rather
  than being treated as convergence parity. See
  `validation/results/mmff94-rdkit-availability-oracle-v1.0.10.json`.
- [x] Measure the #337 D2-root hypothesis with an all-root diagnostic; the
  six residual molecules have identical D2/all-root candidate sets, so root
  enumeration is not treated as the fix.
- [x] Probe equal-length macrocycle edge exchanges for #337; all 6 residuals
  yield same-size GF(2)-exchangeable alternatives. The candidate population
  is intentionally not promoted until a permutation-invariant representative
  tie-break is validated.
- [x] Specify the bounded #337 relevant-cycle selector boundary and fail-closed
  acceptance gates in `docs/rfcs/mmff94_relevant_cycle_selector.md`.
- [x] Expose the bounded #337 selector's cap outcome as a typed diagnostic;
  `find_symmetrized_sssr_with_diagnostics()` returns the complete Horton basis
  with `CapExhausted` instead of exposing a partial candidate family.
- [x] Make the #337 candidate tie-break bond-order aware across eight WL
  rounds and collect D2 roots against one immutable bond graph; retain the
  RDKit representative-family parity gate as a separate, still-open completion
  item.
- [x] Record the #337 local determinism boundary as a machine-readable result;
  `validation/results/mmff94-issue337-determinism-v1.0.10.json` captures six
  fixtures and 1,536 seeded relabeling checks while keeping RDKit representative
  family parity explicitly open.
- [x] Measure deterministic ensemble diversity, class-level failure rates,
  symmetry-aware RMSD/TFD, and energy sanity. The local deterministic-diversity
  slice is now recorded for 58/58 same-seed reproductions, 58/58 adjacent-seed
  non-aliasing checks, and three flexible molecules across eight seeds;
  the status-class failure-rate slice is also recorded as 58/58 new-embedder
  successes across ten categories; the bounded force-field energy-sanity slice
  is recorded for four arms with 61/61 finite and non-increasing successful
  outcomes each. The pinned 250-molecule RDKit TFD measurement now validates
  201 finite paired values plus 49 explicitly rigid/not-applicable probes;
  `scripts/check_tfd_oracle_evidence.py` records that boundary. This closes
  the bounded measurement task; the broader P5 quality item remains open
  because these measurements are not an independent force-field oracle.
  A six-case automorphism-aware RMSD comparison against RDKit's `GetBestRMS`
  now has 5/5 non-gap agreements within 0.001 Å; the acetate
  `symmetrizeConjugatedTerminalGroups` difference remains explicitly disclosed.
- [x] Record the local deterministic-diversity slice as a reproducible
  v1.0.10 benchmark artifact; the checked-in report preserves the exact gate
  command, seed count, alignment method, and measured RMSD summary.
- [x] Add the bounded automorphism-aware torsion-distance diagnostic and local
  invariants for symmetric relabeling and changed-torsion detection. A six-
  molecule corpus gate now records six valid measurements and bounded,
  change-sensitive distances. This is a reusable local foundation; it does
  not replace the pinned RDKit TFD measurement below.
- [x] Validate the pinned 250-molecule RDKit TFD measurement boundary. The
  checked-in result contains 201 finite paired values and 49 explicitly
  rigid/not-applicable probes; `scripts/check_tfd_oracle_evidence.py` verifies
  unique corpus identity, finite values, and the rigid-probe boundary. This
  is measurement evidence, not independent force-field correctness or full
  conformer-quality parity.
- [x] Add a dependency-free validator for the bounded P5 quality-evidence
  bundle. `scripts/check_3d_quality_evidence.py` verifies the pinned diversity,
  failure-rate, energy, RMSD, and torsion records; it does not promote them to
  independent force-field or full TFD parity evidence.

## P6 — Ecosystem durability

- [x] Publish compatibility, provenance, licensing, security, benchmark, and
  migration documents that state unsupported scope explicitly.
- [x] Record the open-issue triage and distinguish implemented bounded APIs
  from unresolved correctness, safety, performance, and CI work.
- [x] Add stable extension points and a contributor corpus policy covering
  provenance, licensing, minimization, and oracle versioning.
- [x] Publish a reproducible compatibility dashboard only after a clean
  checkout can regenerate it.

## Security and release gates

## Dependency classification for autonomous work

The autonomous local objective covers the first two rows only when their
inputs are available in the current checkout. The last row is intentionally
never marked complete from local self-review or a local dry run.

| Class | Examples | Completion evidence |
|---|---|---|
| Internal implementation | parser limits, canonicalization invariants, shared contracts, local benchmarks, Rust/Python/Node/WASM tests, deterministic reports | checked-in code or fixture, focused test output, and reproducible local command |
| Local-environment dependent | browser-engine memory/timing, regenerated Web WASM artifacts, cross-platform measurements, optional installed engines | the required toolchain/artifact is present, then a captured run with platform metadata; otherwise remain open with the missing tool recorded |
| External coordination | independent S5 review, S6 advisory rehearsal with publication/backport, future registry or release publication | external artifact or reviewer evidence; local preparation is not completion |

The complete open-item disposition is maintained in
[`docs/roadmap-open-work.md`](docs/roadmap-open-work.md); it distinguishes
local implementation, environment/toolchain dependency, explicitly abandoned
stretch goals, and external coordination.

- [x] S0/S1: threat model, bounded public input, typed resource-limit errors,
  and parser-wide adversarial coverage.
- [x] S2: unsafe-surface inventory, optional native-InChI FFI isolation, four
  fuzz targets, focused Miri, and ASan/LSan/TSan procedures and evidence.
- [x] S3: shared binding contracts plus Python, Node/WASM, browser, MCP, CLI,
  and filesystem boundary tests.
- [x] S4: dependency/license checks, immutable Actions, SBOM/provenance,
  checksums, release-key custody, and registry verification.
- [x] S5 local gate: the repository-local review packet and executable review
  checks are complete.
- [ ] S5 independent gate: obtain a non-maintainer review or external audit of
  major parser, serialization, and binding boundaries. This is external work
  and must not be marked complete from self-review.
- [ ] S6 continuous maintenance: rehearse advisory intake, fix, backport,
  artifact publication, and supported-version synchronization.

Reproduction instructions are in
[`docs/v1.0-local-release-gate.md`](docs/v1.0-local-release-gate.md),
[`docs/security-review/`](docs/security-review/), and
[`SECURITY.md`](SECURITY.md).

## Compatibility boundaries

The v1.0 contract remains unchanged:

- `canonical_smiles()` is a representation; the stable-key API is fail-closed.
- Aromaticity and CIP expose explicit default/opt-in models and do not promise
  universal RDKit parity.
- Python `RWMol`, CDXML editing, and Markush/polymer expansion are bounded
  subsets.
- 3D/MMFF94 are Experimental.
- Pure-Rust InChI is approximate; exact standard InChI is opt-in native FFI.

See [`docs/compatibility-scope.md`](docs/compatibility-scope.md).

## Execution order

1. Finish the P0 scorecard validator.
2. Extend equivalent format-streaming benchmarks and input-safety fixtures.
3. Close canonical E/Z invariance before broadening identity-dependent APIs.
4. Unify stable cross-binding fixtures and expected results.
5. Stabilize fingerprint/descriptor contracts and held-out reports.
6. Advance reaction, 3D/materials, and ecosystem breadth only behind their
   respective correctness and measurement gates.

Do not start broad feature expansion while a shared primitive has a known
silent-corruption or unbounded-work regression.

## Definition of “ahead”

chematic is ahead of a competitor only for a named use case and a measured
dimension. Every claim must identify versions, corpus, configuration, hardware,
failure policy, and reproduction command. “Faster”, “more accurate”, “smaller”,
and “more compatible” are otherwise hypotheses, not product claims.
