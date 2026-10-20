# Open Babel interchange and semantic preservation plan

Status: proposed execution plan, 2026-10-11. No superiority gate is complete.

The objective is to exceed Open Babel in useful read/write coverage and in
the chemical and document information retained through conversion. Establish
measured preservation on common workflows first, then expand coverage using
the same acceptance contract for each new adapter.

## October 11 implementation progress

The first IO1 slice is implemented in the local source candidate:

- Partial charges no longer become rounded formal charges; `am` reads as
  single, while unsupported bond tokens refuse. Graph-only output declares
  `NO_CHARGES` instead of claiming computed Gasteiger charges.
- Additive immutable `Mol2SourceRecord` APIs retain source text and expose source
  atom/bond IDs, names, types, substructure references, partial charges,
  coordinates, and bond tokens. Extra sections are explicitly opaque. General
  formal-charge/UNITY and geometry-stereo interpretation remain unsupported.
- CLI/Python/WASM MOL2-to-MOL2 uses validated exact-source retention. Python
  and WASM expose the same `chematic.mol2-record.v1` inspection schema. The
  graph-only APIs retain their signatures; malformed counts, duplicate
  references, nonfinite numbers, and multiple records now reject.
- Native/CLI tests passed 803 cases; the source-built CPython 3.13 wheel
  passed 179 cases. Fresh Node/WASM passed the new preservation suite and
  existing cross-binding and format-parity suites. Clippy passed for MOL/CLI.

Evidence and artifact/source hashes are recorded in
[`mol2-preservation-source-2026-10-11.json`](https://github.com/kent-tokyo/chematic/blob/main/validation/results/mol2-preservation-source-2026-10-11.json).
This is development evidence for source retention and the bounded reader
repairs. IO0 comparison, the rest of IO1, general checked conversion in IO2,
released-package acceptance, and Open Babel superiority remain open. Exact
source retention does not establish interpreted/editable opaque semantics.

The initial local IO2 implementation adds `convert_with_report` with one Rust policy and
report implementation exposed by CLI/Python/WASM. Initial routes are MOL2 →
MOL2 (exact bytes, opaque semantics explicitly separate) and MOL2 → SMILES
(field-group losses, strict refusal, opt-in lossy candidate). Candidate graph
verification reparses SMILES and compares mapped atoms/bonds using the new
SMILES writer correspondence API. Formal-charge interpretation, source
hydrogen assignment and stereo are not claimed verified. See the
[checked conversion contract](format-capabilities.md#checked-conversion-reports-initial-adapter)
for status, refusal, mapping and resource-limit semantics. Other IO2 routes
and the paired IO0 comparison remain open. MOL/CLI passed 816 tests; unchanged
SMILES sources passed 316 tests (two unrelated doctests ignored). The fresh
source wheel passed 191 tests, Node/WASM passed four suites, and all 17 shared
CLI/Python/WASM report/error cases agreed. MOL/SMILES/CLI Clippy passed.
Evidence is recorded in
[`conversion-report-source-2026-10-11.json`](https://github.com/kent-tokyo/chematic/blob/main/validation/results/conversion-report-source-2026-10-11.json).

The coordinate-preservation slice extends that adapter to MOL2 → CJSON.
Source x/y/z values and projected formal charges are emitted, then reparsed
and checked through source atom correspondence. Partial charges, source
names/types and opaque sections are still explicit losses; strict policy
still refuses. Signed zero and extreme finite f64 values are included in
precision regressions. The new tests reproduced a pre-existing CJSON numeric
reader error (`-1.4175225513503917` changed by 1 ULP), repaired by enabling
serde_json's precise float parsing feature. The existing CJSON reader benefits
from the same fix. At this milestone, CJSON input, broader MOL2 charge/stereo
interpretation and the paired IO0 comparison remained open. MOL/CLI passed 823 tests, the fresh
source wheel passed 201, Node/WASM passed four suites, 30 shared report/error
cases agreed, and Clippy passed. The precision test covers 2,045 finite f64
vectors from a fixed seed, each with signed-zero and opposite-sign components.
See [`conversion-cjson-source-2026-10-11.json`](https://github.com/kent-tokyo/chematic/blob/main/validation/results/conversion-cjson-source-2026-10-11.json).

The CJSON-input slice adds immutable source records, inspection of opaque
field paths, strict validation of known arrays and bond indices/orders, and
recursive rejection of duplicate JSON keys. Existing permissive CJSON readers
remain compatible. Checked CJSON → CJSON retains exact UTF-8 source;
CJSON → SMILES reports field losses and defaults, verifies candidate element,
formal-charge, isotope and bond correspondence, and refuses under strict policy.
Source IDs in CJSON mappings are zero-based array indices. Atom aromaticity,
coordinate stereo and opaque field semantics remain outside interpretation.
This is a local version-1/absent-version profile, not support for every CJSON
variant. Broader adapters and the paired IO0 comparison remain open.
Local verification passed 835 MOL/CLI tests, 211 source-wheel Python tests,
five Node/WASM suites and 53 shared report/error cases (36 identical reports,
17 shared errors). Both CJSON readers preserve all coordinate bits in the
2,045-vector seeded test. MOL/CLI Clippy passed with `--no-deps`; dependency
linting encountered `needless_range_loop` in concurrently modified FF code,
which this interchange slice did not alter.
See [the CJSON record contract](format-capabilities.md#immutable-chemicaljson-records)
and [`cjson-input-source-2026-10-11.json`](https://github.com/kent-tokyo/chematic/blob/main/validation/results/cjson-input-source-2026-10-11.json)
for source-candidate verification evidence.

The first paired IO0 development runner now freezes eight exposed inputs,
two routes and both artifact identities, checks explicit field expectations,
and retains all 32 terminal outcomes. Common, stress and negative cases remain
separate; mapping ambiguity, limits and refusals cannot count as preservation.
It exposed missing partial-charge/name metadata in checked MOL2 → CJSON,
which the source candidate now retains and verifies. No missing charge is
filled, no method is computed, and formal charges stay separate.

Targeted cross reads reproduced Open Babel 3.2.1 reading floating integer bond
orders as zero; the CJSON writer now emits integer orders as JSON integers.
Open Babel self-read controls retain evidence of remaining name/partial-charge
losses on its default CJSON → MOL2 path. See the
[paired protocol, commands and records](https://github.com/kent-tokyo/chematic/blob/main/validation/interchange/README.md).
This is exposed development evidence; independent acceptance and superiority
remain open.

## Inspected baseline

This audit covers local source at `7f84ae4a98399086ce9f389c099f3fb77b66066e`
(workspace manifest `1.0.20`) with existing documentation/Python-example
changes. It is not an audit of later tags or published packages. The local
`/opt/homebrew/bin/obabel -V` reports Open Babel 3.2.1; its `-L formats`
inventory was inspected. These observations identify candidates for testing,
not a completed comparative measurement.

Open Babel has a broad [format catalog](https://openbabel.org/docs/FileFormats/Overview.html),
including computational chemistry, crystallography, reactions, docking, and
drawing formats. Count format families, dialects, read/write directions, and
workflow support separately; aliases such as `mol`, `mdl`, `sd`, and `sdf`
must not inflate coverage. Archive the actual pinned binary's inventory when
starting the comparison, since the online manual can describe another build.

| Inspected surface | Existing foundation | Gap to address |
|---|---|---|
| [MOL2 reader/writer](https://github.com/kent-tokyo/chematic/blob/main/crates/chematic-mol/src/mol2_tripos.rs) | Graph and coordinate I/O with parse limits | Reader rounds partial charge into formal charge, maps `am` to aromatic, and substitutes single for unknown bond tokens. Writer emits zero charges, simplified atom types, and a `GASTEIGER` label without computing charges. Names/substructures are not retained by the return model. |
| [CLI conversion](https://github.com/kent-tokyo/chematic/blob/main/crates/chematic-cli/src/main.rs) | Common-format bridge; special V3000 self-conversion retains metadata/coordinates and certain unsupported query conversions refuse | Other paths project into a molecule and reconstruct output with default metadata/empty coordinates. There is no common per-field conversion report. |
| [Python conversion](https://github.com/kent-tokyo/chematic/blob/main/crates/chematic-py/python/chematic/__init__.py) | Common graph/coordinate bridge | Carries selected input coordinates; format metadata is not promised, and missing PDBQT charges become zeros. MOL2 input coordinates are not carried by this bridge. |
| [WASM conversion](https://github.com/kent-tokyo/chematic/blob/main/crates/chematic-wasm/src/mol_io.rs) | Explicit topology-only bridge | Drops coordinates and format metadata. It does not share the CLI's V3000 self-conversion/query guard implementation. |
| [V3000 fixtures](https://github.com/kent-tokyo/chematic/blob/main/crates/chematic-mol/tests/v3000_semantic_interop.rs) | SGROUP/COLLECTION ordering, references, and retained attributes | Opaque retention needs its own status; it does not establish editable or query chemistry. |
| [CDXML document](https://github.com/kent-tokyo/chematic/blob/main/crates/chematic-mol/src/cdxml_document.rs), [semantic model](https://github.com/kent-tokyo/chematic/blob/main/crates/chematic-mol/src/semantic.rs), [reaction adapter](https://github.com/kent-tokyo/chematic/blob/main/crates/chematic-mol/src/rxn.rs) | Document retention, typed constructs, and explicit reaction conversion losses | Reuse these contracts in conversion instead of flattening every object to a molecule. |
| [PDBQT](https://github.com/kent-tokyo/chematic/blob/main/crates/chematic-mol/src/pdbqt.rs) | Rigid-body I/O | Flexible ligand torsion trees need a typed representation and reader/writer tests. Open Babel documents torsion-tree and rigid output options in its [PDBQT manual](https://openbabel.org/docs/FileFormats/AutoDock_PDBQT_format.html). |
| [Stereo gate](https://github.com/kent-tokyo/chematic/blob/main/scripts/stereo_format_roundtrip_gate.py), [binding fixtures](https://github.com/kent-tokyo/chematic/blob/main/crates/chematic-mol/tests/format_binding_parity.rs), [CLI timing](https://github.com/kent-tokyo/chematic/blob/main/scripts/benchmarks/bench_openbabel_cli.py) | Correspondence-aware stereo diagnostics, known-value fixtures, and separate CLI timings | Extend into a paired Open Babel preservation gate; timing alone is not semantic validation. |

## Preservation contract

Define success per input record, requested field, format pair, and policy.
Writing parseable output or matching a formula is insufficient.

| Information | Required comparison |
|---|---|
| Atom/bond identity | Element, isotope, formal charge, radical state, explicit/implicit hydrogen policy, atom maps, bond order, and aromaticity under a declared model; compare through atom/bond correspondence. |
| Stereochemistry | Tetrahedral, E/Z, enhanced groups, specified/unspecified state, and supported non-tetrahedral configurations. Compare mapped configurations and query truth tables where applicable. A changed canonical spelling alone is not loss. |
| Coordinates | Dimension, units, atom correspondence, conformer/model count, and predeclared format-specific precision. Missing geometry must stay missing unless generation was requested. |
| Charges and typing | Partial charge values and declared method, separately from formal charge; source atom types, names, residue/substructure ownership, and any explicit retyping provenance. |
| Records and documents | Record order/count, names, ordered/repeated properties, reaction roles/maps, SGROUP/attachment references, and supported document objects. |
| Domain objects | Lattice/periodicity/occupancy, trajectory frames, volumetric axes/units/values, and quantum metadata in their own typed containers. |

Each present source feature must be accounted for as `preserved`,
`normalized`, `inferred`, `opaque_retained`, `dropped`, `unsupported`, or
`unverified`. The initial report implements the states needed by its MOL2 adapter; broader adapter guarantees remain proposed.
Only an independently checked field can be `preserved`. A report must retain
source record/object identity, parse/write/verification stage, reason, and
the atom/bond mapping or explicit correspondence failure.

For example, MOL2 partial charge `-0.70` must not become formal charge `-1`.
Reading `am` must not introduce aromaticity. SDF to plain SMILES must disclose
lost properties/geometry; XYZ to a bonded molecule must disclose any bond
inference. Generated charges/coordinates must carry their method, not be
presented as source facts. Preserved opaque bytes must never be labeled as
interpreted polymer, Markush, or query semantics.

Offer an additive checked API, provisionally `convert_with_report`, with:

- `strict`: refuse when requested information cannot be represented or
  verified. Return the reasons without publishing a successful output file.
- `allow_loss`: return output and all detected losses under an explicit policy.
- `preserve_source`: retain supported source objects for same-format editing;
  report opaque retention separately from interpreted fields.

Keep existing 1.x signatures compatible. Put the conversion policy and
adapters in Rust and expose the same report schema to CLI/Python/Node/WASM.
Do not add required fields to existing public structs merely to build the
new envelope. Reuse molecule, reaction, crystal, grid, and document models;
an ordinary `Molecule` cannot be the universal interchange container.

Loss accounting starts at the source reader, before projection discards
fields. An empty report from an adapter that never inspected unknown source
constructs is not proof of lossless conversion. Unsupported constructs must
be retained opaquely where safe or classified as unsupported/unverified.
Opaque references must be validated or invalidated when edits renumber or
delete their target atoms/objects.

## Execution order and acceptance

The proposed order below applies to interchange work alongside the existing
A2/A4 correctness program. Each row remains open until its evidence exists.

| Step | Deliverable | Acceptance |
|---|---|---|
| IO0 | Freeze a baseline and build the paired preservation harness | Record both implementations/artifact hashes, input hashes, options, supported fields, oracle/hand-derived expectations, resource budgets, and one terminal result for every input. Separate source-built and public-package runs. |
| IO1 | Repair MOL2 charge/bond semantics and add a metadata-bearing record API | Independent fixtures cover negative/positive fractional charges, malformed/nonfinite charges, `am`/`ar`, supported/unknown bond tokens, atom types, names/substructures, and precision. Preserve partial charges separately; derive formal charges only under a declared chemistry policy or refuse ambiguity. Do not label zeros as computed Gasteiger charges. Every new repair has a minimized regression. |
| IO2 | Introduce checked conversion and shared loss accounting | First cover SMILES/CXSMILES, MOL V2000/V3000, SDF, MOL2, CML, and CDXML. Same-format conversion retains representable metadata/geometry. Cross-format paths preserve requested fields or classify/refuse each loss. Include multi-record input, malformed records, and V3000 query guards across bindings. |
| IO3 | Complete major workflow gaps | PDB/mmCIF/PQR/PDBQT retain supported charge/type/residue/model data; implement flexible PDBQT torsion trees with validation. RXN/reaction SMILES and CDXML/MRV/KET use document/semantic adapters. Each unsupported dialect has an explicit boundary. |
| IO4 | Expand read/write coverage by workflow | Inventory gaps against the pinned Open Babel catalog. Prioritize Gaussian output/fchk, GAMESS, ORCA, Q-Chem, GRO, and other requested scientific pipelines after confirming source specifications, actual local support, and user demand. Use crystal/grid/trajectory/QC types as appropriate. Add format families only with fixtures, reports, limits, and end-to-end tests. |
| IO5 | Publish package-level comparisons | Clean-install released Rust/Python/npm artifacts and exercise supported workflows. Publish coverage, preservation, loss detection, refusal, errors, timeouts, throughput, p95, and RSS with full denominators. Retain wins and losses for both engines. |

Begin IO0/IO1 with a small hand-reviewed suite that reproduces the inspected
MOL2 issues and existing V3000 retention boundaries. Build development
coverage from minimized cases before acquiring an independent acceptance
corpus. Exposed fixtures stay development evidence after fixes; they cannot
be relabeled as unused holdout inputs.

## Comparison and claim rules

For each format pair, run self round trips and cross-engine directions:
Open Babel output read by chematic and chematic output read by Open Babel.
Use independent field expectations, hand-reviewed cases, and a pinned third
reader where useful. Neither Open Babel nor chematic is the sole oracle.
Check source fixtures first; source-generation failures remain visible in
the full accounting and cannot become unexplained exclusions.

Test 0D/2D/3D separately. Open Babel's
[MOL manual](https://openbabel.org/docs/FileFormats/MDL_MOL_format.html)
documents parity/coordinate options and a 0D stereo extension; record the
exact flags and distinguish standard and extension lanes. Also test atom
renumbering, equivalent spellings, property order/duplication, multiple
records, unsupported constructs, malformed input, timeout, and limits.

Report these outcomes separately: verified preservation, detected loss with
output, silent loss/wrong meaning, safe refusal, parse/write failure,
timeout/resource limit, and unresolved verification. A refusal protects a
caller but contributes zero successfully preserved coverage. Keep declared
support fixed before the run so that exclusions cannot improve the score.

A first bounded claim requires zero silent losses on the frozen supported
suite and complete requested-field preservation on its valid inputs.
Superiority additionally requires better preservation on a predeclared
paired independent corpus, with a positive lower bound on the paired 95%
accuracy-difference interval and no regression in successfully preserved
coverage on the required workflow cells. Report rare semantic categories
separately so aggregate counts cannot hide stereo/query failures. A small
fixture win is a narrow regression result, not general superiority.

Track breadth separately: compare deduplicated format families and verified
read/write directions, then publish a workflow-weighted score with weights
frozen in advance. Until those comparisons pass, claim only the specific
format pairs/fields that have evidence. Broader Open Babel replacement
requires measured coverage beyond the initial molecular subset.

Measure speed only after correctness, comparing equivalent CLI workloads
with CLI and equivalent API workloads with API. Freeze hardware, threads,
options, startup treatment, and output requirements; preserve failures and
resource use alongside timings. Safe parsing and typed errors need malformed
input tests and measured bounds; Rust alone is not evidence of superiority.

The reverse-adapter slice adds typed root CJSON partial-charge datasets and
`properties.name`, then checked CJSON → MOL2. It requires explicit coordinates
and one representable dataset at most; it refuses ambiguous dataset selection,
missing geometry and unsafe header text. Candidate float values are written
without fixed-precision rounding and checked bit for bit after reparse. Generated
MOL2 atom names/types/IDs, unsupported formal charges/isotopes, opaque fields and
stereo remain visible losses or uncertainty. Existing graph-only writers remain
unchanged. See [the adapter contract](format-capabilities.md#immutable-chemicaljson-records).

The expanded cross-read lane adds six CheMatic MOL2 projections from both
engines' CJSON on the three exposed common inputs. It retains 15 terminal rows
and separately checks field preservation against the existing declared oracle.
All six new projections preserve name, charge method and partial charges;
atom names/types remain losses. The nine historical control outcomes and
field checks are unchanged. Local gates pass 848 MOL/CLI tests, 230 isolated
source-wheel Python tests, five Node/WASM suites, 79 shared report/error cases
and 33 runner tests. See [the indexed evidence](https://github.com/kent-tokyo/chematic/blob/main/validation/results/cjson-mol2-development-summary-2026-10-11.json).


The original v1.1.0 integration preserved the public editable MOL2 record and
stereo-aware SMILES writer. The immutable profile uses `Mol2SourceRecord` /
`parse_mol2_source_record` names to avoid changing that contract. The preceding
v1.0.20 measurements and v1.1.0 integration checks are historical; they are not
current-package results.

The original integration validation is recorded in
`validation/results/main-integration-2026-10-11/summary.json`.
The earlier v1.0.20 development measurements remain historical evidence.
