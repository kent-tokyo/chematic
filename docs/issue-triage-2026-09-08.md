# Open issue triage — 2026-09-08

This snapshot records the bounded completion of issue #507. It does not turn
format-specific safety cases into a throughput, compatibility, or universal
parser claim.

## P2 descriptor binding parity — bounded completion

The 5,000-row descriptor contract was replayed from the committed SMILES
corpus across the Rust crate example, the PyO3 extension, and the Node/WASM
package. All three bindings parsed 5,000/5,000 rows and matched exactly for
MW, TPSA, HBD, HBA, and heavy-atom count (floating-point tolerance 1e-9).
The machine-readable evidence is
`validation/results/descriptor-cross-binding-parity-5000-v1.0.9.json`.
This closes only the same-source binding agreement slice; the broader held-out
fingerprint, topology, torsion, and standardization parity reports remain open,
as does external RDKit accuracy beyond their separately measured reports.

## P2 native ECFP4 binding parity — bounded completion

The native 2,048-bit ECFP4 representation was replayed on the same 5,000-row
corpus through the Rust crate example, PyO3 extension, and Node/WASM package.
All three bindings produced 5,000/5,000 successful records and identical
LSB-first 256-byte values. Evidence is recorded in
`validation/results/ecfp4-cross-binding-parity-5000-v1.0.9.json`.
This does not close RDKit-compatible ECFP4 accuracy, sparse/count provenance,
or the remaining held-out fingerprint operations.

## P2 native MACCS binding parity — bounded completion

The native 166-key MACCS representation was replayed over the same 5,000-row
corpus through Rust, Python, and Node/WASM. The three bindings produced
5,000/5,000 successful records and identical 21-byte LSB-first values. The
machine-readable evidence is
`validation/results/maccs-cross-binding-parity-5000-v1.0.9.json`.
The external RDKit key mapping and remaining fingerprint parity gates stay
separate and open.

## P2 native topological-path binding parity — bounded completion

The native `topo_path` 2,048-bit representation was replayed over the same
5,000-row corpus through Rust, Python, and Node/WASM. All bindings produced
5,000/5,000 successful records with identical LSB-first bytes. Evidence is
recorded in `validation/results/topo-path-cross-binding-parity-5000-v1.0.9.json`.
RDKit-compatible path accuracy and the remaining torsion/standardization gates
remain separate.

## P2 native torsion binding parity — bounded completion

The native 2,048-bit torsion representation was replayed over the same
5,000-row corpus through Rust, Python, and Node/WASM. All bindings produced
5,000/5,000 successful records with identical LSB-first bytes. Evidence is
recorded in `validation/results/torsion-cross-binding-parity-5000-v1.0.9.json`.
RDKit-compatible torsion accuracy remains a separate gate.

## P2 standardization binding parity — bounded completion

The explicit shared profile (largest fragment, charge neutralization, explicit
hydrogen removal, and canonical tautomer) was replayed over 5,000 corpus rows
through Rust, Python, and Node/WASM. All bindings produced 5,000/5,000
successful records with identical canonical outputs. Evidence is recorded in
`validation/results/standardization-cross-binding-parity-5000-v1.0.9.json`.
The ten Phase-1 expected-identity holdouts and external standardization quality
remain separate gates.

## P2 RDKit-compatible ECFP4 binding parity — bounded completion

The RDKit-compatible ECFP4 implementation was replayed on 5,000 shared
SMILES through Rust, Python, and Node/WASM. The report compares exact
2,048-bit values for successful preprocessing and compares failures by typed
status rather than unstable error text. Evidence is recorded in
`validation/results/rdkit-ecfp4-cross-binding-parity-5000-v1.0.9.json`.
Independent RDKit accuracy and sparse/count/bitInfo parity remain separate
gates.

## Issues #149 and #503 — aromatic E/Z stash residual boundary

The current joint carrier resolver and canonical-fidelity partition remain
green on the held-out residual suite, but the three measured aromatic-stash
families still produce two valid canonical spellings and therefore remain
fail-closed by `canonical_smiles_stable_key()`. The coupled-carrier audit was
re-run against the committed 5,000-line corpus: 31 coupling components were
found, all of size 2, with no cycle or zero-private-substituent component.
This confirms the residual is not an unbounded coupling-graph case. A safe
fix still needs a representation-independent aromatic carrier/stereo traversal
rule; no index-based winner was introduced, and both issues remain open.

The follow-up experiment also treated literal `/` and `\\` carrier edges
adjacent to a double bond as non-discriminating during orbit coloring, while
pinning the structural alkene endpoint. The focused residual suite (6/6) and
the `chematic-smiles` library suite (221/221) stayed green, but the three
held-out families still retained their two output variants. The experiment was
therefore reverted: it did not provide a convergence rule, and changing the
literal-carrier semantics would expand the risk without satisfying #503.

Additional solver tracing distinguishes the remaining failure mode from a
rank tie: in two of the held-out spellings, every joint carrier assignment is
rejected because the available non-shared carriers are on the DFS ring
close-side (and the shared assignment conflicts or is likewise not writable).
Temporarily allowing close-side candidates did not make the three families
converge. The close-side safety rule is therefore retained; the next viable
fix must choose a representation-independent traversal/open-side arrangement
for shared carrier bonds, with a round-trip geometry gate.

A third experiment biased canonical DFS toward topology-derived shared carrier
bonds so that those bonds would more often remain on the spanning tree. The
smiles residual suite stayed green only after accepting a newly formatted
variant, while the held-out boundary remained fail-closed and the private
two-way residual probe did not converge. The bias was reverted because it
changed canonical spellings without establishing representation-independent
stereo resolution.

A fourth experiment removed the index-derived uniqueness pin from structural
alkene endpoints adjacent to aromatic direction stashes. It did not reduce
the three held-out residuals, and it regressed the existing ring-closure stash
round-trip: the canonical output reparsed with the opposite aromatic ring
direction (`N(/C)=c\\1...` versus `N(/C)=c/1...`). The endpoint pin is therefore
retained and the experiment was reverted; a future fix must provide an
intrinsic stereo color or a traversal proof rather than simply unpinning the
endpoint.

Evidence:

- `cargo test -p chematic-smiles --lib --offline` — 221 passed.
- `cargo run -p chematic-smiles --release --example
  ez_shared_carrier_component_audit -- scripts/descriptor_census_corpus.smi`
  — 5,000 parsed, 31 size-2 components, maximum size 2.
- `cargo test --workspace --offline` — all workspace unit, integration, and
  doctests passed; known ignored tests remain unchanged.

## Issue #70 — Criterion gate local runner boundary

The process-level block arithmetic, ABBA/BAAB metadata, routing fixtures, and
strict malformed-input checks pass with `bash scripts/test_criterion_gate.sh`.
The runner now treats a sandbox-denied macOS `sysctl vm.loadavg` read as
`loadavg: "unavailable"` instead of aborting the measurement block, while
preserving the metadata schema. This is local portability evidence only;
hosted +5%, +10%, and contamination calibration remains required before #70
can close.

The same local contract was re-run on 2026-09-08 after the schema-v2 artifact
and environment fields were checked in: both ABBA and BAAB fixtures emitted
`measurement_unit=criterion_process_point_estimate`, timestamps, execution
order, load average, CPU model, and `/proc/stat` steal-tick fields, and the
synthetic routing/incident fixtures returned the expected route or no-route
decisions. This confirms the local shell contract, not the hosted sensitivity
calibration or required-adjacent trust gate.

## Issue #462 — polymer repeat editing boundary

The typed semantic command surface now includes
`SetPolymerRepeatCount { unit_id, repeat_count }`. A polymer unit with
`repeat_count: null` is accepted as a valid unselected editing state, remains
non-expandable until selected, and can then be updated through the shared JSON
command used by Rust, Python, and WASM. The expansion mapping remains
deterministic: the three-repeat `[*]CC[*]` fixture maps six generated atoms to
the source unit.
JSON decoding now uses checked conversions for `u32` repeat counts and
platform-sized selected-alternative indices, rejecting out-of-range numeric
values instead of allowing integer truncation.
Polymer `repeat_endpoint_atoms` now applies the same checked `u32` conversion,
so oversized endpoint indices are rejected rather than truncated.

This advances only the edit-to-expansion boundary. Nested Markush choices,
polymer contraction, and the broader typed R-group/polymer/biomolecule API
surface remain open, so #462 is not closed.

Evidence:

- `cargo test -p chematic-mol --offline` — semantic and polymer regression
  tests passed.
- Python and Node/WASM cross-binding contract tests cover the same null-state
  JSON and explicit repeat-count command.

## Issue #460 — rich RXN document validation boundary

The existing typed reaction-document foundation now validates each component
as one molecule with the molecule parser directly. This keeps the component
boundary distinct from the three-section reaction parser and reports the
component ID in the typed parse error. A regression test rejects a reaction
payload placed where a single component SMILES is required.

The RXN V2000 adapter now also rejects files whose header count differs from
the actual `$MOL` block count. Missing and extra blocks are surfaced as a
typed parse error rather than being silently synthesized or discarded.
Malformed or negative reactant/product counts are likewise rejected instead
of being coerced to zero.
The block scanner accepts both LF and CRLF marker lines while only scanning
the section after the header count.
Typed condition and provenance records now reject missing identifying fields,
so JSON-authored metadata cannot enter the model in an ambiguous state.
Condition keys are also unique within each step, preventing ordered JSON
records from silently overwriting or ambiguously merging the same condition.
The public `ReactionDocument::from_json_str` constructor now combines serde
deserialization with model validation. Python and WASM RXN-document writers
use this checked constructor, so binding callers cannot bypass component,
coefficient, ID, or metadata validation by supplying JSON directly.
Components now also expose optional `atom_maps` identities containing the
map number and component-local atom index. Derived documents populate these
from their SMILES, while authored documents that provide the field must match
the serialized molecule exactly; this makes atom-map identity explicit without
breaking older JSON that omits the additive field.
The V2000 MOL atom-map field is now read and written for three-digit values,
and the RXN adapter returns a typed loss for larger map numbers instead of
producing a shifted fixed-width field.

This is validation hardening only. Full upstream-backed RXN dialect support,
including any format-specific metadata beyond the current loss-aware V2000
adapter, remains open.

Evidence:

- `cargo test -p chematic-rxn document --offline` — 9 passed.
- `cargo test -p chematic-mol rxn --offline` — 11 passed.
- `git diff --check` — passed.

## Issue #461 — CDXML document attribute safety boundary

The loss-preserving CDXML document adapter now applies its configured
`max_attribute_bytes` limit consistently to document, page, and object
attributes. Oversized presentation metadata is rejected with the existing
typed resource-limit error instead of being parsed without the advertised
bound.
Document edits also reject empty or non-element replacement/insertion
payloads before mutating the source representation.
The parser now rejects unmatched closing pages and nested page elements with
typed invalid-document errors, keeping the ordered page model aligned with
the source structure rather than silently dropping or overwriting a page.
Attribute-edit commands now validate attribute names before serialization, so
quotes, whitespace, and markup characters cannot create malformed XML.
Edits preserve the source's CRLF/LF convention and whether the source ended
with a newline, avoiding unrelated representation churn.
Documents parsed with explicit resource limits retain those limits for later
edits and revalidation, so an edit cannot silently fall back to defaults.
If a document contains duplicate page IDs, edits now reject that ambiguous
target instead of silently modifying the first matching page.
The structural scanner also recognizes minified one-line CDXML while keeping
the original source representation byte-for-byte available through `write()`.
The edit path uses the same logical tag records for minified input, preserving
its compact layout while applying page and object edits.
The parser now requires exactly one opened and closed `<CDXML>` root with an
exact tag-name boundary. Garbage input, `<CDXMLFoo>`, and an unterminated root
are rejected as invalid documents rather than being accepted as empty
documents; pages outside the root are rejected as well.

This is a resource-safety slice only. Full ChemDraw presentation semantics,
including every style, geometry, grouping, and annotation dialect, remain
outside the current adapter contract.

Evidence:

- `cargo test -p chematic-mol cdxml_document --offline` — 16 passed.
- `git diff --check` — passed.

## Issue #337 — typed symmetrized-ring cap outcome

The bounded symmetrized-SSSR path now exposes
`find_symmetrized_sssr_with_diagnostics()`. It returns the selected ring set
with `SymmetrizedSssrStatus::Complete` or `CapExhausted`; cap exhaustion keeps
the complete Horton basis and never exposes a partial candidate family. The
existing `find_symmetrized_sssr()` API remains a compatibility wrapper.

This completes the fail-closed diagnostic boundary from
`docs/rfcs/mmff94_relevant_cycle_selector.md`, but does not close #337: the
permutation-invariant relevant-cycle representative policy and the six-fixture
MMFF94/RDKit parity gate remain open.

The upstream boundary is now explicit: current RDKit's default
`GetSymmSSSR` selects `atomRelevantCycles()` from its RDL ring-family
calculation; its legacy Figueras replacement path is separate. The current
schematic selector is a bounded, deterministic legacy-style expansion, so its
`Complete` diagnostic means local candidate enumeration completed, not that
the selected representative family is RDKit-equivalent. A production fix
still requires an independent relevant-cycle implementation or a validated
equivalence proof.

The candidate tie-break refinement is now bond-order aware and runs for eight
WL rounds, reducing accidental ties between chemically distinct paths while
retaining the existing permutation-stability contract. This is a deterministic
selection improvement, not an RDKit parity claim; the observed six-fixture
macrocycle counts remain the pinned boundary until the representative policy
is independently validated.

The latest candidate audit confirms the remaining mismatch is representative
selection rather than failure to enumerate a cycle family. In particular,
`chembl_tier_b_0028` returns four 32-member macrocycles from the bounded
selector, while the pinned RDKit result returns three. The extra candidate is
not removed by the existing GF(2) independence check because it can replace a
same-size basis member. A future selector therefore needs the RDL
representative-family rule, not an arbitrary lexicographic deletion or a
reduced independence test. Current RDKit's `FindRings.cpp` obtains SSSR cycles
through `RDL_getSSSR`; its legacy Figueras symmetrization path is separate.

Evidence:

- `cargo test -p chematic-perception --lib --offline` — 204 passed, 1 ignored.
- `cargo test -p chematic-ff --lib --offline` — 202 passed.

A follow-up removed raw `BondIdx` from the ordering of direct replacement
candidates and kept it only for exact edge-set identity and GF(2) rank
calculation. The six-fixture boundary and canonical macrocycle family remain
stable across 64 seeded atom relabelings per fixture. An attempted replacement of edge
identity itself with canonical rank keys was rejected immediately because it
collapsed distinct symmetric cycles (`chembl_tier_b_0023` changed from four
to two representatives); that change was reverted. The remaining identity
versus representative-orbit distinction is therefore explicit, but the RDKit
representative-family parity gate is still open.

## Issue #227 — MMFF94 coverage audit boundary

The checked-in 265-molecule audit separates classification errors from final
parameter resolution and must not be summarized as a single coverage number.
The current production stretch-bend path has no final unresolved rows in the
audit, while its type-only diagnostic still records 427 routing candidates
and 1,680 genuine table-gap rows. Bond, angle, and torsion gaps remain
separate deferred axes; the dominant aromatic typing residual is coupled to
the #337 symmetrized-ring/aromaticity boundary.

The audit therefore remains evidence for follow-up work, not a closure of
#227. In particular, a context-blind numeric-type substitution is explicitly
rejected because the checked-in negative simulation regresses furan. Any next
typing change requires a coordinated C/N/O/S oracle-parity gate.

The current audit also caught and fixed an independent P-typing collision:
the generic phosphorus path returned numeric type 20, whose registry entry is
the carbon-only CR4R type. The corrected path returns registry type 26 for
tricoordinate P and type 75 for P=C; the constructed phosphonium-ylide probe
now passes the semantic-compatibility invariant. The re-run removes the
typing error and reduces the bond/angle gate-would-fail count from 2 to 1,
but leaves one final unresolved angle and 24 torsion misses, so #227 remains
open.

Evidence: `validation/results/mmff94_coverage_227_term_audit_summary.json`,
`validation/results/mmff94_coverage_227_root_cause_classification.json`, and
the provenance decision in `scripts/mmff94_provenance/PROVENANCE.md`.

The workspace verification on this checkout also passes after the P-typing
change. A fresh Tier B audit still reports one final unresolved angle,
`(angle_type=0, type_i=43, type_j=18, type_k=63)`. It is intentionally
fail-closed: type 63 is absent from the checked-in eqLevel definition for the
required substitution path, so adding an inferred fallback would be an
unvalidated energy-model change rather than coverage work. The audit summary
is therefore pinned as `total=265`, `bond+angle-gate-would-fail=1`,
`bonds_final_unresolved=0`, `angles_final_unresolved=1`,
`torsions_missing=24`, and `stbn_final_unresolved=0`.

## Issue #303 — bounded structural slice completed

The explainable reactivity API now retains the existing epoxide, aziridine,
and Michael-acceptor findings and adds a deterministic bifunctional-electrophile
heuristic. It pairs deduplicated, non-overlapping motif sites, reports the
minimum heavy-atom bond distance, and caps output at 64 site pairs. The result
is explicitly structural triage only: it is not a genotoxicity classifier,
biological prediction, 3D geometry estimate, or cross-engine validation.

`validation/genotox_structural_fixtures.json` records three source-referenced
PubChem structures (Mitomycin C, Aflatoxin B1 8,9-epoxide, and anti-BPDE) with
retrieved isomeric SMILES and expected structural motifs. The test verifies
that the checked-in structures parse and exercise the expected motif rules;
the manifest carries provenance rather than a biological ground-truth label.

## Issue #303 evidence

- `cargo test -p chematic-chem --lib genotox --offline` — 6 tests passed,
  including all three source-referenced structures and the bifunctional case.
- The broader genotoxicity predictor, additional structural categories,
  biological fixtures, and licensing review remain outside this bounded slice.

## Issues #255 and #256 — completed placement slice

`generate_coords` now routes through the connectivity-ordered engine. The
engine repairs fused-ring seam placement and ring-chain-ring bridges, while
retaining deterministic new-island direct-bond anchoring. The existing
33-molecule evaluation records raw soundness 33/33, deterministic output
33/33, and UFF-only success 33/33.

Evidence from the current checkout:

- `cargo test -p chematic-3d --lib generate_coords_ --offline` — 30 passed.
- `cargo test -p chematic-3d --lib --offline` — 589 passed, 10 ignored.
- The remaining experimental MMFF94/UFF force-field residuals are separate
  from this coordinate-placement slice.

## Issue #507 — completed

The checked-in `validation/streaming_format_safety_cases.json` corpus now has
ten malformed-input cases for each supported runner format: SDF, V2000 MOL,
XYZ, Extended XYZ, V3000, MOL2, CML, CDXML, mmCIF, and PDB. The dependency-free
gate checks all 100 cases and requires exactly one failed record for each case.

The same gate retains ten oversized-input rejections and 20 gzip controls: one
valid decompressed control and one post-decompression input-limit rejection for
every format. CML, CDXML, and PDB remain explicitly line-limit safety checks
because their current readers are intentionally lenient about unknown, empty,
or non-record input.

## Evidence

- `python3 scripts/check_streaming_format_limits.py` — 100 negative cases,
  10 oversized cases, and 20 gzip cases passed.
- The corpus schema and exact ten-case-per-format count are validated before
  any runner invocation.

This closes the bounded malformed-corpus expansion in #507 and the local
Extended XYZ benchmark/safety extension. Cross-language
streaming parity, equivalent cross-engine throughput, and broader parser
semantics remain separate open roadmap gates.

## P1 cross-language streaming parity — completed local Python slice

The Python binding now exposes `iter_xyz_batched` and `iter_extxyz_batched`,
matching the existing SDF file-backed batch surface: bounded batch sizes,
lazy input-order emission, batch-boundary cancellation, and a deterministic
`manifest_json()` status envelope. Plain XYZ and Extended XYZ both return the
existing frame-dict shape used by `from_extxyz`, including coordinates,
lattice, per-atom properties, and frame metadata. Malformed frames are
counted as `rejected_frames` and do not escape as partial records.

Evidence from the current checkout:

- `cargo check -p chematic-py --offline` — passed.
- `cargo test -p chematic-py --offline --no-run` — passed.
- Python regression coverage added for SDF batch cancellation, XYZ ordering
  and batch boundaries, and Extended XYZ metadata plus cancellation.

This is a local Python parity slice only. WASM/browser streaming bindings,
shared cross-language fixture execution, and recovery after a malformed frame
remain open gates.

## P1 same-input Extended XYZ cross-engine contract — completed local slice

The existing cross-engine checker now accepts Extended XYZ and compares the
same checked-in two-frame fixture against RDKit's `MolFromXYZBlock` lane. The
contract requires equal valid-frame counts, zero failures, and identical
source-byte accounting while keeping the Rust file-backed `BufRead` boundary
separate from Python frame splitting.

Evidence:

- `python3 scripts/check_streaming_cross_engine.py --format extxyz --repeats 20`
  — chematic 40/0 and RDKit 40/0, both 5,360 source bytes.
- Report: `benchmarks/2026-09-08-streaming-cross-engine-extxyz.md`.

This does not close equivalent-operation throughput or same-process parity;
those remain open as required by the P1 completion contract.

## P1 gzip Extended XYZ cross-engine contract — completed local slice

The deterministic gzip contract now accepts Extended XYZ in addition to SDF
and XYZ. It feeds the same decompressed two-frame fixture to RDKit's XYZ block
parser while chematic consumes the gzip stream through its file-backed reader.
The report keeps compressed and decompressed byte stages explicit and checks
only record/failure accounting.

Evidence:

- `python3 scripts/check_streaming_gzip_contract.py --format extxyz --repeats 20`
  — 40/0 for chematic and RDKit.
- Report: `benchmarks/2026-09-08-streaming-gzip-extxyz-contract.md`.

This is a bounded gzip/Extended XYZ slice; equivalent-operation throughput,
same-process parity, and the broader strict malformed corpus remain open.

## P2 exact MACCS binding fixture — completed local slice

The shared binding manifest now carries exact 2048-bit ECFP4 bit positions and
166-bit MACCS bytes for the four small fingerprint fixtures. Rust, Python, and
Node/WASM compare the same least-significant-bit-first representation rather
than only checking length and non-emptiness.

Evidence:

- Rust `descriptor_contract` test, Python `test_cross_binding_contract.py`,
  and Node/WASM `cross_binding_contract.test.mjs` consume the same
  `ecfp4_bits` and `maccs_hex` values.
- This is a four-fixture cross-binding contract; it does not close the 5,000
  molecule held-out parity report or establish RDKit semantic parity.

The same manifest now also covers native `topo_path` bit positions. WASM has a
new `topo_path_bitvec` entry point so the native operation is directly
observable in all three binding test suites; this remains distinct from the
RDKit-compatible `rdkit_rdk` held-out lane.

The same exact-fixture contract now covers native topological torsion bits as
well. The native torsion operation is distinct from the RDKit-compatible
hashed torsion held-out report, whose 5,000-molecule parity remains open.

The shared fixture also now checks the RDKit-compatible hashed torsion bits;
WASM exposes this separately as `rdkit_torsion_bitvec`. This four-fixture
binding contract does not change the held-out 5,000-molecule result or its
known residual boundaries.

## P3 shared Extended XYZ fixture — completed local slice

`validation/cross_binding_contract.json` now contains a versioned
`extxyz_contract`. Rust `parse_extxyz`, Python `from_extxyz`, and Node/WASM
`extxyz_frame_json` consume the same input and check the same coordinates,
lattice, typed per-atom properties, and frame metadata.

Evidence:

- `cargo test -p chematic-mol --test cross_binding_adversarial --offline` — 2
  passed, including the shared Extended XYZ contract.
- `node crates/chematic-wasm/tests/cross_binding_contract.test.mjs` — passed.
- A clean offline maturin wheel smoke imported the current 1.0.9 binding and
  executed the shared RXN/Extended XYZ assertions successfully. Full pytest
  collection remains affected by the known pytest/pytest-asyncio environment
  mismatch.

This is one shared stable-operation fixture, not completion of the broader
all-stable-operation Rust/Python/Node/WASM manifest.

## P3 shared semantic expansion fixture — completed local slice

`validation/cross_binding_contract.json` now owns the Markush R-group and
polymer-repeat semantic cases, including the selected command result and the
complete source-to-expanded atom mapping. Rust `SemanticModel`, the Python
semantic functions, and the Node/WASM semantic functions consume the same two
cases.

Evidence:

- `cargo test -p chematic-mol --test cross_binding_adversarial --offline` — 3
  passed, including both semantic cases.
- `node crates/chematic-wasm/tests/semantic_expansion_contract.test.mjs` —
  passed after rebuilding the current-source Node/WASM package.
- The Python semantic assertions are wired to the same fixture; the focused
  pytest invocation remains affected by the known pytest/pytest-asyncio
  collection mismatch in the local environment.

This is a bounded shared-contract slice. It does not close the broader
all-stable-operation manifest item or the rich Markush/polymer/biomolecule
Issue #462.

## P2 standardization cross-binding fixture — completed local slice

The shared `standardization_contract` now contains all ten committed Phase 1
holdouts. Rust, Python, and Node/WASM apply the explicit largest-fragment
profile to the same inputs and compare the same returned SMILES. The external
held-out report remains separate: eight of ten expected fragment identities
match, while the two acetate cases intentionally disclose the charge-
neutralization difference.

Evidence:

- `cargo test -p chematic-chem --test descriptor_contract --offline` — 4
  passed.
- `node crates/chematic-wasm/tests/cross_binding_contract.test.mjs` — passed.
- Current-source Python extension direct check — 10/10 shared fixtures passed.

This closes only the ten-fixture cross-binding contract slice; the broader
5,000-row standardization parity and its external identity policy remain open.

## P1 same-input SDF cross-engine accounting — completed local slice

The installed Open Babel CLI was added to the checked-in 20-repetition SDF
record-accounting run. The identical 633-byte fixture produced 40 records and
zero failures in chematic, RDKit, and Open Babel. The report keeps chematic's
Rust `BufRead`, RDKit's Python supplier, and Open Babel's per-repetition CLI
startup as separate boundaries, so it does not claim same-condition speed.
The report also records the Open Babel executable version in machine-readable
`tool_versions` metadata.

Evidence: `benchmarks/2026-09-08-streaming-cross-engine-openbabel.md`.

The same runner now covers V3000 and MOL2 single-record fixtures, including
Open Babel's singular `molecule converted` output and pinned executable
version. The new reports are
`benchmarks/2026-09-08-streaming-cross-engine-openbabel-v3000.md` and
`benchmarks/2026-09-08-streaming-cross-engine-openbabel-mol2.md`.
The machine-readable runner now also records RDKit `2025.09.3` beside the
Open Babel version. The runner preserves the established implicit SDF Open
Babel comparison while keeping V3000/MOL2 comparisons explicit.

## Issue #460 — shared typed RXN document contract slice

The versioned cross-binding manifest now also owns the minimal authored
reaction-document fixture. Rust `ReactionDocument::from_json_str`, Python's
`to_rxn_document_json`/`from_rxn_document_json`, and Node/WASM's
`rxn_document_to_rxn`/`rxn_document_from_rxn` all consume the same document and
assert the same ordered reactant/product roles and SMILES after the V2000
round-trip.

Evidence:

- `cargo test -p chematic-rxn --test cross_binding_contract --offline` — 1
  passed.
- `node crates/chematic-wasm/tests/rxn_document_contract.test.mjs` — passed.
- A clean offline maturin wheel smoke executed the shared RXN assertions;
  full pytest collection remains affected by the known pytest/pytest-asyncio
  environment mismatch.

This strengthens the binding contract but does not close #460: full upstream
rich RXN dialect coverage, multi-step preservation, and broader loss fixtures
remain open.

## Issue #185 — UFF soundness observability

The existing fail-closed `sound` result is now accompanied by
`worst_bond_length` in `UffMinimizeResult`, Python `Mol.minimize_uff()`, and
WASM `minimize_uff_json()`. The value is computed from the same final geometry
used by the soundness gate, so callers do not need to duplicate the bond-length
calculation or mistake `converged` for geometrical validity.

This is diagnostic/safety-surface work only. UFF torsion and out-of-plane terms
remain unimplemented, so the fused-aromatic stationary-point residual and #185
itself remain open.

The `rejected_unsound_step` diagnostic now records only energy-decreasing
proposals rejected by the UFF geometry soundness gate. A caller-supplied
constraint may still reject a proposal, but that is no longer misreported as a
UFF unsoundness event; the distinction is covered by the constrained-minimizer
regression test.

## Issue #372 — canonical Boc/tBu symmetry performance

The exact twin/orbit path was re-run with the checked-in Tier A/B harness and
the canonical-search instrumentation feature. The run had zero old/new
correctness mismatches, zero search-budget exhaustions, and an 8.19x Tier A
geometric-mean speedup (Tier B negative control: 2.47x). Across Tier A the
exhaustive engine visited 6,186 leaves while the orbit-pruned engine wrote 13
leaves, visited 62 nodes, and performed 46 orbit tests; the repeated multi-Boc
and multi-pivaloyl fixtures each collapsed to one leaf from 432 exhaustive
leaves.

This is local proxy evidence only. The exact RENKIN witness and its preferred
2x acceptance target are external to this checkout, so #372 remains open and
no downstream throughput claim is made.

Evidence: `cargo run --release -p chematic-smiles --features
canonical-search-instrumentation --example canonical_orbit_perf`.

## P3 Explorer browser adversarial smoke — completed local slice

Added `scripts/explorer_browser_smoke.mjs` and attached it to the existing
browser compatibility workflow. The smoke loads the Explorer over HTTP in a
real headless browser, submits valid and malformed pasted records, checks the
stable `2 loaded, 1 failed` status and result count, then starts a bounded
  2,001-record parse and verifies that the 2,000-record display cap is reported,
  then the Cancel action stops before the capped workload completes
without page or console errors.

Evidence:

- `node scripts/explorer_browser_smoke.mjs chromium` — passed.
- `.github/workflows/browser-compat.yml` now runs the same smoke for the
  existing Chromium, Firefox, and WebKit matrix entries.
- A local HTTP in-app browser replay also verified the exact empty-paste error
  (`No SMILES found in the pasted text.`) and recovery to `1 molecule loaded.`;
  the direct Playwright Chromium smoke also passes after waiting for the
  asynchronous Cancel button transition.

This closes only the local Explorer cancellation/malformed-record/display-limit
and stable empty-paste recovery slice. It does not claim local Firefox/WebKit
execution, full error-envelope coverage, or completion of the broader agent
adversarial matrix.
