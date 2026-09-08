# Open issue triage — 2026-09-08

This snapshot records the bounded completion of issue #507. It does not turn
format-specific safety cases into a throughput, compatibility, or universal
parser claim.

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

This is validation hardening only. Full upstream-backed RXN dialect support,
including any format-specific metadata beyond the current loss-aware V2000
adapter, remains open.

Evidence:

- `cargo test -p chematic-rxn document --offline` — 4 passed.
- `cargo test -p chematic-mol rxn --offline` — 9 passed.
- `git diff --check` — passed.

## Issue #461 — CDXML document attribute safety boundary

The loss-preserving CDXML document adapter now applies its configured
`max_attribute_bytes` limit consistently to document, page, and object
attributes. Oversized presentation metadata is rejected with the existing
typed resource-limit error instead of being parsed without the advertised
bound.
Document edits also reject empty or non-element replacement/insertion
payloads before mutating the source representation.
Attribute-edit commands now validate attribute names before serialization, so
quotes, whitespace, and markup characters cannot create malformed XML.
Edits preserve the source's CRLF/LF convention and whether the source ended
with a newline, avoiding unrelated representation churn.
Documents parsed with explicit resource limits retain those limits for later
edits and revalidation, so an edit cannot silently fall back to defaults.
If a document contains duplicate page IDs, edits now reject that ambiguous
target instead of silently modifying the first matching page.

This is a resource-safety slice only. Full ChemDraw presentation semantics,
including every style, geometry, grouping, and annotation dialect, remain
outside the current adapter contract.

Evidence:

- `cargo test -p chematic-mol cdxml_document --offline` — 11 passed.
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
XYZ, V3000, MOL2, CML, CDXML, mmCIF, and PDB. The dependency-free gate checks
all 90 cases and requires exactly one failed record for each case.

The same gate retains nine oversized-input rejections and 18 gzip controls:
one valid decompressed control and one post-decompression input-limit rejection
for every format. CML, CDXML, and PDB remain explicitly line-limit safety
checks because their current readers are intentionally lenient about unknown,
empty, or non-record input.

## Evidence

- `python3 scripts/check_streaming_format_limits.py` — 90 negative cases,
  9 oversized cases, and 18 gzip cases passed.
- The corpus schema and exact ten-case-per-format count are validated before
  any runner invocation.

This closes the bounded malformed-corpus expansion in #507. Cross-language
streaming parity, equivalent cross-engine throughput, and broader parser
semantics remain separate open roadmap gates.

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
correctness mismatches, zero search-budget exhaustions, and an 8.48x Tier A
geometric-mean speedup (Tier B negative control: 2.33x). Across Tier A the
exhaustive engine visited 6,186 leaves while the orbit-pruned engine wrote 13
leaves, visited 62 nodes, and performed 46 orbit tests; the repeated multi-Boc
and multi-pivaloyl fixtures each collapsed to one leaf from 432 exhaustive
leaves.

This is local proxy evidence only. The exact RENKIN witness and its preferred
2x acceptance target are external to this checkout, so #372 remains open and
no downstream throughput claim is made.

Evidence: `cargo run --release -p chematic-smiles --features
canonical-search-instrumentation --example canonical_orbit_perf`.
