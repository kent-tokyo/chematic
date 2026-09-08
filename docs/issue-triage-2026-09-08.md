# Open issue triage — 2026-09-08

This snapshot records the bounded completion of issue #507. It does not turn
format-specific safety cases into a throughput, compatibility, or universal
parser claim.

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

Evidence:

- `cargo test -p chematic-perception --lib --offline` — 204 passed, 1 ignored.
- `cargo test -p chematic-ff --lib --offline` — 202 passed.

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

Evidence: `validation/results/mmff94_coverage_227_term_audit_summary.json`,
`validation/results/mmff94_coverage_227_root_cause_classification.json`, and
the provenance decision in `scripts/mmff94_provenance/PROVENANCE.md`.

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
