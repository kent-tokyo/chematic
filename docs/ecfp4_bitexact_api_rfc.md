# RFC: ECFP4 bit-exactness parameter matrix + API rollout design

**Status:** diagnosis/design (below) merged as-is; the rollout it designed has since
shipped on `feat/ecfp4-bitexact-stable-api` (see `CHANGELOG.md`'s `[Unreleased]`
entry): Python (`Mol.rdkit_ecfp4`/`rdkit_ecfp4_detail`/`rdkit_ecfp_config`/
`rdkit_ecfp_config_detail`) and WASM
(`rdkit_ecfp4_bitvec`/`rdkit_ecfp4_detail_json`/`rdkit_ecfp_config_bitvec`/
`rdkit_ecfp_config_detail_json`) bindings for the exact config diagnosed here (§2.3/
§2.4), plus a generalized, independently-oracle-reverified `(radius, fpSize)` matrix
(`chematic_fp::rdkit_morgan_fingerprint`/`RdkitMorganConfig`,
`crates/chematic-fp/src/rdkit_morgan_config.rs`) and a shared cross-language fixture
corpus (`validation/ecfp4_rdkit_stable_api_fixtures.json`, §2.6). `useChirality`/
`useBondTypes=false`/alternative atom invariants remain `architecturally_unimplemented`
(§1.4-1.6), unchanged and out of scope for that PR too. The measurement/design content
below is otherwise historical and unedited.

**Branch:** `diag/ecfp4-bitexact-api`, forked from `main` at
`659baca221f71f135ce0e1780e71245d8770f132`.

**Files touched (all new, none under `crates/*/src/**`):**
- `crates/chematic-fp/examples/rdkit_ecfp4_bitexact_matrix_dump.rs` — dumps the
  production API's output plus a `diagnostics`-feature radius sweep for the frozen
  fixture corpus.
- `crates/chematic-fp/Cargo.toml` — one added `[[example]]` entry for the dump above
  (`required-features = ["diagnostics"]`, matching the existing 4 entries).
- `scripts/ecfp4_bitexact_matrix_fixtures.csv` (new) — 33-row fixture corpus, single
  source of truth for both the Rust dump and the Python oracle script.
- `scripts/ecfp4_bitexact_matrix_diagnosis.py` (new) — RDKit-oracle cross-check,
  fail-closed (non-empty `unclassified` exits 1).
- `validation/results/ecfp4_bitexact_matrix_dump.jsonl` (new) — frozen chematic-side
  dump, reproducible via the example above.
- `validation/results/ecfp4_bitexact_matrix_summary.json` (new) — the diagnosis
  script's machine-readable output.
- `docs/ecfp4_bitexact_api_rfc.md` (this file, new).

**Explicitly out of scope / not touched:** anything under `crates/*/src/**` (no
production code modified, no re-implementation of `rdkit_morgan_ecfp4_experimental`
or its supporting hash port); `feat/io-mrv`, `feat/io-tdt`,
`feat/io-smiles-supplier-writer`, `fix/smiles-bracket-implicit-h`,
`diag/stereo-reader-integration-boundary`, `feat/stereo2d-local-parity`,
`diag/aromaticity-rdkit-parity`, `diag/canonical-smiles-residual`,
`diag/etkdg-3d-gap`, `diag/accurate-cip-audit` (other agents' branches — not read
except the aromaticity RFC via `git show`, not rebased onto, not merged from).

**Deliverables:** this RFC (measurement matrix + design sketch); the fixture corpus;
a Rust dump example; a fail-closed Python diagnosis script; a frozen JSON summary.

**Done condition:** every one of the ~180 measured matrix cells (33 fixtures ×
molecule-shape/radius/nBits/count axes, plus 3 oracle-only capability checks)
classifies into a named bucket with zero `mismatch_unclassified` rows (verified: the
diagnosis script's own `--summary-out` reports `"unclassified_count": 0` and exits 0),
and the design sketch in Part 2 covers all six items the task specified.

## 0. Headline, up front

The production API, `chematic_fp::rdkit_morgan_ecfp4_experimental`
(`crates/chematic-fp/src/rdkit_morgan_ecfp4.rs`), exposes exactly **one point** in the
full RDKit Morgan/ECFP parameter space: radius = 2 (ECFP4), fpSize = 2048,
`includeRedundantEnvironments = false`, `useChirality = false`, `useBondTypes = true`,
RDKit's default `GetConnectivityInvariants` atom invariant. There is no `EcfpConfig`-
style struct, no radius/nBits/chirality argument — the function signature is
`fn(&Molecule) -> Result<RdkitMorganEcfp4, RdkitMorganError>`, full stop.

Given that, **most of the requested parameter matrix is not a pass/fail measurement —
it's an "exposed vs. not exposed" classification**, and reporting it as a blended
percentage would misrepresent what was actually tested. This RFC uses four states
instead, per the advisor's framing, confirmed against the module source before use
(not assumed):

| State | Meaning |
|---|---|
| `verified_bit_exact` | Driven end-to-end against a live RDKit oracle at the production API's one real config; matches. |
| `verified_reachable_via_internal_math` | The port's underlying hash machinery (`expand_one_pass`/`rdkit_morgan_raw_trace`, `diagnostics` feature) generalizes to this cell (e.g. radius ≠ 2) even though no production entry point exposes it; matches the oracle when driven directly. |
| `verified_via_postfold_of_public_data` | Derived by folding the production API's already-public `sparse_counts` field in Python (`raw_id % N`) — no source touched — and matching the oracle's real fold. |
| `architecturally_unimplemented` | No code path in `rdkit_morgan_hash.rs` can express this option at all (confirmed by reading `connectivity_invariant`/`checked_bond_invariant`), independent of whether it would match RDKit if it existed. The oracle check for these cells only demonstrates that real RDKit *does* distinguish the option — i.e. the gap is consequential, not theoretical. |
| `mismatch_aromaticity_kekulize_hardfail` | Traced to the `diag/aromaticity-rdkit-parity` PR's already-identified `kekulize()` hard-fail classes (tropylium cation, imidazolium, pyridinium, pyrylium, tellurophene, phosphole). |
| `mismatch_unclassified` | Anything else. **Zero rows landed here** (see §1.6) — the diagnosis script exits non-zero if this bucket is ever non-empty. |

Methodology note: M4-A0's own 5,046/5,046-molecule corpus validation (cited in the
module docs) already established bit-exactness at the one supported config at scale.
This RFC does **not** re-run that corpus — it builds a small (33-fixture), deliberately
constructed corpus targeting the axes M4-A0 never swept (radius, nBits, count,
chirality, bond-types, alternative invariants, disconnected/isotope/charged/stereo
molecule shapes) and cross-checks each against a live oracle (`rdkit==2026.03.3`,
matching the module's own pinned commit `8afba32ec539dcb2369bc84549d802aca3f7eb39`).

## 1. Measurement matrix

### 1.1 Molecule-shape axis (production's one real config: r=2, 2048 bits, binary)

Every fixture run through `rdkit_morgan_ecfp4_experimental` directly, compared against
an RDKit oracle built with `rdFingerprintGenerator.GetMorganGenerator(radius=2,
fpSize=2048, includeChirality=False, useBondTypes=True,
includeRedundantEnvironments=False)` on three independent signals at once: the
unfolded `(atom, radius) → raw_id` pairs, the unfolded sparse counts, and the folded
2048-bit on-bit set.

| Shape bucket | Fixtures | `verified_bit_exact` | `mismatch_aromaticity_kekulize_hardfail` | `mismatch_unclassified` |
|---|---|---|---|---|
| baseline (benzene, ethane, propane) | 3 | 3 | 0 | 0 |
| disconnected (multi-fragment SMILES) | 3 | 3 | 0 | 0 |
| isotope-labeled (¹³C, ²H) | 4 | 4 | 0 | 0 |
| charged, kekulizes OK (acetate, ammonium, nitro, sulfate…) | 6 | 6 | 0 | 0 |
| charged, kekulize hard-fails (tropylium/imidazolium/pyridinium/pyrylium/Te-phene/phosphole) | 6 | 0 | 6 | 0 |
| aromatic vs. Kekulé input, same molecule (pyridine/naphthalene/furan/thiophene pairs) | 8 | 8 | 0 | 0 |
| stereo, tetrahedral (L-/D-alanine) | 2 | 2 | 0 | 0 |
| stereo, E/Z (2-butene) | 2 | 2 | 0 | 0 |
| **Total** | **33 (27 non-kekulize-fail + 6 kekulize-fail)** | **27/27** | **6/6** | **0** |

**Disconnected molecules, isotopes, charged atoms (the kekulizable ones), and stereo-
marked atoms (both tetrahedral and E/Z) are all bit-exact at the one supported
config** — 27/27, no exceptions, no silent drops. **Aromatic vs. Kekulé input for the
same molecule** also matches 8/8: `rdkit_morgan_ecfp4_experimental` internally calls
`apply_aromaticity_rdkit_parity_experimental` regardless of the input molecule's
existing flags, so both spellings normalize to the same result and both match RDKit
(RDKit's own `Chem.MolFromSmiles` always sanitizes before hashing too, so this is
matching RDKit's real precondition, not a chematic-specific accommodation).

The 6 charged-aromatic fixtures that fail are **exactly** the six molecule classes the
`diag/aromaticity-rdkit-parity` PR already identified as `chematic_core::kekulize()`
hard-fail cases (verified independently here, not merely cited: this diagnosis's own
Rust dump reproduces the identical error, e.g. fixture 24 (tropylium,
`c1ccc[cH+]cc1`) reports `"rdkit-exact ecfp4: aromaticity: rdkit-parity aromaticity:
kekulization failed: atom 6 (C) cannot be assigned a double bond"` — the same atom
index and reason the other RFC's §1a table gives). This is the production API's
*intended* behavior on this input class, not a new bug: the module's own doc comment
states there is "no fallback to another aromaticity engine anywhere in this module's
public path," specifically because a silent Hückel fallback would invalidate the
bit-exactness claim exactly where a caller most needs to know it doesn't hold (see the
module's `hueckel_fallback_would_be_detectable_if_silently_reintroduced` test). So this
bucket is correctly attributed to the other RFC's root cause, not double-counted as a
new ECFP4-specific defect — but it is a real, currently-unfixed **coverage gap**: any
caller whose corpus contains tropylium/imidazolium/pyridinium/pyrylium/tellurophene/
phosphole-type rings gets `Err` from the bit-exact path today, with no bit-exact
fallback available (see Part 2 §2 for what "explicit non-silent fallback" should mean
here).

### 1.2 Radius axis (0, 1, 2, 3)

Production only computes radius 2. To test whether the *underlying hash machinery*
(not any production entry point) generalizes, this diagnosis drives
`chematic_fp::diagnostics::rdkit_morgan_raw_trace` (the `diagnostics`-feature function
`rdkit_morgan_ecfp4_experimental` itself reuses internally for its one fixed radius)
directly at `max_radius = 3`, on the same aromaticity engine the production path uses
(`apply_aromaticity_rdkit_parity_experimental`, not plain Hückel — chosen specifically
so a mismatch here can't be misattributed to a different aromaticity engine), and
compares the `raw_identifier_default` lifecycle only (RDKit's real
`includeRedundantEnvironments=false` mode — the `full` lifecycle is a distinct
diagnostic surface, not compared).

| Radius | Fixtures tested (aromaticity succeeded) | Match | N/A (aromaticity failed upstream) | Fixtures with a *non-empty* radius environment |
|---|---|---|---|---|
| 0 | 27 | 27/27 | 6 | 27 (unconditional, every atom) |
| 1 | 27 | 27/27 | 6 | 27 (every non-degree-0 atom) |
| 2 | 27 | 27/27 | 6 | 17 |
| 3 | 27 | 27/27 | 6 | 13 |

**Honesty check on the "27/27" figures at radius ≥ 2**: an empty-vs-empty comparison
(an atom fully degree-0-dead or fully suppressed by that radius on both sides) counts
as a trivial match, not evidence of generalization — most of this corpus's small
molecules (ethane, propane, the charged-ok ions) legitimately have nothing left to
compute by radius 3. Counted directly from the frozen dump: **13 of 27** fixtures have
at least one genuinely non-empty radius-3 `(atom, radius) → raw_id` pair, spanning
benzene, both spellings each of pyridine/naphthalene/furan/thiophene, both disconnected
fixtures, deuterated benzene, and nitrobenzene — not just one large ring system
(naphthalene alone contributes 9 of the 33 total radius-3 pairs, but the other 12
fixtures contribute real, independently-checked non-empty pairs too). So the radius-3
match is real, non-vacuous evidence across a real spread of ring sizes and molecule
shapes, not a single lucky case — but it is thinner evidence than radius 0/1 (which are
richly non-empty for literally every fixture) or the production-locked radius 2 (17/27
non-empty). The architectural argument (`expand_one_pass` has no radius-specific logic,
only a loop bound) is what actually carries the "radius as a parameter" conclusion in
Part 2 §1 — this data corroborates it without contradiction, but a future radius-4+
rollout should widen the corpus with more polycyclic fixtures before treating radius
generalization as fully proven past radius 3.

**Classification: `verified_reachable_via_internal_math`.** Radius 0, 1, and 3 all
match the oracle bit-for-bit on every fixture where aromaticity perception succeeded
(the same 27/33, including disconnected and isotope-labeled fixtures — the radius
sweep is not restricted to simple rings). This directly contradicts the module's own
"never compared radius 3 against the oracle" disclaimer (`rdkit_morgan_ecfp4.rs`'s
doc comment, line 7-9) as a statement about the *underlying port* — it was true about
the *production wrapper* (which never called the hash port at radius ≠ 2) and remains
true that no production API can be asked for radius 3 today, but the port itself
(`expand_one_pass`, `connectivity_invariant`, `checked_bond_invariant`) has no radius
dependence beyond the loop bound — nothing in the hash formulas is ECFP4-specific.
This is a real, evidence-backed API-design finding: exposing radius as a parameter is
very unlikely to require new hash logic, only a new thin wrapper (see Part 2 §1).

### 1.3 nBits axis (128, 256, 512, 1024, 2048) and count vs. binary

Production hardcodes `fpSize = 2048` and only ever returns a binary `BitVec2048` —
there is no folded *count* fingerprint in the public API at all (only the unfolded
`sparse_counts` field, radius ≤ 2, one raw-id → count map). Before testing anything,
the fold convention itself was verified empirically (not assumed) on naphthalene: for
every `fpSize` in {128, 256, 512, 1024, 2048}, `{raw_id % fpSize for raw_id in
sparse_ids} == set(GetFingerprint(mol).GetOnBits())`, and the analogous count-fold
equality holds for `GetCountFingerprint`. Both held exactly, confirming RDKit's real
Morgan-generator fold is plain modulo, not some other bucketing scheme.

Given that, both axes are testable **without touching any production code**, by
folding the already-public `sparse_counts` field in the diagnosis script itself:

| nBits | Binary fold match | Count fold match |
|---|---|---|
| 128 | 27/27 | 27/27 |
| 256 | 27/27 | 27/27 |
| 512 | 27/27 | 27/27 |
| 1024 | 27/27 | 27/27 |
| 2048 | 27/27 | 27/27 |

**Classification: `verified_via_postfold_of_public_data`.** Every fixture that
succeeds at the one supported config would also be bit-exact at any of these 5 sizes,
in both binary and count representations, if a production entry point folded
`sparse_counts` this way. No new hashing is needed for this axis either — folding is
pure post-processing of data the API already exposes.

### 1.4 useChirality on

**Classification: `architecturally_unimplemented`.** `connectivity_invariant`
(radius-0 atom invariant) and `checked_bond_invariant` (bond invariant) — the two
functions `rdkit_morgan_ecfp4_experimental` calls for every atom/bond — have no
chirality byte, no chirality branch, and no way to receive a chirality flag at all.
This is not an oversight to be diagnosed as a bug; the module's own doc comment states
it directly ("the atom-CIP chirality re-fold is likewise out of scope,
`includeChirality=false` pinned throughout this workstream"). Confirmed the option is
real and consequential (not a theoretical one) via the oracle directly on the L-/D-
alanine fixture pair: RDKit's own `includeChirality=False` (matching chematic's only
mode) gives identical fingerprints for the enantiomer pair (`True`), but
`includeChirality=True` gives *different* fingerprints for the same pair (`False`) —
i.e., a real caller who needs stereo-sensitive ECFP4 (a common cheminformatics
requirement, e.g. distinguishing enantiomers in a similarity search) cannot get it from
this API today, silently or otherwise — it's simply not there to ask for.

### 1.5 useBondTypes off

**Classification: `architecturally_unimplemented`.** `checked_bond_invariant`
hard-codes exactly one branch (`useBondTypes=true`); there is no `false` branch, no
parameter to select one. Confirmed consequential the same way: RDKit's own
`useBondTypes=True` vs `False` generators produce different fingerprints for the same
molecule (pyridine, `True` differs from `False`), so this is a real, requestable RDKit
option chematic's bit-exact path cannot express.

### 1.6 Alternative atom invariants (e.g. FCFP-style feature invariants)

**Classification: `architecturally_unimplemented`.** `connectivity_invariant` hard-
codes exactly RDKit's default `GetConnectivityInvariants` component set
(`[atomicNum, totalDegree, totalNumHs, formalCharge, deltaMass, ring?]`); there is no
alternative invariant function, and no parameter to select one, anywhere in
`rdkit_morgan_hash.rs`. (Chematic's separate, non-bit-exact legacy path,
`chematic_fp::ecfp::EcfpInvariantMode`, does offer a second mode, `RdkitMorgan` — but
that mode is explicitly documented as *partition*-parity only, not raw-hash parity;
its own doc comment: "not a claim of RDKit fingerprint bit-compatibility." It's a
different code path with a different, FNV-1a-based hash function, so it cannot serve
as a bit-exact FCFP-style option either.) Confirmed consequential: RDKit's default
invariant generator vs. `GetMorganFeatureAtomInvGen()` (its FCFP-style feature
invariant) give different fingerprints for phenol.

### 1.7 Denominator discipline and the fail-closed check

Per this task's requirement (no silent drops, no single blended percentage): the
diagnosis script's summary JSON reports `unclassified_count` explicitly, and the
script's own exit code is non-zero if that count is ever nonzero. On the frozen run
committed here, **`unclassified_count == 0`** — every one of the 33 fixtures × every
axis it participates in landed in a named bucket, with the fixed six kekulize-hardfail
molecules correctly and exclusively populating the `mismatch_aromaticity_*` bucket (not
pooled into `verified_bit_exact`, not silently dropped from the denominator).

## 2. API design sketch (planning only — nothing below is implemented)

### 2.1 Rust default vs. opt-in

**Recommendation: keep it opt-in, do not make it the Rust default**, for two reasons
independent of each other:

1. **It's fallible.** `rdkit_morgan_ecfp4_experimental(&Molecule) -> Result<_, _>`
   already has a different signature shape than `ecfp4(&Molecule) -> BitVec2048`
   (infallible). Making the fallible path "the default" would either force every
   existing caller of `ecfp4()`/`mol.ecfp4()` to start handling `Result`, or require a
   silent `.unwrap_or_else(|| legacy_ecfp4())` fallback — which is precisely the
   fallback-pooling shape the module's own tests exist to catch (§0's
   `hueckel_fallback_would_be_detectable_if_silently_reintroduced`). §1.1 shows this
   isn't a corner case either — 6/33 hand-built fixtures (charged aromatic
   heterocycles, a real and not-uncommon structural class) fail today.
2. **`_experimental` is in the name for a reason.** Promoting an `_experimental`-named
   function to be *the* default silently changes what every existing caller's `ecfp4()`
   call means without them opting in, which is exactly the kind of silent-replacement
   the task explicitly rules out in item (2) below.

Keep `rdkit_morgan_ecfp4_experimental` as a distinctly-named, explicitly-opted-into
function (rename off `_experimental` only after a version or two of real-world use, per
normal chematic convention — not part of this RFC's scope to decide the exact name/
timing).

### 2.2 Legacy path stays available, explicitly, never silently replaced

This is already true today and should stay true: `chematic_fp::ecfp`/`ecfp4`/`ecfp6`/
`ecfp_with_invariant_mode` (the FNV-1a-based, non-RDKit-hash-exact path,
`crates/chematic-fp/src/ecfp.rs`) are a structurally separate module from
`rdkit_morgan_ecfp4.rs`/`rdkit_morgan_hash.rs`, with no shared hashing code (confirmed
by reading both — the only genuinely hash-independent code either module reuses from
the other is `BondSet`/degree/H-count/isotope-delta *value* computations, not any
hashing or invariant-assembly logic). Concretely, going forward:

- Never repoint `ecfp4()`/`mol.ecfp4()` at the RDKit-exact path internally — the byte
  layout, hash function (FNV-1a vs. RDKit's `hash_combine`), and bond-invariant
  encoding (`aromatic=4` vs. RDKit's `12`) are deliberately different by design (see
  `rdkit_morgan_hash.rs`'s own doc comment on this exact collision risk), so any stored
  fingerprint computed with today's `ecfp4()` must keep meaning the same bits after any
  future change.
- Any new bit-exact entry point gets its own name (`rdkit_ecfp4`,
  `morgan_fingerprint_rdkit_compatible`, or similar — exact naming out of scope here),
  never an overload or a config flag silently bolted onto the existing `ecfp()`
  function signature.

### 2.3 Python exposure

**Currently not exposed at all.** Verified directly: no `rdkit_morgan_ecfp4` /
`RdkitMorganEcfp4` / `rdkit_morgan` string appears anywhere under
`crates/chematic-py/src/`. Only the legacy path is bound (`mol.ecfp4()`,
`mol.ecfp4_chiral()`, `mol.ecfp4_numpy()`, `bulk.ecfp4()` — all `chematic_fp::ecfp4`/
`ecfp`, `crates/chematic-py/src/mol_methods.rs` and `bulk.rs`).

Rollout sketch, given a `Result`-returning Rust function needs a Python-idiomatic
fallible shape:
- New method, distinctly named from `ecfp4()` (e.g. `mol.rdkit_ecfp4()` or
  `mol.morgan_fingerprint_rdkit_exact()`), returning either the folded bit vector
  (`bytes`/`numpy` array, matching `ecfp4_numpy()`'s existing convention) on success or
  raising a distinct Python exception class (mapping `RdkitMorganError`'s three variants
  — `Aromaticity`, `UnsupportedBondOrder`, `InternalInvariantViolation` — to distinct,
  catchable exception types or at minimum a structured `.args` payload, not a bare
  string) rather than returning `None`/a sentinel on failure — silent `None` on failure
  reintroduces exactly the "caller can't tell bit-exactness didn't hold" problem the
  Rust module works to avoid.
- Expose `sparse_counts`/`raw_bit_info`/`folded_bit_info` as a second method or an
  optional-detail flag (mirroring `ecfp_with_bitinfo`'s existing split from `ecfp4()`),
  not bolted onto the main return type — most callers want just the fingerprint.
- `bulk.rdkit_ecfp4()` (batch numpy array), matching `bulk.ecfp4()`'s existing shape,
  needs an explicit per-row failure story: either a boolean success mask alongside the
  array, or skip-with-warning is explicitly rejected (silent drop) — this needs its own
  design pass, flagged here rather than resolved.

### 2.4 WASM exposure

**Currently not exposed at all**, same as Python — verified directly: no
`rdkit_morgan`/`RdkitMorganEcfp4` reference anywhere under `crates/chematic-wasm/src/`.
`crates/chematic-wasm/src/mol_fingerprints.rs` only binds the legacy path
(`ecfp4_bitvec`, `ecfp4_bitvec_with_chirality`, `ecfp6_bitvec`, `ecfp_bitvec_custom`,
etc. — all thin wrappers over `chematic_fp::ecfp`/`ecfp4`/`ecfp6`).

Rollout sketch: `wasm-bindgen` has no native `Result` ergonomics as clean as PyO3's
exception mapping — the existing WASM bindings that can fail already use
`Result<T, JsValue>` (see `ecfp_bitvec_custom`, `virtual_screen_ecfp4_json`-adjacent
functions elsewhere in the crate) and construct a `JsValue` error via
`JsValue::from_str` or similar. Follow that existing convention: new function returns
`Result<Vec<u8>, JsValue>`, with the `JsValue`'s string payload carrying
`RdkitMorganError`'s `Display` text (already implemented, see §0) rather than a generic
"failed" message, so a JS caller can at least log/branch on *why*.

### 2.5 Documentation / versioning / changelog needs

- **CLAUDE.md's `chematic-fp` row** (`ECFP/FCFP, MACCS, MAP4, Tanimoto`) doesn't
  currently distinguish the two ECFP4 implementations at all; once a Python/WASM entry
  point exists, it needs a one-line pointer to both (bit-exact vs. legacy) and *why*
  two exist, mirroring this RFC's §2.2.
- **`docs/validation.md`** (the 19-descriptor RDKit-parity breakdown CLAUDE.md already
  references) is the natural home for this RFC's §1 matrix once any of it becomes a
  real, shipped API surface — not before, to avoid documenting a capability no caller
  can reach yet.
- **Changelog entry, whenever any of §2.3/§2.4 ships**, must state the exact supported
  config explicitly (radius=2 only, 2048 bits only, no chirality, no alternative
  invariants) — not just "RDKit-bit-exact ECFP4," which this diagnosis's own §1.4-1.6
  findings show would overclaim.
- **This RFC's own six-fixture kekulize-hardfail gap (§1.1)** should be cross-linked
  from whatever tracking issue exists for the `diag/aromaticity-rdkit-parity` PR's
  `kekulize()` fix — fixing that root cause automatically closes this gap too; it
  should not be independently re-implemented here.

### 2.6 Cross-language identical-bit test suite

What it needs to check, concretely, given everything measured in §1:

1. **One shared, versioned fixture corpus** (SMILES + expected raw ids/folded bits),
   generated once against the pinned RDKit oracle and checked into the repo (not
   regenerated ad hoc per language) — this RFC's `scripts/ecfp4_bitexact_matrix_fixtures.csv`
   plus its dump/oracle pair is a template for the shape such a suite would need, not a
   claim that today's 33-fixture corpus is sufficient for a permanent regression gate
   (it was built to span *diagnosis* axes, not to be an exhaustive corpus).
2. **Per-language driver that calls the real public entry point**, not an internal
   function — Rust via the crate directly, Python via the PyO3-bound method (§2.3),
   WASM via a headless JS runtime (e.g. `wasm-bindgen-test` or Node) calling the bound
   function (§2.4) — so the test catches a binding-layer bug (wrong byte order, off-by-
   one in the exception mapping, wrong numpy dtype) that a Rust-only test cannot see.
3. **Identical fixture ids and identical expected values across all three drivers** —
   one canonical JSON/CSV expectations file all three read, not three independently
   maintained copies that can silently drift (the same "single source of truth" this
   RFC's own fixtures CSV already does for its two consumers).
4. **Explicit failure-config fixtures included on purpose** (the six kekulize-hardfail
   molecules from §1.1), asserting all three bindings raise/return the *same* error
   category — not just that the success cases match. A cross-language suite that only
   tests success paths would miss a binding silently swallowing the `Result`'s `Err`
   arm into a wrong/empty success value in exactly one language.
5. **A version-pin check**: since the whole guarantee is pinned to a specific RDKit
   release (`Release_2026_03_4`), the suite should fail loudly (not silently drift) if
   that pin is ever bumped without re-validating — this project already has the
   pattern (`gen_ecfp_rdkit_environment_oracle.py`'s SHA-pinning discipline) to extend.
