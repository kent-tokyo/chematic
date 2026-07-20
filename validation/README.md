# chematic Validation

Documented evidence that chematic's descriptors agree with industry-standard tools.

## Corpora

### 175-mol drug-like corpus

A curated set of 175 drug-like molecules covering common scaffolds (benzoic acid derivatives,
heterocycles, amino acids, steroids, macrolides). Used for per-descriptor regression testing.

- **File:** `scripts/rdkit_ref_properties.tsv` (175 rows)
- **Columns:** name, smiles, mw, logp, tpsa, hac, hbd, hba
- **Reference tool:** RDKit 2026.03.3
- **How to regenerate:** `python scripts/gen_rdkit_reference.py`

### 4,999-mol ChEMBL subset

A random sample from ChEMBL used for large-scale agreement testing on HBA, HBD, and aromatic ring count.

- **File:** external (not committed; requires download)
- **Reproduce:** `python scripts/bench5k.py ~/Downloads/SMILES.csv`

### Morgan/ECFP RDKit environment-parity diagnostic (5,000-mol corpus + 41 fixtures)

Locates, per molecule, the first stage at which chematic's Morgan/ECFP expansion diverges
from RDKit's (radius-0/1/2 invariants, redundant-environment suppression, sparse counts,
2048-bit folding, bitInfo). Diagnostic only -- production `ecfp4()`/`ecfp6()`/
`morgan_fp_counts()` are unchanged.

- **Files:** `ecfp_rdkit_environment_parity_manifest.json`, `_summary.json`,
  `_rows.jsonl` (41 edge-fixture molecules), `_first_divergence.tsv` (full 5,041-input run)
- **Reference tool:** RDKit 2026.03.3 (`rdFingerprintGenerator.GetMorganGenerator`)
- **How to regenerate:** `python scripts/gen_ecfp_rdkit_environment_oracle.py` +
  `cargo run -p chematic-fp --release --features diagnostics --example
  morgan_rdkit_environment_trace` + `python scripts/ecfp_rdkit_environment_parity.py`
  (see each script's docstring for exact invocation)

### Morgan/ECFP RDKit environment-suppression parity (Phase B, 5,041+4-input set)

Whether chematic's RDKit-equivalent redundant-environment suppression
(`crates/chematic-fp/src/morgan_environment.rs`, `SuppressRdkitRedundant` mode
-- additive/experimental; production `ecfp4()`/`ecfp6()`/`morgan_fp_counts()`
unchanged) emits the same set of `(atom_idx, radius)` environments, and the
same raw-identifier sparse-count *shape*, as RDKit's own default
(`includeRedundantEnvironments=False`) generator, on PR #120's original
5,041-input set plus 4 pinned representative-swap fixtures (5,045 total; see
`scripts/ecfp_rdkit_suppression_representative_swap_fixtures.csv`).

Implementation verified directly against RDKit's real C++ source: commit
[`0062b670640352ab63d6256be608615e87e1af53`](https://github.com/rdkit/rdkit/blob/0062b670640352ab63d6256be608615e87e1af53/Code/GraphMol/Fingerprints/MorganGenerator.cpp),
`MorganEnvGenerator<OutputType>::getEnvironments` -- a specific commit SHA,
not a mutable `master` reference.

**Results:**

| Metric | Result |
|---|---|
| Emitted `(atom_idx, radius)`-pair-set exact match | 5,032/5,045 (99.74%) |
| Raw-identifier sparse-count *shape* exact match (multiset of per-id emission counts) | 5,044/5,045 (99.98%) |
| `sparse_count_mismatch` fixtures (8, from the Phase A diagnostic) now shape-resolved | 8/8 |
| Tanimoto-vs-RDKit Pearson r, before (`ecfp4_rdkit_invariants`) → after (`ecfp4_rdkit_environment_experimental`) | 0.9479 → 0.9547 (Δ+0.0068, improved; n=300 sample, seed=42, 44,850 pairs -- non-gating reference) |
| Full-corpus (5,000 mol) wall time, baseline → suppression (median of 5 independent process runs) | 1.315s → 1.508s (1.146x) |
| Full-corpus peak RSS, baseline → suppression (median of 5 runs, `/usr/bin/time -l`) | 18.7 MiB → 19.3 MiB (1.030x) |

Pair-set mismatches: 13 of 5,045 (the same 9 residuals from the original
5,041-input run, plus their 4 pinned duplicates), all single-pair swaps at
the same radius -- **not** a claim that the swapped atoms are chemically
equivalent or near-equivalent, and **not** a claim that the two candidates
provably compute the identical cumulative bond environment (that would
require diagnosing raw bond-index-sets directly, which this validation
doesn't do). What's actually measured, precisely: two different atoms
produce the same *raw identifier*, and the selected representative differs
because RDKit and chematic currently order those candidates using different
hash values (FNV-1a vs RDKit's own hash never match by construction -- same
"not bit-compatible, partition/set-only" scope as every other RDKit-parity
mode in this crate). See the pinned fixtures for concrete cases: `CC(=O)NO`
(atoms 1 vs 3, not a symmetric pair), an isotope-labeled methyl pair, a
steroid-like fused-ring epoxide, and a large polycyclic aromatic -- each
verified to be *exactly* a 1-pair swap with total-emitted-count,
sparse-count shape, and unique-*raw-identifier*-count (deliberately not
called "unique bond-environment count" -- a raw identifier can in principle
be shared by two structurally different environments via hash collision, as
the pyridine case below demonstrates) all preserved; only which atom
represents one shared identifier differs. **These 4 fixtures, plus the 8
`sparse_count_mismatch` fixtures, plus every "both"-bucket mismatch anywhere
in the input, are hard GATES in `scripts/ecfp_rdkit_suppression_parity.py`
(nonzero exit on any regression) -- not just reported numbers.**

Sparse-count-shape mismatch: 1 of 5,045, a pair-set *exact match*
(`C1=CC=NC=C1`, Kekulé pyridine) whose count multiplicities still differ --
traced to accidental cross-radius hash collisions that differ between
FNV-1a and RDKit's hash for this molecule's structurally-symmetric ring
carbons, not a suppression-algorithm defect (the underlying emission
*decision*, i.e. which atoms survive at which radii, is provably identical
between the two implementations for this molecule).

- **File:** `ecfp_rdkit_suppression_parity_summary.json`,
  `ecfp_rdkit_suppression_tanimoto_summary.json`
- **Reference tool:** RDKit 2026.03.3 (same oracle rows as the Phase A
  diagnostic above -- `default` variant's `sparse_bit_info`/`sparse_counts`/`folded_on_bits`)
- **How to regenerate:** `cargo run -p chematic-fp --release --features
  diagnostics --example morgan_suppression_dump` +
  `cargo run -p chematic-fp --release --example morgan_suppression_tanimoto_dump`
  + `python scripts/gen_ecfp_rdkit_environment_oracle.py` +
  `python scripts/ecfp_rdkit_suppression_parity.py` +
  `python scripts/ecfp_rdkit_suppression_tanimoto.py` (see each script's
  docstring for exact invocation). Performance record (not a merge gate):
  `cargo run -p chematic-fp --release --example morgan_suppression_benchmark`.

### Morgan M4-A0: RDKit-exact raw-identifier hash port (diagnostic, 5,048-input set)

Diagnostic-only, source-verified port of RDKit's actual Morgan hash-combine
machinery (`crates/chematic-fp/src/rdkit_morgan_hash.rs`) -- unlike every
prior Morgan-parity mode in this crate, which only claims *partition*
agreement (which atoms are chemically equivalent) and *lifecycle* agreement
(who wins/dies under suppression), this compares actual 32-bit hash VALUES
against real RDKit, atom by atom and radius by radius. Not wired into any
production API at the time of this milestone (later promoted to production
by the "Phase B" section below); `ecfp4()`/`ecfp6()`/`ecfp4_rdkit_invariants()`/
`ecfp4_rdkit_environment_experimental()` are all unchanged (production
snapshot verified byte-identical, see Results).

Ported directly from RDKit's real C++ source, commit
[`8afba32ec539dcb2369bc84549d802aca3f7eb39`](https://github.com/rdkit/rdkit/blob/8afba32ec539dcb2369bc84549d802aca3f7eb39/Code/GraphMol/Fingerprints/MorganGenerator.cpp)
(the true resolution of tag `Release_2026_03_4`, independently verified via
the GitHub tags API this session -- see `THIRD_PARTY_NOTICES.md` for the
attribution and the note on two other, imprecise SHAs already in this
project's history under the same tag label, since fixed in place).

**Two aromaticity-preprocessing paths were compared, and their results are
NOT pooled into one number** -- an earlier draft of this section did pool
them (a Hueckel fallback silently substituted for 2 rows where the
RDKit-parity engine failed, then counted as an RDKit-parity "match"),
which is exactly the measurement accident
[`apply_aromaticity_rdkit_parity_experimental`]'s own `Result` contract
exists to prevent. Corrected: every row is tagged with an explicit
`aromaticity_status`, and the RDKit-parity exact-match rate is computed
ONLY over rows where that engine actually succeeded.

**Results (5,048-input set = 5,000-mol corpus + PR #120's 41 fixtures + PR
#123's 4 fixtures + `ecfp_rdkit_m4a0_hash_fixtures.csv`'s 3):**

| Path | Metric | Result |
|---|---|---|
| Production Hueckel aromaticity | Radius-0 numeric exact match | 5,048/5,048 (100%) |
| Production Hueckel aromaticity | Full numeric exact match (radius 0-2, representative selection, sparse counts, folded bits, bitInfo) | 4,989/5,048 (98.83%) |
| RDKit-parity aromaticity (`apply_aromaticity_rdkit_parity_experimental`, no fallback) | Preprocessing succeeded | 5,046/5,048 |
| RDKit-parity aromaticity | Preprocessing failed (`KekulizationFailed`, both pinned as fixtures -- see below) | 2/5,048 |
| RDKit-parity aromaticity | Full numeric exact match **among the 5,046 successful rows** | **5,046/5,046 (100%)** |
| RDKit-parity aromaticity | Non-exact among successful rows | 0 |
| Hueckel control on JUST the 2 error rows (non-gating -- answers "does the OLD path agree with RDKit here", not "does RDKit-parity work here") | Exact match | 2/2 |
| PR #123's 9 unique representative-selection residuals | Resolve to `exact_match` under the RDKit-exact hash alone (no aromaticity-engine swap needed) | 9/9 |
| PR #123's Kekule-pyridine sparse-count-shape mismatch (`C1=CC=NC=C1`) | Resolves to `exact_match`; confirms the documented root cause (FNV-1a-specific hash collision, not a suppression defect) | resolved |
| Production API byte-identical (`ecfp_regression_snapshot`, before/after SHA-256) | confirmed | 0 change |
| Oracle regeneration determinism (`--verify-determinism`, full 5,048 input) | byte-identical across two runs | confirmed |
| Positive controls (radius-0/1 identifier, bond invariant, 32-bit-wrapping removal, representative swap, folded-bit, dropped row, duplicate row ID) | all correctly cause non-zero exit; reverted, never committed | 8/8 |

**Cross-referencing the 59 Hueckel-path residuals against the RDKit-parity
path, row by row (not just comparing aggregate counts):** all 59 had
RDKit-parity preprocessing succeed, and all 59 became `exact_match` under
it -- `resolved_by_rdkit_parity: 59, not_evaluable_due_to_aromaticity_error:
0, still_mismatching: 0`. Neither of the 2 RDKit-parity error rows
overlaps with the 59 Hueckel residuals (they were already exact matches
under Hueckel).

The 59-row residual under production Hueckel aromaticity traces to ONE
mechanism: chematic's Hueckel-based aromaticity *perception* disagreeing
with RDKit's own aromaticity model on specific fused/macrocyclic ring
systems (e.g. `C[Si](C)(C)c1ccc(C2=Cc3ccccc3C3=NCCCN23)cc1`) -- not a hash
defect. `apply_aromaticity_rdkit_parity_experimental`
(`crates/chematic-perception/src/rdkit_parity.rs`, built for exactly this
kind of disagreement, in an earlier milestone) resolves it.

**The 2 RDKit-parity preprocessing failures are pinned as permanent
fixtures**, not just recorded in a JSON summary --
`scripts/ecfp_rdkit_m4a0_rdkit_parity_kekulization_gap_fixtures.csv`, plus
`chematic-perception::rdkit_parity::tests::known_kekulize_gap_protonated_pyridinium`
(the second is new; `Cc1cn2c(=O)c3ncn(COCCO)c3nc2n1C` was already a pinned
gap case, `production_api_does_not_mutate_input_on_failure`) -- so a future
kekulization fix is verified against the actual engine returning `Ok`, not
just a number changing. Neither fixture is a corpus/fixture duplicate
(each SMILES appears exactly once across the whole 5,048-input set).

A real bug in this diagnostic's own trace logic was found and fixed during
this milestone (not a hash defect either): an early version shared one
`dead`-atom array between RDKit's two `includeRedundantEnvironments`
lifecycles, so an atom suppressed under the `default` lifecycle silently
stopped being *computed* under the `full` lifecycle too in later rounds --
caught by the full-corpus comparator (hand verification on a handful of
fixtures missed it, since it only checked value equality on entries present
on both sides, not entry *count*). Fixed by running two fully independent
passes and merging by `(atom, radius)` key; regression-pinned in
`rdkit_morgan_hash.rs`'s own test suite.

- **Files:** `ecfp_rdkit_raw_identifier_parity_summary.json` (production
  Hueckel aromaticity run), `ecfp_rdkit_raw_identifier_parity_aromaticity_variant_summary.json`
  (RDKit-parity aromaticity engine, full corpus, honest success/error
  denominators), `ecfp_rdkit_raw_identifier_parity_oracle_manifest.json`
- **Reference tool:** RDKit 2026.03.3 (`rdFingerprintGenerator.GetMorganGenerator`,
  `includeRedundantEnvironments` True and False variants; same pinned option
  set as every other Morgan-parity mode in this crate)
- **How to regenerate:** `python scripts/gen_ecfp_rdkit_environment_oracle.py`
  + `cargo run -p chematic-fp --release --features diagnostics --example
  rdkit_morgan_hash_dump` + `python scripts/ecfp_rdkit_raw_identifier_parity.py`
  (see each script's docstring for exact invocation). RDKit-parity-engine
  comparison (no Hueckel fallback): `cargo run -p chematic-fp --release
  --features diagnostics --example rdkit_morgan_hash_dump_aromaticity_variant`
  + `python scripts/ecfp_rdkit_raw_identifier_parity_aromaticity_variant.py`.

**Implemented as Phase B, same day (2026-07-20)** -- see the section
immediately below for the production API and its own full corpus results.

### Phase B: `rdkit_morgan_ecfp4_experimental` -- production, fallible, RDKit-bit-exact ECFP4

Promotes M4-A0's reference engine (`crates/chematic-fp/src/rdkit_morgan_hash.rs`)
to a real public API in `crates/chematic-fp/src/rdkit_morgan_ecfp4.rs`:
`pub fn rdkit_morgan_ecfp4_experimental(mol: &Molecule) -> Result<RdkitMorganEcfp4, RdkitMorganError>`.
Scope is intentionally narrow, matching exactly what M4-A0 verified numerically:

- **radius = 2 (ECFP4) only** -- not ECFP6/radius = 3; M4-A0 never compared radius 3
  against the oracle, so claiming bit-exactness there would be unverified.
- Uses `apply_aromaticity_rdkit_parity_experimental` internally as a fallible `Result`
  step -- **no Hueckel fallback anywhere in the public path.** No entry point accepts a
  pre-aromatized `Molecule` (would let a caller bypass the engine and silently lose the
  bit-exactness guarantee).
- `RdkitMorganError`: `Aromaticity(AromaticityError)` (kekulization/internal-invariant
  failure, wrapping `rdkit_parity.rs`'s own error), `UnsupportedBondOrder { bond_idx, order }`
  (a `BondOrder` with no real RDKit `Bond::BondType` counterpart -- only chematic's
  SMARTS-query-only variants, which cannot occur for SMILES-parsed input; confirmed via a
  programmatically-built `Molecule` test since it can't be reached via `parse()`),
  `InternalInvariantViolation { reason }`.
- One shared computation, not independent per-field loops: a single pass over RDKit's
  `includeRedundantEnvironments=false` ("default") lifecycle populates all four
  `RdkitMorganEcfp4` fields (`fingerprint`, `sparse_counts`, `raw_bit_info`,
  `folded_bit_info`) at once.

**Results (same 5,048-input M4-A0 corpus, fresh dump + comparison against the same RDKit
oracle rows, not re-derived from M4-A0's own numbers):**

| Metric | Result |
|---|---|
| Preprocessing succeeded (`status: "success"`) | 5,046/5,048 |
| Preprocessing failed (`rdkit_parity_kekulization_failed`, the same 2 pinned fixtures as M4-A0) | 2/5,048 |
| Full exact match (default-lifecycle raw pairs, sparse counts, folded on-bits, folded bitInfo) among the 5,046 successful rows | **5,046/5,046 (100%)** |
| Hermetic equivalence to `rdkit_morgan_raw_trace`'s already-oracle-validated `raw_identifier_default` output, same already-aromatized molecule | confirmed (unit test, 4 representative fixtures) |
| Non-regression: `ecfp_regression_snapshot` (10 existing entry points: `ecfp4`, `ecfp6`, `ecfp` chiral, `ecfp_with_bitinfo`, `morgan_fp_counts`, `ecfp4_rdkit_invariants`, `ecfp6_rdkit_invariants`, `ecfp4_rdkit_environment_experimental`, `ecfp6_rdkit_environment_experimental`, `ecfp_with_bitinfo_rdkit_environment_experimental`), full 5,048-input corpus, SHA-256 before/after (git-worktree baseline at the pre-Phase-B commit) | byte-identical, 0 change |
| Unsupported-bond-order path (`BondOrder::QueryAny`, programmatically built -- cannot arise from `parse()`) | explicit `Err(UnsupportedBondOrder)`, confirmed by test |
| Positive control: a silently reintroduced Hueckel fallback would be numerically detectable (Hueckel perceives `c1cc[nH+]cc1`'s ring as fully aromatic where the real path correctly errors) | confirmed by test |

**Performance vs. `ecfp4_rdkit_environment_experimental` baseline** (5 independent process
runs each, full 5,048-corpus, median wall time, `/usr/bin/time -l` for peak RSS -- not a
Criterion-registered benchmark, see `feedback_criterion_gate_pseudo_replication`):

| | Baseline | Candidate | Ratio |
|---|---|---|---|
| Median wall time (5 runs) | 4.862s | 9.734s | **2.00x** |
| Peak RSS | ~20.2 MB | ~20.4 MB | 1.01x |

The ~2x ratio is fully attributable, not an unexplained regression: the baseline reads
whatever aromatic flags are already on the input `Molecule` and never calls an aromaticity
engine, while the candidate performs its own kekulization + RDKit-parity aromaticity
perception on every call (a per-molecule breakdown shows preprocessing is 46-56% of the
candidate's time on aromatic-ring-heavy molecules like benzene/aspirin/a steroid-like
fused system, dropping to 19-20% on large acyclic alkanes with no rings to perceive).
Per the acceptance-gate policy (stop and explain, don't silently tune), this is reported
as measured rather than optimized against.

Any `BondOrder` this engine cannot map to a real RDKit `BondType` (verified against
`Bond.h` during M4-A0: SINGLE=1, DOUBLE=2, TRIPLE=3, QUADRUPLE=4, AROMATIC=12, DATIVE=17,
ZERO=21; only chematic's SMARTS-only `Query*` variants have no RDKit equivalent) is an
explicit `Err`, never an implicit/guessed mapping.

- **Files:** `crates/chematic-fp/src/rdkit_morgan_ecfp4.rs`,
  `crates/chematic-fp/examples/rdkit_morgan_ecfp4_dump.rs`,
  `crates/chematic-fp/examples/rdkit_morgan_ecfp4_benchmark.rs`,
  `scripts/ecfp_rdkit_morgan_ecfp4_parity.py`.
- **How to regenerate:** `cargo run -p chematic-fp --release --example rdkit_morgan_ecfp4_dump
  -- <SMILES.csv> <out.jsonl>` + `python scripts/ecfp_rdkit_morgan_ecfp4_parity.py --chematic
  <out.jsonl> --rdkit-oracle <gen_ecfp_rdkit_environment_oracle.py output>`. Self-test:
  `python scripts/ecfp_rdkit_morgan_ecfp4_parity.py --self-test`.

### IO-1: SMILES table file I/O (`.smi`/`.smiles`/`.csv`/`.tsv`/`.txt`)

New streaming `SmilesRecordReader`/`SmilesRecordWriter` in
`crates/chematic-mol/src/smiles_table.rs`, built from a source-cited audit of
RDKit's `SmilesMolSupplier`/`SmilesWriter` (RDKit commit
`8afba32ec539dcb2369bc84549d802aca3f7eb39`, the true resolution of tag
`Release_2026_03_4`) rather than guessed behavior. `chematic.SmilesMolSupplier`/
`chematic.SmilesWriter` (Python, pyo3) match RDKit's constructor signatures;
`rdkit_compat.py`'s own wrapper classes of the same names were rewritten to
delegate to these (previously a separate, non-streaming, whole-file pure-Python
implementation that used no Rust parser at all).

**Oracle methodology (deliberately avoids the chematic/RDKit canonical-SMILES
divergence, which is a known, separately-tracked, unrelated issue):** never
compares chematic-canonical vs. RDKit-canonical SMILES strings directly.
Instead: (1) exact string equality of extracted `name`/property values
(pure tokenization, no chemistry); (2) each tool's *own* self-consistency
against the fixture's known ground-truth SMILES, canonicalized only within
that same tool — proves each tokenizer extracted the right substring
without ever comparing the two canonicalizers against each other.

**Results (235 rows, 8 scenarios covering space/tab/comma delimiters,
header/no-header, name/no-name column, extra properties, quoted CSV, blank
lines, comments, malformed SMILES, isotopes, charges, disconnected
fragments, stereochemistry):**

| Metric | Result |
|---|---|
| Status parity (success vs. unparseable, per row) | 235/235 (100%) |
| Known-malformed rows correctly rejected by both tools | 5/5 |
| Name/property exact-match (excluding 2 documented divergences below) | 100% |
| Chematic self-consistency vs. known ground truth | 230/230 (100%) |
| RDKit self-consistency vs. known ground truth (non-gating) | 230/230 (100%) |

**Two deliberate, documented divergences from RDKit found via this oracle
(not bugs — see `smiles_table.rs`'s module doc comment for the full
citation):**
1. `name_column=None` (RDKit's `nameColumn=-1`): RDKit falls back to the
   physical line number as `_Name`; chematic's `MoleculeRecord::name` is
   simply empty. Judged a low-value RDKit implementation detail not worth
   reproducing.
2. **RDKit's `SmilesMolSupplier` has no CSV-quote-awareness at all** for its
   comma-delimiter mode — confirmed via this oracle, not merely inferred: a
   quoted field like `"has, a comma"` is split into extra raw columns by
   RDKit's literal comma-splitting. Chematic implements a real RFC 4180
   *subset* (quoted fields, doubled-quote escaping, no multi-line quoted
   fields) — a genuine improvement, not a matched behavior.

CXSMILES is not recognized in the SMILES column (parsed via
`chematic_smiles::parse`, which has no CXSMILES support) — this matches
RDKit's *own* default for `SmilesMolSupplier`, which explicitly disables
CXSMILES for this entry point too.

**Performance** (10,000-record synthetic corpus, ~2% deliberately malformed rows to also
measure invalid-record recovery throughput; 5 independent process runs, `/usr/bin/time -l`):

| | chematic (`SmilesRecordReader`) | RDKit (`SmilesMolSupplier`, Python) |
|---|---|---|
| Records/sec (median of 5 runs) | ~137,000 | ~4,200 |
| Peak RSS | ~2.3 MB | ~45 MB |
| Success/error split | 9,800/200 | 9,800/200 (identical) |

The ~30x throughput difference is Python-interpreter-call overhead, not a controlled
same-language comparison — reported as reference/informational only, per the
"performance is never traded for correctness, and cross-language numbers aren't a gate"
policy. Both tools agree exactly on which 200/10,000 rows are malformed.

**Adversarial/fuzz-style coverage:** no `cargo-fuzz`/libfuzzer harness exists anywhere in
this workspace yet, and introducing that toolchain for one text-tokenizer module was judged
disproportionate — instead, 9 deterministic adversarial unit tests (empty input, truncated
input mid-record and mid-quote, a line exceeding `max_line_bytes`, a 500KB property value
within the limit, invalid-UTF-8 byte handling, 5,000-column rows, a 3,000-atom SMILES field,
and a 2,000-iteration seeded random-mutation corpus) assert only "no panic, no hang, no OOM" —
never a specific output — since malformed input must degrade to a clean `Err`, never worse.

- **Files:** `crates/chematic-mol/src/smiles_table.rs`, `crates/chematic-mol/examples/smiles_table_dump.rs`,
  `crates/chematic-mol/examples/smiles_table_benchmark.rs`,
  `scripts/gen_smiles_table_fixtures.py`, `scripts/gen_rdkit_smiles_table_oracle.py`,
  `scripts/smiles_table_io_parity.py`.
- **Reference tool:** RDKit 2026.03.3.
- **How to regenerate:** `python scripts/gen_smiles_table_fixtures.py --out-dir <dir> --corpus
  <SMILES.csv> --manifest-out <manifest.json>` + `cargo run -p chematic-mol --release --example
  smiles_table_dump -- <manifest.json> <dir> <out.jsonl>` + `python
  scripts/gen_rdkit_smiles_table_oracle.py --manifest <manifest.json> --fixtures-dir <dir> --out
  <oracle.jsonl>` + `python scripts/smiles_table_io_parity.py --chematic <out.jsonl> --rdkit-oracle
  <oracle.jsonl> --manifest <manifest.json>`. Self-test:
  `python scripts/smiles_table_io_parity.py --self-test`. The generated fixture files themselves
  are not committed (regenerable byte-for-byte from the corpus + script); only the generator
  scripts and this summary are.

## Summary results

See [rdkit/README.md](rdkit/README.md) for per-descriptor breakdowns.

| Metric | Corpus | Agreement |
|--------|--------|-----------|
| HBA / HBD / ARC | 4,999 mol | **100%** |
| MW, HAC | 175 mol | **100%** |
| TPSA | 175 mol (drug-like) | **100%** (±0.1 Å²) |
| TPSA | 4,999 mol | 93.3% (±0.1 Å²) |
| LogP | 175 mol | ~99% (±0.3) |

## Methodology

- Reference values are generated by RDKit Python API (rdkit-sys ≥ 2024.x)
- chematic values are computed via `chematic.from_smiles(smi).descriptors()`
- Agreement = fraction of molecules within the stated tolerance
- TPSA uses Ertl (2000) SMARTS-based approach in both tools
- Scripts are deterministic and pinned to RDKit 2026.03.3

## How to run

```bash
# Fast regression on 175-mol in-repo corpus (no download required)
pip install chematic rdkit pandas
python scripts/rdkit_benchmark.py

# Large-scale agreement on 5k ChEMBL subset
python scripts/bench5k.py ~/Downloads/SMILES.csv
python scripts/bench5k.py ~/Downloads/SMILES.csv --detail   # show mismatches
```
