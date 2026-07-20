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
production API; `ecfp4()`/`ecfp6()`/`ecfp4_rdkit_invariants()`/
`ecfp4_rdkit_environment_experimental()` are all unchanged (production
snapshot verified byte-identical, see Results).

Ported directly from RDKit's real C++ source, commit
[`8afba32ec539dcb2369bc84549d802aca3f7eb39`](https://github.com/rdkit/rdkit/blob/8afba32ec539dcb2369bc84549d802aca3f7eb39/Code/GraphMol/Fingerprints/MorganGenerator.cpp)
(the true resolution of tag `Release_2026_03_4`, independently verified via
the GitHub tags API this session -- see `THIRD_PARTY_NOTICES.md` for the
attribution and the note on two other, imprecise SHAs already in this
project's history under the same tag label).

**Results (5,048-input set = 5,000-mol corpus + PR #120's 41 fixtures + PR
#123's 4 fixtures + `ecfp_rdkit_m4a0_hash_fixtures.csv`'s 3, production
Hueckel aromaticity applied before hashing):**

| Metric | Result |
|---|---|
| Radius-0 numeric exact match | 5,048/5,048 (100%) |
| Full numeric exact match (radius 0-2, representative selection, sparse counts, folded bits, bitInfo) | 4,989/5,048 (98.83%) |
| Residual (59/5,048), re-run with `apply_aromaticity_rdkit_parity_experimental` instead of production Hueckel aromaticity | 5,048/5,048 (100%), 0 new mismatches among the other 4,989, 2/5,048 fell back to Hueckel on kekulization failure |
| PR #123's 9 unique representative-selection residuals | all 9 resolve to `exact_match` under the RDKit-exact hash alone (no aromaticity-engine swap needed) |
| PR #123's Kekule-pyridine sparse-count-shape mismatch (`C1=CC=NC=C1`) | resolves to `exact_match`; confirms the documented root cause (FNV-1a-specific hash collision, not a suppression defect) |
| Production API byte-identical (`ecfp_regression_snapshot`, before/after SHA-256) | confirmed, 0 change |
| Oracle regeneration determinism (`--verify-determinism`, full 5,048 input) | byte-identical across two runs |
| Positive controls (radius-0/1 identifier, bond invariant, 32-bit-wrapping removal, representative swap, folded-bit, dropped row, duplicate row ID) | all 8 correctly cause non-zero exit; reverted, never committed |

The entire 59/5,048 (1.17%) residual under production Hueckel aromaticity
traces to ONE mechanism: chematic's Hueckel-based aromaticity *perception*
disagreeing with RDKit's own aromaticity model on specific fused/macrocyclic
ring systems (e.g. `C[Si](C)(C)c1ccc(C2=Cc3ccccc3C3=NCCCN23)cc1`) -- not a
hash defect. This project already has a dedicated RDKit-parity aromaticity
engine (`apply_aromaticity_rdkit_parity_experimental`,
`crates/chematic-perception/src/rdkit_parity.rs`) built for exactly this
kind of disagreement; feeding it instead of production `apply_aromaticity`
resolves the residual to 0 with no regressions anywhere else in the corpus.

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
  (RDKit-parity aromaticity engine, full corpus), `ecfp_rdkit_raw_identifier_parity_oracle_manifest.json`
- **Reference tool:** RDKit 2026.03.3 (`rdFingerprintGenerator.GetMorganGenerator`,
  `includeRedundantEnvironments` True and False variants; same pinned option
  set as every other Morgan-parity mode in this crate)
- **How to regenerate:** `python scripts/gen_ecfp_rdkit_environment_oracle.py`
  + `cargo run -p chematic-fp --release --features diagnostics --example
  rdkit_morgan_hash_dump` + `python scripts/ecfp_rdkit_raw_identifier_parity.py`
  (see each script's docstring for exact invocation). Residual-mechanism
  confirmation: `cargo run -p chematic-fp --release --features diagnostics
  --example rdkit_morgan_hash_dump_aromaticity_variant`.

**Scope for a future implementation PR** (not started): promote
`rdkit_morgan_hash.rs`'s reference engine to a real, RDKit-compatible raw
Morgan identifier API (`_rdkit_hash_experimental`-scoped, not touching
`ecfp4()`/`ecfp4_rdkit_invariants()`), covering raw ID generation,
representative selection, sparse-count-map, folded-bit, and bitInfo parity
-- gated on first resolving which aromaticity engine such an API should use
(production Hueckel, with its documented ~1% fused/macrocyclic gap, or the
RDKit-parity engine, which is currently `Result`-fallible and unproven
outside this corpus).

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
