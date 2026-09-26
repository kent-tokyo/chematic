# Speed without output change — post-v1.0.25 performance branch (`perf/speed-4`)

This record checks that the `perf/speed-4` branch changes no output relative
to its base (`fix/rebaseline-issues-v1025`, commit `912c4b7d`), and measures
how much faster it is, both in Rust (cold, per molecule) and through the
Python per-operation matrix against RDKit 2026.03.6. The branch targets the
operations that were slower than or at parity with RDKit in the base matrix.
It is a source-level record on a shared 2-vCPU x86_64 VM, not a
published-package or release claim.

## Scope

| Item | Value |
|---|---|
| Base | `912c4b7d81ce7735af13b007b6289c555cde9a8f` (v1.0.25 + rebaseline evidence) |
| Branch | `perf/speed-4` at `5269d69f5716122cf00bc310c64618927532100d` (sections 1–4); follow-up commits up to `98f5774645fa8596160c54ff33463171d4916a12` in section 5 |
| Wheels | `maturin build --release` of `crates/chematic-py` at each revision (SHA-256 in the JSON) |
| RDKit | 2026.03.6 (PyPI) |
| Python / host | 3.11.15, Linux x86_64 VM, 2 vCPU "Intel(R) Xeon(R) Processor @ 2.10GHz" |
| Output-identity script | [`scripts/perf_digest_diff.sh`](../scripts/perf_digest_diff.sh) with [`tools/perf_digest`](../tools/perf_digest) (`HEAD_FEATURES=candidate-apis`) |
| Rust timing script | [`scripts/perf_digest_time_pair.py`](../scripts/perf_digest_time_pair.py) |
| Python timing script | [`scripts/bench_python_op_matrix_vs_rdkit.py`](../scripts/bench_python_op_matrix_vs_rdkit.py) |
| Agreement script | [`scripts/check_rdkit_output_agreement.py`](../scripts/check_rdkit_output_agreement.py) |

## 1. Output identity (base vs branch)

`perf_digest_diff.sh` builds `tools/perf_digest` against both trees and
compares every (corpus, row, operation) output byte for byte over the seven
benchmark corpora (46,736 molecules). The head harness is built with the
`candidate-apis` feature: ops such as `rdkit_ecfp4_bits`,
`find_matches_perceived` or `rdkit_parity_ok` call the new borrowing /
bit-only entry points on the branch and the base's equivalent public path on
the base, so the comparison checks the new APIs against the old results.

| Run | Operations | Rows compared | Differing rows | Raw data |
|---|---|---:|---:|---|
| `5269d69f`, all operations | 49 (SSSR rings and ring sets, ring/aromatic ring counts incl. RDKit-compatible, RDKit-parity aromatic view (full and ok/err), kekulize, TPSA, `rdkit_tpsa`, HBA/HBD, rotatable bonds, ring bundle, LogP, MR, QED, chi (each, 1v, all), kappa, native and RDKit-compatible ECFP4 incl. bit info/bits/prepared, RDKit fingerprint, pattern fingerprint, atom pair, torsion, MACCS, canonical SMILES, stereocenters, InChI, CIP, SMARTS all/non-unique/existence/perceived/atom sets, largest fragment, standardize, PAINS, Brenk) | 2,290,064 | **8** (all `kekulize` error text, see below) | [JSON](2026-09-25-perf-digest-diff-912c4b7d-5269d69f.json) |

The 8 differing rows are all `kekulize` on molecules that cannot be
kekulized: both revisions return the same error, but the base names the first
unmatched atom it meets while iterating a `std::collections::HashSet`, whose
order is randomized per process, so the base's own message changes from run
to run (a second base-vs-branch run over the same corpora gave 0 differing
rows). The branch reports the lowest-index unmatched atom, which is
deterministic. Successful kekulizations (every bond order) are identical.

## 2. Rust cold timing (µs per molecule, ChEMBL 5k)

Median of 3 alternating rounds, each the best of 3 cold passes (fresh
molecule clones, empty derived caches) ([JSON](2026-09-25-perf-digest-time-pair-912c4b7d-5269d69f.json)).

| Operation (`tools/perf_digest` name) | base | branch | speedup |
|---|---:|---:|---:|
| RDKit-parity aromatic view (`rdkit_parity_ok`) | 11.38 | 3.59 | 3.17x |
| aromatic ring count (`aromatic_ring_count`) | 15.66 | 7.62 | 2.06x |
| RDKit-compatible aromatic ring count incl. view (`rdkit_aromatic_ring_count`) | 34.72 | 21.35 | 1.63x |
| ring count (`ring_count`) | 1.70 | 0.34 | 5.00x |
| TPSA (`tpsa`) | 13.51 | 3.68 | 3.67x |
| 15 SMARTS, perceived-aromaticity existence (`has_sub_perceived`) | 53.31 | 37.43 | 1.42x |
| 1 SMARTS `[OH]`, perceived, one query per fresh molecule (`has_sub_perceived_1`) | 11.97 | 3.26 | 3.67x |
| 15 SMARTS, perceived atom sets (`find_matches_perceived`) | 163.72 | 152.32 | 1.07x |
| RDKit fingerprint (`rdkit_rdk_fp`) | 1095.91 | 186.89 | 5.86x |
| RDKit pattern fingerprint (`rdkit_pattern_fp`) | 242.84 | 179.80 | 1.35x |
| RDKit-compatible Morgan r2/2048 bits (`rdkit_ecfp4_bits`) | 31.99 | 13.10 | 2.44x |
| RDKit-compatible atom pair (`rdkit_atom_pair`) | 33.70 | 21.74 | 1.55x |
| RDKit-compatible torsion (`rdkit_torsion`) | 20.18 | 9.13 | 2.21x |
| chi1v (`chi1v`) | 3.92 | 1.97 | 1.99x |
| chi0–chi4, plain and valence (`chi_all`) | 62.00 | 26.21 | 2.37x |
| MACCS (`maccs`) | 534.24 | 508.87 | 1.05x |
| QED (`qed`) | 416.06 | 380.05 | 1.09x |
| kekulize (`kekulize`) | 6.43 | 2.54 | 2.53x |

## 3. Python matrix against RDKit (µs per row, ChEMBL 5k)

Base [912c4b7d matrix](2026-09-25-python-op-matrix-vs-rdkit-912c4b7d.json),
branch [5269d69f matrix](2026-09-25-python-op-matrix-vs-rdkit-branch-5269d69f.json).
Only `checked` operations are listed (outputs compared row by row with
RDKit); "agree" is the branch's row agreement, which equals the base's for
every operation.

| Operation | base | branch | branch speedup | RDKit | base vs RDKit | branch vs RDKit | agree |
|---|---:|---:|---:|---:|---:|---:|---|
| `has_substructure *~*~*~*~*~*` | 14.84 | 1.24 | 11.97x | 3.99 | 0.19x | 3.22x | 2000/2000 |
| `has_substructure [#7;R]` | 14.86 | 1.88 | 7.92x | 3.79 | 0.22x | 2.02x | 2000/2000 |
| `find_matches [#6]~[#7]` | 17.28 | 2.75 | 6.27x | 21.78 | 1.14x | 7.91x | 2000/2000 |
| `rdkit_fp(rdkit-compatible)` | 1229.78 | 202.12 | 6.08x | 1380.05 | 1.10x | 6.83x | 2000/2000 |
| `has_substructure [CX3](=O)[OX2H1]` | 18.50 | 4.97 | 3.73x | 4.42 | 0.19x | 0.89x | 2000/2000 |
| `has_substructure [OH]` | 14.21 | 4.06 | 3.50x | 3.96 | 0.23x | 0.97x | 2000/2000 |
| `morgan_r2_2048(rdkit-compatible)` | 40.00 | 12.18 | 3.28x | 54.58 | 1.24x | 4.48x | 5000/5000 |
| `has_substructure C(=O)N` | 14.51 | 4.52 | 3.21x | 4.75 | 0.23x | 1.05x | 2000/2000 |
| `tpsa` | 13.73 | 4.63 | 2.96x | 29.66 | 1.69x | 6.40x | 3659/5000 |
| `has_substructure c1ccccc1` | 15.12 | 5.30 | 2.85x | 5.17 | 0.29x | 0.98x | 2000/2000 |
| `has_substructure [NX3;H2,H1;!$(NC=O)]` | 13.94 | 4.89 | 2.85x | 6.38 | 0.37x | 1.30x | 2000/2000 |
| `ring_count` | 1.81 | 0.72 | 2.51x | 0.77 | 0.29x | 1.06x | 4946/5000 |
| `parse+tpsa` | 17.55 | 7.57 | 2.32x | 118.49 | 7.00x | 15.65x | 3659/5000 |
| `rdkit_tpsa` | 13.46 | 6.52 | 2.06x | 31.77 | 1.64x | 4.87x | 5000/5000 |
| `torsion(rdkit-compatible)` | 18.24 | 8.84 | 2.06x | 106.00 | 6.35x | 12.00x | 4999/5000 |
| `aromatic_ring_count` | 16.81 | 8.64 | 1.95x | 4.29 | 0.21x | 0.50x | 4940/5000 |
| `chi1v` | 4.33 | 2.88 | 1.50x | 5.15 | 0.99x | 1.79x | 3380/5000 |
| `atom_pair(rdkit-compatible)` | 31.50 | 21.73 | 1.45x | 125.87 | 3.92x | 5.79x | 5000/5000 |
| `pattern_fp(rdkit-compatible)` | 264.96 | 184.20 | 1.44x | 314.60 | 1.10x | 1.71x | 2000/2000 |
| `parse+ring_count` | 4.98 | 3.62 | 1.38x | 117.07 | 22.37x | 32.38x | 4946/5000 |
| `lipinski_bundle(mw,logp,hbd,hba)` | 78.26 | 58.17 | 1.35x | 322.99 | 3.82x | 5.55x | 3641/5000 |
| `has_substructure c1ccc2ccccc2c1` | 18.43 | 14.42 | 1.28x | 11.46 | 0.51x | 0.79x | 2000/2000 |
| `parse+aromatic_ring_count` | 18.37 | 14.49 | 1.27x | 124.59 | 6.24x | 8.60x | 4940/5000 |
| `parse+logp` | 70.84 | 57.50 | 1.23x | 417.31 | 5.76x | 7.26x | 5000/5000 |
| `logp` | 68.44 | 56.14 | 1.22x | 287.83 | 4.10x | 5.13x | 5000/5000 |
| `mr` | 64.21 | 55.15 | 1.16x | 291.47 | 4.09x | 5.29x | 5000/5000 |
| `hba` | 2.99 | 2.63 | 1.13x | 50.93 | 13.80x | 19.36x | 3641/5000 |
| `largest_fragment` | 4.06 | 3.77 | 1.08x | 79.90 | 18.19x | 21.21x | 2000/2000 |
| `maccs` | 519.18 | 487.86 | 1.06x | 771.43 | 1.50x | 1.58x | 4995/5000 |
| `qed` | 395.14 | 372.32 | 1.06x | 1146.78 | 2.91x | 3.08x | 4962/5000 |
| `num_stereocenters` | 36.40 | 35.10 | 1.04x | 232.66 | 5.91x | 6.63x | 4932/5000 |
| `rotatable_bonds` | 2.81 | 2.72 | 1.03x | 101.57 | 31.30x | 37.32x | 5000/5000 |
| `standardize vs Cleanup+Uncharge+Canonicalize` | 1398.80 | 1377.29 | 1.02x | 31407.14 | 20.29x | 22.80x | 241/300 |
| `neutralize` | 4.40 | 4.40 | 1.00x | 60.86 | 11.57x | 13.82x | 2000/2000 |
| `tanimoto_1xN(per target, compatible Morgan)` | 0.08 | 0.08 | 0.98x | 0.09 | 1.09x | 1.11x | 50/50 |
| `murcko_scaffold` | 15.82 | 16.27 | 0.97x | 87.66 | 5.39x | 5.39x | 4998/5000 |
| `mw` | 1.15 | 1.19 | 0.97x | 2.37 | 1.49x | 2.00x | 5000/5000 |
| `canonical_tautomer` | 1229.75 | 1335.90 | 0.92x | 30092.20 | 22.91x | 22.53x | 241/300 |
| `formula` | 1.87 | 2.18 | 0.86x | 6.17 | 1.97x | 2.83x | 4992/5000 |
| `kappa2` | 3.08 | 3.63 | 0.85x | 55.82 | 18.03x | 15.40x | 7/5000 |
| `bertz_ct` | 1.13 | 1.39 | 0.81x | 335.11 | 294.63x | 241.58x | 0/5000 |
| `exact_mass` | 0.91 | 1.29 | 0.70x | 2.16 | 1.65x | 1.67x | 5000/5000 |
| `fsp3` | 0.65 | 1.01 | 0.64x | 2.20 | 2.56x | 2.18x | 5000/5000 |
| `hbd` | 0.56 | 1.09 | 0.51x | 25.57 | 32.80x | 23.49x | 5000/5000 |

Summary: 44 counted operations in both runs. Faster than RDKit: base 33,
branch 39; faster with full row agreement: base 17, branch 21. The base's
`has_substructure` rows (0.19x–0.51x of RDKit: every call built the
RDKit-parity aromatic view) are now at parity or faster except
`[CX3](=O)[OX2H1]` (0.89x) and the naphthalene query (0.79x); `ring_count`
flips to 1.06x and `chi1v` to 1.79x.

The rows at the bottom (`hbd`, `fsp3`, `exact_mass`, `bertz_ct`, `kappa2`,
`formula`, `mw`) are sub-4 µs calls whose code paths the branch does not
touch (`kappa2` and `bertz_ct` share no code with the chi change); they are
dominated by Python call overhead and move by this much between runs on this
shared VM (the v1.0.23 record measured `hbd` at 0.94 µs and 1.09 µs for the
same code). The Rust cold timings above are the per-change evidence.

## 4. RDKit output agreement

[`check_rdkit_output_agreement.py`](../scripts/check_rdkit_output_agreement.py)
on the base and branch wheels (base [JSON](2026-09-25-rdkit-agreement-912c4b7d.json),
branch [JSON](2026-09-25-rdkit-agreement-branch-5269d69f.json)): for all 15 operations on all three corpora
(ChEMBL 5k, NCI 5k, descriptor census 5k) the branch's counts are identical
to the base's — e.g. ChEMBL: RDKit-compatible Morgan and atom pair 5000/5000,
torsion 4999, MACCS 4995, QED 4962, Murcko 4998, PAINS 5000, Brenk 2565. The
speed work kept accuracy unchanged. (The base wheel was rebuilt for this run;
its SHA-256 differs from the wheel of the base matrix because the build is not
bit-reproducible, the source revision is the same.)

## 5. Follow-up: `98f57746` (small-component ring bases, SMARTS degree pruning)

Three further commits on the same branch:

- `c8f4bf7a`: cyclic components of cycle rank 1–2 whose minimum cycle basis
  is unique (one ring; a spiro pair; a theta graph whose shortest path is
  strictly shortest) get their rings straight from the graph instead of the
  Horton selection, in the ring-set helpers and in the aromatic ring count
  (a unique minimum basis is what every SSSR selection returns). VF2 prunes
  target atoms with fewer neighbours than the query atom has query bonds
  (exact; not under a visit budget, and in the map-returning search only for
  queries of at most 7 atoms, like the existing cycle pruning). Search plans
  use incremental counts (same order), and `has_match_perceived` runs the
  size/element screen before building the aromatic view.
- `d41e1957`: the Python SMARTS pattern cache uses Fx hashing.
- `98f57746`: the aromatic ring count returns 0 for acyclic Kekulé input
  without running perception, and seeds the perceived copy's ring data.

**Output identity** ([JSON](2026-09-25-perf-digest-diff-912c4b7d-98f57746.json)):
50 operations x 46,736 molecules against `912c4b7d`, 2,336,800 rows; the only
differing rows (9) are again `kekulize` error text on molecules that cannot be
kekulized (the base's reported atom follows hash order, so the count moves
between runs: 8 in section 1, 9 here). RDKit agreement on the three corpora:
identical to the base for all 15 operations
([JSON](2026-09-25-rdkit-agreement-branch-98f57746.json)).

**Rust cold timing** (µs per molecule, ChEMBL 5k, median of 3 alternating
rounds of best-of-3; [JSON](2026-09-25-perf-digest-time-pair-912c4b7d-98f57746.json)):

| Operation | base `912c4b7d` | `5269d69f` (section 2) | `98f57746` | speedup vs base |
|---|---:|---:|---:|---:|
| aromatic ring count (`aromatic_ring_count`) | 15.93 | 7.62 | 6.08 | 2.62x |
| 15 SMARTS, perceived existence (`has_sub_perceived`) | 52.02 | 37.43 | 34.40 | 1.51x |
| 1 SMARTS `c1ccccc1`, perceived (`has_sub_perceived_q`) | 13.75 | – | 4.57 | 3.01x |
| 1 SMARTS `[OH]`, perceived (`has_sub_perceived_1`) | 10.43 | 3.26 | 3.29 | 3.17x |
| 15 SMARTS, `has_match_bounded` with rings (`has_sub`) | 39.13 | – | 29.99 | 1.30x |
| RDKit pattern fingerprint (`rdkit_pattern_fp`) | 249.46 | 179.80 | 131.09 | 1.90x |
| MACCS (`maccs`) | 535.03 | 508.87 | 445.85 | 1.20x |
| QED (`qed`) | 408.25 | 380.05 | 285.36 | 1.43x |
| RDKit-parity aromatic view (`rdkit_parity_ok`) | 11.54 | 3.59 | 3.96 | 2.91x |

(`rdkit_parity_ok` is unchanged by these commits; 3.59 vs 3.96 is run-to-run
noise.) By instruction count (callgrind, 2,000 ChEMBL molecules) the
naphthalene query's matching fell from ~89k to ~31k instructions per molecule
and the aromatic ring count from ~71k to ~57k.

**Python matrix** ([JSON](2026-09-25-python-op-matrix-vs-rdkit-branch-98f57746.json);
all 44 checked operations keep the base's row agreement). The clearest
change is `has_substructure c1ccc2ccccc2c1`: 14.4 µs → 7.3 µs (RDKit 10.3,
0.79x → 1.40x). The other single-query `has_substructure` rows (4–6 µs per
call) sit at parity with RDKit and change sign between runs on this VM: in
an A/B of the `5269d69f` and `98f57746` wheels, two alternating rounds each,
`[OH]` measured 5.04/5.60 and 3.28/5.09 µs for the same code paths, against
RDKit's 3.6–3.8. The Rust timings above are the per-change evidence.
`aromatic_ring_count` on prepared input stays at about 0.45–0.57x of RDKit
(RDKit counts rings perceived inside `MolFromSmiles`).

## What changed

- **RDKit-parity aromatic view** (`chematic-perception`): the shortcut
  decision (identity / clear flags / full path) is memoized per molecule and
  computed from shared per-component facts of the cyclic subgraph (ring-bond
  flags and 2-edge-connected component labels come from one bridge DFS).
  Explicit aromatic input is proven unchanged either without kekulization
  (exact components, rank-1 Hückel test, mixed-component bridge test) or by
  computing the RDKit-parity verdict of its Kekulé form directly on the
  explicit graph: each flagged atom's donor type and candidacy follow from its
  matched status, and the SSSR ring set is used only when it is decidable
  without canonical ranks. Molecules where the bridgehead-N fallback,
  exocyclic aromatic bonds or rank-dependent ring choices could matter still
  take the full path, which now kekulizes one copy instead of rebuilding
  twice. Callers borrow the view (`with_rdkit_parity_view`) instead of
  cloning it.
- **Rings**: component-restricted, order-free SSSR ring sets for the Hückel
  verdict and aromatic ring counts (ranks skipped when a whole same-length
  group is independent or dependent); `sssr_ring_count` in one union-find
  pass.
- **Kekulization**: index-based augmenting-path matcher, same matching as
  before; deterministic error text (see section 1).
- **SMARTS**: queries whose atoms and bonds cannot see aromaticity skip the
  aromatic view; molecules whose view is the identity are matched in place;
  element-count pre-filter in the VF2 search; `Mol.find_matches` builds the
  sorted atom sets straight from embeddings.
- **Fingerprints**: RDKit fingerprint walks subgraphs as a stream (the
  previous implementation stays as a test oracle); pattern fingerprint
  visits embeddings without result maps; RDKit-compatible Morgan emits bits
  directly (u128 neighbourhood masks for small molecules); atom pair, torsion
  and MACCS borrow the shared view.
- **Chi**: per-atom deltas computed once; the DFS order and summation are
  unchanged.

## Remaining gaps and notes

- `aromatic_ring_count` on prepared inputs (about 0.5x) and `ring_count`
  (parity) compare against RDKit values computed inside `MolFromSmiles`; the
  parse-inclusive rows are about 11x and 30x faster.
- Single-query `has_substructure` rows other than the naphthalene query are at
  parity with RDKit (within this VM's noise) after section 5. They are
  bounded by deciding the aromatic view, about 20k instructions (~2.5 µs) per
  fresh molecule: shared cyclic-component facts (~8.5k), the rank-1 / bridge
  proofs (~10k for a third of the molecules) and, for ~5% of molecules, the
  Kekulé-form verdict (~150k, mostly the Horton SSSR of large fused
  systems). Queries that cannot see aromaticity skip it and run 1.5–3.5x
  faster than RDKit.
- `tanimoto_1xN` stays at parity: the popcount kernel needs POPCNT, which the
  baseline x86-64 wheel does not enable and `#![forbid(unsafe_code)]` rules
  out runtime dispatch.
