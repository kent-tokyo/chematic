# Speed without output change — post-v1.0.23 performance branch

This record checks that the `perf/speed-3` branch changes no output relative
to v1.0.23, and measures how much faster it is, both in Rust (cold, per
molecule) and through the Python per-operation matrix against RDKit
2026.03.6. It is a source-level record on a shared 2-vCPU x86_64 VM, not a
published-package or release claim.

## Scope

| Item | Value |
|---|---|
| Base | v1.0.23, commit `6f8645efaaccc9fc14b9a63736edaadc3aafd515` |
| Branch | `20b2cbc3` (main change set) and `6bd3fae0` (pruning warm-up fix); timings below are for `6bd3fae0` unless marked |
| Wheels | `maturin build --release` of `crates/chematic-py` at each revision (SHA-256 in the JSON) |
| RDKit | 2026.03.6 (PyPI) |
| Python / host | 3.11.15, Linux x86_64 VM, 2 vCPU "Intel(R) Xeon(R) Processor @ 2.10GHz" |
| Output-identity script | [`scripts/perf_digest_diff.sh`](../scripts/perf_digest_diff.sh) with [`tools/perf_digest`](../tools/perf_digest) |
| Rust timing script | [`scripts/perf_digest_time_pair.py`](../scripts/perf_digest_time_pair.py) |
| Python timing script | [`scripts/bench_python_op_matrix_vs_rdkit.py`](../scripts/bench_python_op_matrix_vs_rdkit.py) |
| Agreement script | [`scripts/check_rdkit_output_agreement.py`](../scripts/check_rdkit_output_agreement.py) |

## 1. Output identity (v1.0.23 vs branch)

`perf_digest_diff.sh` builds `tools/perf_digest` against both trees and
compares every (corpus, row, operation) output byte for byte over the seven
benchmark corpora (46,736 molecules).

| Run | Operations | Rows compared | Differing rows | Raw data |
|---|---|---:|---:|---|
| `20b2cbc3`, all operations | 29 (SSSR rings, ring/aromatic ring counts, RDKit-parity aromatic view, TPSA, `rdkit_tpsa`, HBA/HBD, rotatable bonds, ring bundle, LogP, MR, QED, native and RDKit-compatible ECFP4 incl. bit info, canonical SMILES, stereocenters, InChI, CIP, SMARTS all/non-unique/existence, largest fragment, standardize, MACCS, PAINS, Brenk) | 1,402,080 | **0** | [JSON](2026-09-24-perf-digest-diff-v1.0.23-20b2cbc3.json) |
| `6bd3fae0`, SMARTS-affected operations | 12 (SMARTS all/non-unique/existence/first-embedding, MACCS, QED, PAINS, Brenk, LogP, MR, standardize, RDKit-parity aromatic view) | 560,832 | **0** | [JSON](2026-09-24-perf-digest-diff-v1.0.23-6bd3fae0-smarts.json) |

In the `20b2cbc3` run the operation list contained the name `rdkit_ecfp4`
twice (the fingerprint alone and a variant with sorted bit info), so that
name shows 93,472 rows; both variants were identical. The tool now calls the
second one `rdkit_ecfp4_bitinfo`. The JSON's `base_revision` is the v1.0.23
tag object `2ab467b8`, which points at commit `6f8645ef`; the script now
records the commit.

## 2. Rust cold timing (µs per molecule, ChEMBL 5k)

Median of 3 alternating rounds, each the best of 3 cold passes ([JSON](2026-09-24-perf-digest-time-pair-v1.0.23-6bd3fae0.json)). PAINS was left out of this run to stay within the session's time box.

| Operation (`tools/perf_digest` name) | v1.0.23 | branch | speedup |
|---|---:|---:|---:|
| SSSR rings (formatted) (`sssr`) | 16.63 | 12.17 | 1.37x |
| ring count (`ring_count`) | 15.21 | 1.66 | 9.16x |
| aromatic ring count (`aromatic_ring_count`) | 19.85 | 15.37 | 1.29x |
| RDKit-parity aromatic view (`rdkit_parity_view`) | 23.44 | 15.07 | 1.56x |
| TPSA (N,O,S,P) (`tpsa`) | 22.16 | 13.89 | 1.60x |
| `rdkit_tpsa` (`rdkit_tpsa`) | 21.98 | 13.85 | 1.59x |
| HBA (`hba`) | 16.54 | 2.22 | 7.45x |
| rotatable bonds (`rotb`) | 16.48 | 2.46 | 6.70x |
| Crippen LogP (`logp`) | 78.09 | 68.30 | 1.14x |
| QED (`qed`) | 508.79 | 423.80 | 1.20x |
| MACCS (`maccs`) | 655.28 | 531.07 | 1.23x |
| 15 SMARTS, `has_match_bounded` with rings (`has_sub`) | 47.43 | 38.58 | 1.23x |
| 15 SMARTS, first-embedding `find_matches_with_config` (`has_sub_first`) | 44.93 | 41.73 | 1.08x |

## 3. Python matrix against RDKit (µs per row, ChEMBL 5k)

Base [v1.0.23 matrix](2026-09-24-python-op-matrix-vs-rdkit-v1.0.23-6f8645ef.json),
branch [6bd3fae0 matrix](2026-09-24-python-op-matrix-vs-rdkit-branch-6bd3fae0.json)
(the superseded [20b2cbc3 matrix](2026-09-24-python-op-matrix-vs-rdkit-branch-20b2cbc3.json)
is kept because it exposed the single-query `has_substructure` regression
fixed in `6bd3fae0`). Only `checked` operations are listed (outputs compared
row by row with RDKit); "agree" is the branch's row agreement, which equals
the base's because the outputs are identical.

| Operation | v1.0.23 | branch | branch speedup | RDKit | v1.0.23 vs RDKit | branch vs RDKit | agree |
|---|---:|---:|---:|---:|---:|---:|---|
| `ring_count` | 16.98 | 2.24 | 7.59x | 0.76 | 0.04x | 0.34x | 4946/5000 |
| `hba` | 18.05 | 3.18 | 5.67x | 79.74 | 4.17x | 25.05x | 3641/5000 |
| `rotatable_bonds` | 18.42 | 3.52 | 5.24x | 98.99 | 5.11x | 28.16x | 5000/5000 |
| `parse+ring_count` | 20.67 | 5.42 | 3.81x | 116.30 | 5.40x | 21.45x | 4946/5000 |
| `torsion(rdkit-compatible)` | 31.37 | 19.21 | 1.63x | 120.54 | 3.85x | 6.27x | 4999/5000 |
| `rdkit_tpsa` | 25.15 | 15.63 | 1.61x | 31.93 | 1.19x | 2.04x | 5000/5000 |
| `tpsa` | 23.82 | 15.28 | 1.56x | 31.94 | 1.20x | 2.09x | 3659/5000 |
| `parse+tpsa` | 27.39 | 17.91 | 1.53x | 124.40 | 4.51x | 6.95x | 3659/5000 |
| `lipinski_bundle(mw,logp,hbd,hba)` | 96.75 | 69.31 | 1.40x | 386.83 | 3.72x | 5.58x | 3641/5000 |
| `has_substructure c1ccccc1` | 2.98 | 2.17 | 1.37x | 5.26 | 1.91x | 2.42x | 2000/2000 |
| `atom_pair(rdkit-compatible)` | 44.00 | 33.24 | 1.32x | 125.66 | 3.01x | 3.78x | 5000/5000 |
| `largest_fragment` | 5.02 | 3.81 | 1.32x | 80.32 | 15.67x | 21.05x | 2000/2000 |
| `logp` | 84.79 | 68.54 | 1.24x | 293.13 | 3.56x | 4.28x | 5000/5000 |
| `parse+aromatic_ring_count` | 23.75 | 19.39 | 1.23x | 120.16 | 4.86x | 6.20x | 4940/5000 |
| `morgan_r2_2048(rdkit-compatible)` | 50.07 | 40.99 | 1.22x | 54.19 | 1.11x | 1.32x | 5000/5000 |
| `mr` | 82.85 | 68.23 | 1.21x | 301.33 | 3.53x | 4.42x | 5000/5000 |
| `parse+logp` | 86.76 | 71.49 | 1.21x | 399.30 | 4.57x | 5.59x | 5000/5000 |
| `aromatic_ring_count` | 21.41 | 17.86 | 1.20x | 4.59 | 0.20x | 0.26x | 4940/5000 |
| `qed` | 498.75 | 421.32 | 1.18x | 1191.01 | 2.40x | 2.83x | 4962/5000 |
| `canonical_tautomer` | 1635.51 | 1393.96 | 1.17x | 31077.93 | 19.06x | 22.29x | 241/300 |
| `tanimoto_1xN(per target, compatible Morgan)` | 0.10 | 0.08 | 1.16x | 0.10 | 1.11x | 1.17x | 50/50 |
| `maccs` | 625.69 | 540.88 | 1.16x | 770.41 | 1.25x | 1.42x | 4995/5000 |
| `has_substructure *~*~*~*~*~*` | 1.71 | 1.53 | 1.12x | 3.70 | 2.22x | 2.42x | 2000/2000 |
| `has_substructure C(=O)N` | 1.84 | 1.69 | 1.09x | 4.69 | 2.60x | 2.77x | 2000/2000 |
| `murcko_scaffold` | 17.62 | 16.35 | 1.08x | 89.51 | 5.20x | 5.48x | 4998/5000 |
| `chi1v` | 5.97 | 5.58 | 1.07x | 4.86 | 0.76x | 0.87x | 3380/5000 |
| `has_substructure [OH]` | 1.73 | 1.63 | 1.06x | 3.98 | 2.32x | 2.44x | 2000/2000 |
| `find_matches [#6]~[#7]` | 3.74 | 3.53 | 1.06x | 25.69 | 6.74x | 7.28x | 2000/2000 |
| `standardize vs Cleanup+Uncharge+Canonicalize` | 1646.69 | 1583.23 | 1.04x | 29896.11 | 18.51x | 18.88x | 241/300 |
| `has_substructure [NX3;H2,H1;!$(NC=O)]` | 1.77 | 1.72 | 1.03x | 6.40 | 3.85x | 3.72x | 2000/2000 |
| `kappa2` | 3.43 | 3.43 | 1.00x | 74.37 | 21.17x | 21.68x | 7/5000 |
| `exact_mass` | 1.82 | 1.85 | 0.98x | 2.08 | 1.13x | 1.12x | 5000/5000 |
| `pattern_fp(rdkit-compatible)` | 257.29 | 264.99 | 0.97x | 301.28 | 1.22x | 1.14x | 2000/2000 |
| `mw` | 1.31 | 1.35 | 0.97x | 2.25 | 1.71x | 1.67x | 5000/5000 |
| `rdkit_fp(rdkit-compatible)` | 1229.70 | 1277.52 | 0.96x | 1437.89 | 1.16x | 1.13x | 2000/2000 |
| `formula` | 2.26 | 2.37 | 0.95x | 4.41 | 1.97x | 1.86x | 4992/5000 |
| `bertz_ct` | 1.24 | 1.31 | 0.95x | 345.56 | 278.53x | 264.52x | 0/5000 |
| `has_substructure [CX3](=O)[OX2H1]` | 1.84 | 1.96 | 0.94x | 4.68 | 2.47x | 2.39x | 2000/2000 |
| `num_stereocenters` | 34.57 | 38.46 | 0.90x | 255.16 | 6.65x | 6.63x | 4932/5000 |
| `has_substructure c1ccc2ccccc2c1` | 9.22 | 10.67 | 0.86x | 10.94 | 1.08x | 1.03x | 2000/2000 |
| `hbd` | 0.94 | 1.09 | 0.86x | 23.21 | 25.10x | 21.28x | 5000/5000 |
| `has_substructure [#7;R]` | 3.14 | 3.68 | 0.85x | 3.75 | 1.21x | 1.02x | 2000/2000 |
| `fsp3` | 1.03 | 1.23 | 0.84x | 2.31 | 2.35x | 1.89x | 5000/5000 |
| `neutralize` | 4.40 | 5.28 | 0.83x | 41.61 | 13.09x | 7.89x | 2000/2000 |

Summary: both runs have 44 counted operations, 41 of them faster than RDKit
and 25 faster with full row agreement — the same 25 operations. The branch
widens existing margins (ring/HBA/rotatable-bond counts, TPSA, RDKit-compatible
atom-pair/torsion/Morgan, LogP/MR, QED, MACCS) rather than flipping a checked
operation. On this shared 2-vCPU VM, run-to-run noise is roughly ±15% for
microsecond-scale rows: the unchanged `hbd`, `fsp3` and `neutralize` paths
move by that much, and `has_substructure [#7;R]` measured 2.60 µs in the
`20b2cbc3` run and 3.68 µs here with no code change on its path in between.
The two ring-query rows (`c1ccc2ccccc2c1`, `[#7;R]`) stay at parity with RDKit
(1.03x, 1.02x).

## 4. RDKit output agreement

[`check_rdkit_output_agreement.py`](../scripts/check_rdkit_output_agreement.py)
on the `6bd3fae0` wheel ([JSON](2026-09-24-rdkit-agreement-branch-6bd3fae0.json))
gives, for all 15 operations on all three corpora (ChEMBL 5k, NCI 5k, descriptor
census 5k), exactly the counts of the previous
[accuracy-branch record](2026-09-24-rdkit-agreement-accuracy-branch.md)
(`ea967d47`, whose code v1.0.23 released) — e.g. ChEMBL: RDKit-compatible Morgan
and atom pair 5000/5000, MACCS 4995, QED 4962, Murcko 4998. The speed work kept
accuracy unchanged.

## What changed

- **SSSR** (`chematic-perception`): the cyclic-subgraph bridge flags use a
  flat adjacency; Horton candidates test path simplicity and build the bond
  mask in O(1) from per-root path masks; roots of ring systems that are a
  single cycle are skipped after the first (they can only rediscover the same
  cycle); GF(2) pivots live in an array; the canonical tie-break ranks are
  computed only when a same-length candidate group must be sorted, and only
  for atoms within eight bonds of a ring. `find_sssr_horton_reference` stays
  the differential oracle (new unit test).
- **Counts without rings**: `ring_count` is the cycle rank
  (`sssr_ring_count`), and HBA / rotatable bonds take ring bonds from the
  linear bridge pass (an SSSR is a cycle basis, so its ring bonds are exactly
  the non-bridge bonds).
- **RDKit-parity aromatic view**: no redundant molecule clones, the
  re-perception verdict is reused, and the Kekulé-form SSSR/perception pass is
  skipped when a bridge test on the aromaticity-candidate subgraph proves the
  verdict cannot extend the explicit aromatic input (new unit test checks the
  skipped cases are subsets).
- **SMARTS**: `has_match_with_config` runs the existence search without
  result maps; recursive `$(...)` sub-queries reuse per-search plans and
  buffers; exact pruning of acyclic target atoms for ring query atoms (after a
  256-check warm-up, and only for queries of at most 7 atoms in the
  map-returning search so result-map layout cannot change); an element-count
  screen. Pruning and the screen are disabled when a visit budget is set.

## Remaining gaps and notes

- `ring_count` and `aromatic_ring_count` on prepared inputs still lose to
  RDKit, which perceives rings and aromaticity inside `MolFromSmiles`; the
  parse-inclusive rows are 20x and 6x faster.
- `embed+minimize_mmff94` and `standardize vs Cleanup` are not equivalent
  operations and are unchanged here.
- `tanimoto_1xN` is bounded by the popcount: the wheel is built for baseline
  x86-64 (no POPCNT). A local micro-benchmark shows ~2x on the 2048-bit
  kernel with `-C target-cpu=x86-64-v2`; that is a distribution decision and
  was not changed.
