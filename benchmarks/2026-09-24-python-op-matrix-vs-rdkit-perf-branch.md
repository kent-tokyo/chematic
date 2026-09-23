# Python per-operation timing vs RDKit — perception/SMARTS performance branch

This record measures the `perf/rdkit-speed-sssr-cache` source branch against
its base (`main` at `01a86b62`, v1.0.21) and against RDKit, one operation at a
time, from Python. It also records the differential output check that backs
the branch's "outputs unchanged" statement. It is a source-level
measurement on a small shared cloud VM, **not** a published-package or
release claim, and it supersedes the unrecorded "54/62 operations" figure
that was reported in conversation before the script and raw data were
under version control.

## Scope

| Item | Value |
|---|---|
| Base source | `01a86b621f401d22fc9a0e90f85e4c85b6dd49fe` (`main`, v1.0.21) |
| Branch source | `8cc01365c3e3829478609fe5e7bab4b7d947906b` (crates identical at `e1cfa3d8`, which adds the harnesses) |
| Wheels | `maturin build --release` of `crates/chematic-py` at each revision; SHA-256 base `cb4885df…f776d7`, branch `235a3477…3c5d7d` (full hashes in the JSON) |
| RDKit | 2026.03.6 (PyPI wheel) |
| Python | 3.11.15, CPython |
| Host | Linux x86_64 VM, 2 vCPU "Intel(R) Xeon(R) Processor @ 2.10GHz"; shared cloud host |
| Corpus | `scripts/chembl_accuracy_corpus_4999.smi`, SHA-256 `1c47371dcbe37f4e0a141bf545b72bf238de2761fa3894fa251a552d84728d3e`; 5,000/5,000 rows parsed by both libraries |
| Script | [`scripts/bench_python_op_matrix_vs_rdkit.py`](../scripts/bench_python_op_matrix_vs_rdkit.py) |
| Raw data | [base JSON](2026-09-24-python-op-matrix-vs-rdkit-base-01a86b62.json), [branch JSON](2026-09-24-python-op-matrix-vs-rdkit-branch-8cc01365.json) |

## Protocol

- Every operation is timed in both libraries in the same process,
  single-threaded: all chematic repeats, then all RDKit repeats.
- 3 repeats; the reported value is the median of (repeat wall time / rows).
  Every repeat time is kept in the JSON.
- Inputs are rebuilt before each repeat and construction is not timed:
  chematic re-parses the SMILES; RDKit re-parses with `MolFromSmiles`
  (`Chem.Mol(m)` would carry over cached Crippen/TPSA values).
- RDKit perceives rings and aromaticity inside `MolFromSmiles`; chematic does
  it lazily. Prepared-input rows therefore charge ring/aromaticity perception
  to chematic only; `parse+…` rows include parsing on both sides.
- Exceptions are counted per repeat and never timed as wins (all rows below
  had zero failures on both sides).
- The base and branch runs were made back to back on the same VM with no
  other job running; ratios below ~1.1x are within this host's run-to-run
  noise.

### Equivalence classes

| Class | Meaning | Counted? |
|---|---|---|
| `checked` | Outputs compared row by row in an untimed pass; agreement recorded | Yes, and reported separately with and without full agreement |
| `same-task` | Same task, legitimately different representation (e.g. canonical spelling) | No |
| `not-equivalent` | Different algorithm, scope or definition | Never |

## Result summary

| | Base `01a86b62` | Branch `8cc01365` |
|---|---:|---:|
| Checked operations | 43 | 43 |
| Checked and faster than RDKit (median ratio > 1) | 14 | **39** |
| …with full output agreement | 6 | **21** |
| …with full agreement and every chematic repeat faster than every RDKit repeat | 6 | **18** |

Output agreement with RDKit is identical in the base and branch runs for
every operation (the branch does not change outputs; see the differential
check below). Where agreement is below 100% the gap is pre-existing
behaviour of chematic's default API relative to that RDKit call (for example
TPSA includes S/P by default, HBA/QED/Murcko definitions differ), so those
rows are shown but not claimed as equivalent wins.

## Checked operations

Ratio = RDKit time / chematic time (>1 means chematic faster). µs per row.

| Operation | Rows | main chematic µs | main ratio | branch chematic µs | branch ratio | RDKit µs (branch run) | Output agreement |
|---|---:|---:|---:|---:|---:|---:|---|
| `mw` | 5000 | 1.02 | 1.67x | 1.13 | **1.49x** | 1.68 | 5000/5000 (full) |
| `exact_mass` | 5000 | 1.17 | 1.40x | 0.94 | **1.68x** | 1.57 | 5000/5000 (full) |
| `logp` | 5000 | 838.27 | 0.32x | 85.32 | **3.06x** | 260.86 | 5000/5000 (full) |
| `mr` | 5000 | 820.29 | 0.33x | 87.78 | **2.86x** | 250.97 | 5000/5000 (full) |
| `tpsa` | 5000 | 112.82 | 0.22x | 25.59 | **0.96x** | 24.54 | 3659/5000 |
| `hbd` | 5000 | 0.83 | 9.91x | 0.69 | **7.87x** | 5.41 | 5000/5000 (full) |
| `hba` | 5000 | 261.31 | 0.14x | 19.19 | **3.17x** | 60.84 | 3641/5000 |
| `rotatable_bonds` | 5000 | 246.42 | 0.35x | 20.35 | **4.19x** | 85.34 | 5000/5000 (full) |
| `fsp3` | 5000 | 0.86 | 1.85x | 0.79 | **2.00x** | 1.58 | 5000/5000 (full) |
| `ring_count` | 5000 | 242.17 | 0.00x | 20.72 | **0.02x** | 0.48 | 4946/5000 |
| `aromatic_ring_count` | 5000 | 255.17 | 0.01x | 24.13 | **0.14x** | 3.34 | 4940/5000 |
| `qed` | 5000 | 3993.74 | 0.28x | 459.16 | **2.36x** | 1084.73 | 3201/5000 |
| `kappa2` | 5000 | 3.07 | 23.34x | 3.25 | **21.26x** | 69.05 | 7/5000 |
| `chi1v` | 5000 | 4.89 | 4.40x | 6.15 | **3.70x** | 22.72 | 3380/5000 |
| `bertz_ct` | 5000 | 1.03 | 305.62x | 1.12 | **296.38x** | 333.36 | 0/5000 |
| `formula` | 5000 | 1.86 | 1.90x | 2.11 | **1.57x** | 3.31 | 4992/5000 |
| `num_stereocenters` | 5000 | 161.82 | 1.23x | 34.57 | **5.35x** | 184.82 | 4932/5000 |
| `lipinski_bundle(mw,logp,hbd,hba)` | 5000 | 1087.91 | 0.29x | 95.52 | **3.17x** | 302.80 | 3641/5000 |
| `parse+logp` | 5000 | 848.74 | 0.44x | 92.39 | **3.98x** | 367.35 | 5000/5000 (full) |
| `parse+tpsa` | 5000 | 121.19 | 0.98x | 24.92 | **4.14x** | 103.19 | 3659/5000 |
| `parse+ring_count` | 5000 | 252.72 | 0.42x | 19.67 | **5.60x** | 110.06 | 4946/5000 |
| `parse+aromatic_ring_count` | 5000 | 265.84 | 0.45x | 26.89 | **3.79x** | 101.95 | 4940/5000 |
| `morgan_r2_2048(rdkit-compatible)` | 5000 | 143.10 | 0.37x | 51.81 | **0.97x** | 50.05 | 5000/5000 (full) |
| `maccs` | 5000 | 41458.32 | 0.02x | 512.96 | **1.41x** | 722.60 | 4714/5000 |
| `atom_pair(rdkit-compatible)` | 5000 | 18.15 | 6.84x | 14.50 | **8.44x** | 122.48 | 4249/5000 |
| `torsion(rdkit-compatible)` | 5000 | 4.15 | 24.79x | 3.69 | **27.03x** | 99.76 | 4248/5000 |
| `rdkit_fp(rdkit-compatible)` | 2000 | 1183.70 | 1.21x | 1175.94 | **1.18x** | 1392.57 | 2000/2000 (full) |
| `pattern_fp(rdkit-compatible)` | 2000 | 4174.27 | 0.06x | 225.23 | **1.32x** | 297.00 | 2000/2000 (full) |
| `tanimoto_1xN(per target, compatible Morgan)` | 50 | 0.57 | 0.26x | 0.06 | **1.33x** | 0.08 | 50/50 (full) |
| `has_substructure c1ccccc1` | 2000 | 242.99 | 0.02x | 1.67 | **2.33x** | 3.89 | 2000/2000 (full) |
| `has_substructure [OH]` | 2000 | 239.94 | 0.01x | 0.83 | **4.06x** | 3.39 | 2000/2000 (full) |
| `has_substructure C(=O)N` | 2000 | 264.26 | 0.02x | 1.51 | **2.74x** | 4.12 | 1997/2000 |
| `has_substructure [#7;R]` | 2000 | 251.38 | 0.01x | 2.99 | **1.10x** | 3.28 | 2000/2000 (full) |
| `has_substructure c1ccc2ccccc2c1` | 2000 | 344.49 | 0.03x | 8.14 | **1.17x** | 9.53 | 2000/2000 (full) |
| `has_substructure [CX3](=O)[OX2H1]` | 2000 | 267.55 | 0.02x | 1.62 | **2.27x** | 3.67 | 2000/2000 (full) |
| `has_substructure [NX3;H2,H1;!$(NC=O)]` | 2000 | 249.72 | 0.02x | 1.57 | **3.35x** | 5.26 | 2000/2000 (full) |
| `has_substructure *~*~*~*~*~*` | 2000 | 231.51 | 0.01x | 1.26 | **2.29x** | 2.89 | 2000/2000 (full) |
| `find_matches [#6]~[#7]` | 2000 | 261.00 | 0.08x | 2.54 | **7.30x** | 18.58 | 2000/2000 (full) |
| `murcko_scaffold` | 5000 | 221.45 | 0.39x | 22.12 | **3.53x** | 78.11 | 1591/5000 |
| `largest_fragment` | 2000 | 233.62 | 0.34x | 3.63 | **23.29x** | 84.56 | 2000/2000 (full) |
| `neutralize` | 2000 | 4.45 | 12.12x | 4.20 | **11.59x** | 48.72 | 2000/2000 (full) |
| `standardize vs Cleanup+Uncharge+Canonicalize` | 300 | 32037.67 | 0.88x | 1278.96 | **21.16x** | 27056.36 | 241/300 |
| `canonical_tautomer` | 300 | 28479.34 | 1.02x | 1592.06 | **17.05x** | 27141.84 | 241/300 |

## Same-task operations (not counted)

| Operation | Rows | main ratio | branch ratio | Why not counted |
|---|---:|---:|---:|---|
| `parse_smiles` | 5000 | 25.39x | 33.01x | RDKit sanitizes (rings/aromaticity) at parse; chematic defers perception |
| `canonical_smiles(prepared)` | 5000 | 0.74x | 1.45x | different canonical spelling |
| `parse+canonical_smiles` | 5000 | 2.34x | 3.81x | — |
| `mol_block_read` | 2000 | 0.47x | 3.42x | each library reads its own writer's block |
| `add_hydrogens` | 5000 | 6.48x | 4.72x | — |
| `remove_hydrogens` | 5000 | 17.04x | 13.80x | — |

## Not-equivalent operations (never counted)

| Operation | Rows | main ratio | branch ratio | Why not counted |
|---|---:|---:|---:|---|
| `mol_block_write` | 2000 | 206.07x | 209.40x | RDKit computes 2D coordinates; chematic writes without layout |
| `inchi` | 1000 | 0.45x | 1.85x | chematic default Mol.inchi is its own non-standard layering, not IUPAC Standard InChI (standard_inchi needs the native-inchi build); agreement is recorded, never counted (agreement 0/1000) |
| `inchikey` | 1000 | 0.35x | 2.16x | see inchi (agreement 0/1000) |
| `morgan_r2_2048(native ecfp4)` | 5000 | 0.19x | 1.67x | chematic native hash differs from RDKit Morgan |
| `standardize vs Cleanup` | 300 | 0.02x | 0.37x | chematic default also neutralizes and canonicalizes the tautomer |
| `brics_fragments` | 500 | 8.39x | 70.61x | chematic returns fragment Mols without RDKit's dummy-labelled SMILES set |
| `cip_labels` | 1000 | 1.63x | 46.98x | chematic default LegacyFast vs RDKit new CIP labeler |
| `sssr_rings` | 5000 | 0.13x | 1.35x | SSSR vs symmetrized SSSR |
| `svg_depiction` | 300 | 3.76x | 15.94x | different layout and renderer |
| `2d_layout` | 500 | 1.68x | 11.97x | different layout algorithms |
| `embed_3d` | 100 | 22.48x | 189.27x | no geometry-quality check; see the dedicated 3D records |
| `embed+minimize_mmff94` | 10 | 0.10x | 0.10x | legacy numeric-gradient chematic path; no quality check |
| `mcs_pair` | 50 | 2.08x | 2.24x | different defaults/timeouts; inputs not rebuilt |

`Mol.inchi` in the default build is chematic's own layering and does not
match IUPAC Standard InChI on any compared row; an InChI speed win must use
`standard_inchi` from a `native-inchi` build and is not claimed here.
`minimize_mmff94` is the legacy finite-difference path and remains ~10x
slower than RDKit's MMFF optimization; the analytic/stereo-safe MMFF94 lanes
have their own records.

## Output preservation (base vs branch)

[`scripts/perf_digest_diff.sh`](../scripts/perf_digest_diff.sh) builds
[`tools/perf_digest`](../tools/perf_digest/src/main.rs) against both
revisions and compares full per-row outputs byte for byte.

| Check | Rows (molecules × ops) | Differences | Record |
|---|---:|---:|---|
| 23 ops (SSSR, ring counts, LogP, MR, TPSA, HBA/HBD, rotatable bonds, ring bundle, ECFP4, RDKit-compatible Morgan, canonical SMILES, stereocenters, InChI, SMARTS match maps with and without uniquify, has-match, largest fragment, CIP, PAINS, Brenk) on 7 corpora, 46,736 molecules | 1,074,928 | 0 | [full](../validation/results/perf-digest-diff-01a86b62-8cc01365-full.json) |
| QED, MACCS, standardize on ChEMBL 5k + NCI 5k (9,999 molecules) | 29,997 | 0 | [heavy ops](../validation/results/perf-digest-diff-01a86b62-8cc01365-heavy.json) |

The two full 23-op digests have the same SHA-256
(`cb2a3f32f1a656419d4af5185ea81b4c1388361d96d4cb963c977c012e76bec6`).
The Rust test suites of the touched crates and the 985 Python binding tests
pass on the branch; `cargo check --workspace --all-targets` passes.

## Reproduction

```bash
# timing (install the wheel built from the revision under test first)
python3 scripts/bench_python_op_matrix_vs_rdkit.py \
  --corpus scripts/chembl_accuracy_corpus_4999.smi --limit 5000 --repeats 3 \
  --chematic-revision "$(git rev-parse HEAD)" --chematic-artifact path/to/chematic.whl \
  --output-json out.json

# differential outputs against the base revision
ONLY=qed,maccs,standardize scripts/perf_digest_diff.sh 01a86b62 heavy.json \
  scripts/chembl_accuracy_corpus_4999.smi scripts/nci_first_5k_smiles_only.smi
```

## Limits

- One 2-vCPU shared VM; no paired-bootstrap confidence intervals and no
  Apple-silicon rerun. Treat ratios near 1 as ties.
- One corpus for timing; substructure, fingerprint-pattern and transform
  rows use the first 300–2,000 rows as listed.
- WASM/RDKit.js and published packages were not measured.
