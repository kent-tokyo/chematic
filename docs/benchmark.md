# Benchmark

## v1.0.1 benchmark preparation status

The repository now contains a machine-readable competitive benchmark protocol
at [`validation/competitive_benchmark_manifest.json`](../validation/competitive_benchmark_manifest.json).
It defines the three in-scope engines (chematic, RDKit, and Open Babel),
corpora, operations, required provenance, fairness
rules, and result statuses needed for a defensible comparison.

This is preparation evidence, not a new performance claim. The dated numbers
below remain historical snapshots until the protocol is executed again on the
current release. Validate the protocol offline with:

```bash
python3 scripts/validate_competitive_benchmark_manifest.py
```

An advantage claim is publishable only when the result includes the exact
source revision, engine versions, hardware, corpus hash, configuration,
failure counts, and reproduction command. Throughput, accuracy, memory,
startup latency, deployment size, and feature coverage must be reported as
separate dimensions.

The execution coordinator is resumable per operation. It writes an atomic
`state.json` and one log per operation; completed operations are skipped and
failed or interrupted operations can be retried:

```bash
python3 scripts/run_competitive_benchmark.py --dry-run
python3 scripts/run_competitive_benchmark.py --resume
```

An interrupted run must be resumed with the same result directory and state
file. The runner stops on the first failed operation while preserving its
status and log.

The runner also fails before timing if the imported Python package version does
not match the manifest target version. This prevents a stale installed wheel
from being reported as a current-workspace result.

## Latest measured snapshot

### 2026-09-03 — chematic 1.0.1

The current protocol run completed for chematic and RDKit on macOS arm64.
Open Babel was not installed, so no Open Babel comparison is claimed. The
results are operation-specific: chematic led import, SMILES parsing, and
ECFP4, while RDKit led canonical SMILES and SDF read/write. See the complete
record and hashes in [`benchmarks/2026-09-03-competitive.md`](../benchmarks/2026-09-03-competitive.md).

The follow-up SDF fast-path measurement separates graph/property parsing from
optional diagnostics and separates serialization from automatic 2D layout.
On the same 365-record corpus, the fast read path is 6.8× faster and
serialization-only write is 8.0× faster than RDKit. See [`benchmarks/2026-09-04-sdf-fast-path.md`](../benchmarks/2026-09-04-sdf-fast-path.md).

The older snapshot below is retained as historical evidence and must not be
combined with this run.

This page reports the latest reproducible benchmark run, not necessarily the
current package version. The benchmark snapshot below was measured before the
workspace moved to v0.31.0; use the version and commit shown here when
reproducing these numbers.

**chematic v0.18.0** (commit `24a9239`, pre-`v0.19.0` bump) · **RDKit 2026.03.4**
· Python 3.13.6 · Apple M4 (10-core, 16 GB RAM) · macOS 26.5.2 · measured 2026-08-23.

Full methodology, raw numbers, and what changed since the last snapshot:
[`benchmarks/2026-08-23.md`](../benchmarks/2026-08-23.md). Snapshot history:
[`benchmarks/README.md`](../benchmarks/README.md).

| Metric | chematic | RDKit |
|--------|----------|-------|
| Import time | **14.6 ms** | 98.1 ms (6.7×) |
| SMILES parse — 5,000 mol (repeated fixture) | **1.0 µs/mol** | 17.7 µs/mol (17.7×) |
| ECFP4 — 10,000 mol (repeated fixture) | **6.76 µs/mol** | 44.48 µs/mol (6.6×) |
| ECFP4 — 5,000 mol (diverse ChEMBL corpus) | **54.7 µs/mol** | 94.3 µs/mol (1.7×) |
| Descriptor accuracy vs RDKit | **20 metrics ≥98.6%, 16 at 100%** (4,999-mol corpus) | baseline |
| CIP R/S/E/Z label agreement | **99.74–99.78%** | baseline |
| Install | `pip install chematic` | `pip install rdkit` (official prebuilt wheels) or conda |
| C/C++ dependencies | **Zero**, even building from source | Not required for the prebuilt wheel; required (Boost, CMake) building from source |
| WASM binary size | **2.98 MB raw / 1.11 MB gzip** | 6.91 MB raw / 2.06 MB gzip (RDKit.js) |

The ECFP4 speedup is fixture-dependent: a small set of simple molecules repeated
to fill the batch shows a larger ratio (6.6×) than a large, structurally diverse
corpus (1.7×) — **never blend these two numbers**. See
[`benchmarks/2026-08-23.md`](../benchmarks/2026-08-23.md) for why, and for the
full record of what moved since the previous (v0.4.29-era) snapshot.

---

## 1. Startup Time (import latency)

Cold-process import time measured by spawning a fresh Python subprocess per
sample (10 samples, median reported). No module-cache warm-up.

| | chematic | RDKit |
|---|---|---|
| `import` only | **14.6 ms** | 98.1 ms |
| `import` + first parse | **14.6 ms** | 103.3 ms |
| **Speedup** | **6.7×** | baseline |

**Why chematic is faster**: chematic is a single PyO3 extension module with no
transitive Python dependencies. RDKit initialises multiple C++ modules, reads
SMARTS definition files, and triggers Boost data-structure setup on first import.

### How to reproduce

```bash
python scripts/bench_startup.py --rdkit --runs 10 --json
```

---

## 2. SMILES Parse Throughput

Timed on the built-in 20-molecule diverse corpus repeated to N total parses.
Warm-up pass excluded.

| N | chematic | RDKit | Speedup |
|---|---|---|---|
| 1,000 | 1.11 µs/mol | 16.38 µs/mol | 14.9× |
| 5,000 | **1.00 µs/mol** | **17.74 µs/mol** | **17.7×** |
| 10,000 | 0.98 µs/mol | 16.00 µs/mol | 16.3× |

**chematic**

```python
import chematic
mols = [chematic.from_smiles(s) for s in smiles_list]
# or batch:
mols = chematic.from_smiles_list(smiles_list)
```

**RDKit**

```python
from rdkit import Chem
mols = [Chem.MolFromSmiles(s) for s in smiles_list]
```

### How to reproduce

```bash
python scripts/bench_smiles_parse.py --n 5000 --rdkit --json
```

**Related micro-benchmarks** (same corpus-tiering convention, not folded into the table
above): `python scripts/bench_canonical.py --rdkit` measures canonical SMILES generation
throughput; `python scripts/bench_smarts.py --rdkit` measures SMARTS substructure-match
throughput (pairs with `scripts/rdkit_benchmark.py`'s `bench_smarts_match`). Not
re-measured this cycle.

---

## 3. Speed — ECFP4 Fingerprint Generation (batch)

Rayon parallelism across all CPU cores.

**Repeated small fixture** (20 small drug-like SMILES cycled to fill each N):

| Molecules (N) | chematic (`bulk.ecfp4`) | RDKit (Python loop, deprecated `GetMorganFingerprintAsBitVect`) | Speedup |
|---------------|------------------------|---------------------|---------|
| 100 | 9.16 µs/mol | 46.98 µs/mol | 5.1× |
| 1,000 | 6.75 µs/mol | 44.86 µs/mol | 6.6× |
| 10,000 | 6.76 µs/mol | 44.48 µs/mol | **6.6×** |

**Diverse corpus** (5,000 unique molecules,
`scripts/descriptor_census_corpus.smi`):

| | chematic | RDKit | Speedup |
|--|----------|-------|---------|
| ECFP4 | **54.7 µs/mol** | 94.3 µs/mol | **1.7×** |

**Never blend these two numbers.** Unique, structurally diverse molecules
consistently show a smaller margin than a small fixture repeated many times —
see [`benchmarks/2026-08-23.md`](../benchmarks/2026-08-23.md#notes) for the
mechanism this project has partially traced (SSSR ring-perception cost is a
non-trivial share of `ecfp4`'s per-molecule cost), not fully root-caused.

**chematic**

```python
import chematic
fps = chematic.bulk.ecfp4(smiles_list)  # (N, 2048) uint8 numpy array
```

**RDKit**

```python
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors
fps = [rdMolDescriptors.GetMorganFingerprintAsBitVect(
           Chem.MolFromSmiles(s), 2, 2048)
       for s in smiles_list]
```

### How to reproduce

```bash
# repeated small fixture:
python scripts/benchmark_vs_rdkit.py --rdkit --n 10000 --repeats 5
# diverse corpus:
python scripts/benchmark_vs_rdkit.py --rdkit --corpus scripts/descriptor_census_corpus.smi --repeats 5
```

---

## 4. Descriptor Accuracy vs RDKit

Tested on a 4,999-molecule ChEMBL-derived SMILES corpus
(`scripts/chembl_accuracy_corpus_4999.smi`, committed to this repo). See
[Validation](validation.md) for the full per-metric breakdown, the
stereocenter oracle-calibration detail, and the separate CIP R/S/E/Z
label-agreement measurement.

| Descriptor | Agreement | Tolerance |
|-----------|-----------|-----------|
| Molecular weight | 99.82% | ±0.01 Da |
| Heavy atom count | **100%** | exact |
| H-bond donors (HBD) | **100%** | exact |
| H-bond acceptors (HBA) | **100%** | exact |
| TPSA | **100%** | ±0.1 Å² |
| LogP (Crippen) | **100%** | exact* |
| MR (molar refractivity) | **100%** | ±0.01 |
| Fsp3 | **100%** | ±0.001 |
| Aromatic ring count | **100%** | exact |
| Aliphatic ring count | **100%** | exact |
| Saturated ring count | **100%** | exact |
| Rotatable bonds | **100%** | exact |
| Num heteroatoms | **100%** | exact |
| Num spiro atoms | **100%** | exact |
| Num bridgehead atoms | **100%** | exact |
| Num amide bonds | **100%** | exact |
| Aromatic/aliphatic heterocycles | **100%** | exact |
| Num stereocenters (legacy)  | **99.96%** | exact |
| Num stereocenters (new CIP) | 98.6% | exact |
| [nH] SMARTS match | **100%** | precision/recall |
| CIP R/S/E/Z label agreement | 99.74–99.78% | exact (separate metric — see below) |

20 of 20 descriptor/count metrics reach ≥98.6% on the 4,999-molecule ChEMBL
corpus (RDKit 2026.03.4). Molecular weight is a **new check this release** —
`bench5k.py` never actually measured it before (only monoisotopic exact mass);
see [Validation](validation.md) for that fix and the corrected stereocenter
oracle-disagreement accounting.

**CIP R/S/E/Z label agreement is a distinct question from stereocenter
*count* agreement above**: given a stereocenter both engines agree exists,
does chematic assign the same R/S/E/Z label? 99.74% vs a modern
`rdCIPLabeler` oracle, 99.78% vs the legacy `_CIPCode` oracle — both a large
improvement over a prior snapshot's 96.30–96.83%, reflecting CIP-engine work
landed in the interim. Full detail: [Validation](validation.md).

### How to reproduce

```bash
python scripts/bench5k.py scripts/chembl_accuracy_corpus_4999.smi --detail
```

---

## 5. Installation & Deployment

| | chematic | RDKit |
|---|----------|-------|
| Python | `pip install chematic` | `pip install rdkit` (official prebuilt wheels) or `conda install -c conda-forge rdkit` |
| C/C++ compiler | Not required, even building from source | Not required for the prebuilt wheel; required (Boost) building RDKit itself from source |
| Docker image size delta | ~4 MB (approximate; not independently re-measured this cycle) | ~200 MB+ (approximate; not independently re-measured this cycle) |
| GitHub Actions | Single pip line | Separate conda setup step |
| JavaScript / WASM | `npm install @kent-tokyo/chematic` (2.98 MB raw / 1.11 MB gzip) | `npm install @rdkit/rdkit` (RDKit.js, a separate community project — 6.91 MB raw / 2.06 MB gzip) |
| Browser deployment | Yes | Yes, via RDKit.js |

---

## 6. Feature Comparison

| Feature | chematic | RDKit |
|---------|----------|-------|
| pKa prediction | Built-in, rule-based screening — not for clinical use (23 SMARTS rules) | External tool required |
| ADMET profile (BBB, Caco-2, hERG, CYP3A4) | Built-in, rule-based screening — not for clinical use | External tool required |
| MCP server (AI agent integration) | 20 tools (stdio only) | Not available |
| LSH approximate nearest-neighbour index | Built-in | Not available |
| IUPAC name generation | Built-in (offline) | Not available |
| Browser / WASM deployment | Yes (2.98 MB raw / 1.11 MB gzip) | Yes, via RDKit.js (a separate community project — 6.91 MB raw / 2.06 MB gzip) |
| ECFP4 batch speed | 1.7× faster (diverse corpus), 6.6× faster (repeated fixture) | Baseline |
| SMARTS atom map `:N` | Yes | Yes |
| Retrosynthesis (template-based) | 60 retro-SMIRKS built-in | External tool |
| File formats | 20+ | 100+ |
| 3D conformer generation | **Experimental** — distance geometry + torsion-aware pipeline (`embed_pipeline_v2`); the dated 265-molecule comparison is summarized in the benchmark snapshot above, see [`rdkit-migration.md`](rdkit-migration.md) for the limitations | Mature — ETKDGv3 with ML-assisted torsion corrections |
| Community & publications | Growing | Established (20+ years) |

---

## 7. Batch Descriptor Computation

`chematic.bulk.descriptors` returns 55+ descriptors per molecule including ADMET and pKa — all in parallel.
`chematic.bulk.descriptors_array` returns selected columns as numpy arrays (~25% faster for column-oriented access).

| N | chematic (`bulk.descriptors`) | Descriptors per call |
|---|-------------------------------|---------------------|
| 100 | 63 ms (629 µs/mol) | 55+ (incl. pKa, ADMET) |
| 1,000 | 605 ms (605 µs/mol) | 55+ (incl. pKa, ADMET) |

`descriptors_array(smiles, columns)` runs this same full pipeline regardless of
the requested column subset — selecting fewer columns does not reduce the
computation, only the returned dict.

```python
import chematic
import pandas as pd

# list-of-dicts (general purpose)
df = pd.DataFrame(chematic.bulk.descriptors(smiles_list))

# columnar numpy arrays (faster for specific columns)
result = chematic.bulk.descriptors_array(smiles_list, ["mw", "logp", "tpsa", "hba"])
df = pd.DataFrame(result)   # float64 / bool arrays, no per-molecule dict overhead
```

### Compound screening

```python
# One call bundles lipinski / veber / pains / brenk / qed / sa_score:
results = chematic.screen(smiles_list, profile="druglike")
passing = [r for r in results if r["overall_pass"]]
```

### Large SDF files (streaming)

```python
# iter_sdf() streams one record at a time — no full-file load:
for rec in chematic.iter_sdf("large.sdf"):
    print(rec.smiles, rec.get("Activity"))

# batch pipeline:
for batch in chematic.iter_sdf_batched("large.sdf", batch_size=1000):
    descs = chematic.bulk.descriptors([r.smiles for r in batch])
```

RDKit's `rdkit.Chem.Descriptors.CalcMolDescriptors` covers ~200 descriptors but does not include pKa or ADMET.

### How to reproduce

```bash
python scripts/benchmark_vs_rdkit.py --rdkit --n 10000 --repeats 5
```

---

## Reproduction & Raw Results

Every number on this page traces to one of these commands, run against
commit `24a9239` (Python 3.13.6, RDKit 2026.03.4, Apple M4, macOS 26.5.2,
2026-08-23):

```bash
# Environment (fix once if the editable install is stale):
.venv/bin/maturin develop --release -m crates/chematic-py/Cargo.toml

# Startup / import
python scripts/bench_startup.py --rdkit --runs 10 --json

# SMILES parse throughput
python scripts/bench_smiles_parse.py --n 1000 --rdkit --json
python scripts/bench_smiles_parse.py --n 5000 --rdkit --json
python scripts/bench_smiles_parse.py --n 10000 --rdkit --json

# ECFP4 + bulk.descriptors (repeated small fixture)
python scripts/benchmark_vs_rdkit.py --rdkit --n 10000 --repeats 5

# ECFP4 (diverse corpus)
python scripts/benchmark_vs_rdkit.py --rdkit --corpus scripts/descriptor_census_corpus.smi --repeats 5

# Descriptor accuracy (4,999-mol ChEMBL corpus)
python scripts/bench5k.py scripts/chembl_accuracy_corpus_4999.smi --json validation/results/bench5k_latest.json
python3 scripts/gen_validation_report.py validation/results/bench5k_latest.json

# CIP R/S/E/Z label agreement
cargo run -p chematic-cip --release --example corpus_snapshot -- \
    --candidate scripts/chembl_accuracy_corpus_4999.smi /tmp/candidate.tsv
python scripts/cip_accurate_full_corpus_report.py \
    /tmp/candidate.tsv /tmp/candidate.tsv scripts/chembl_accuracy_corpus_4999.smi

# WASM bundle size
cd crates/chematic-wasm && rm -rf pkg
wasm-pack build --target web --release
wasm-opt -O3 -o pkg/chematic_wasm_bg.wasm pkg/chematic_wasm_bg.wasm
ls -lh pkg/chematic_wasm_bg.wasm
gzip -k pkg/chematic_wasm_bg.wasm && ls -lh pkg/chematic_wasm_bg.wasm.gz
```

Raw output and full methodology notes (including everything that changed
since the v0.4.29-era snapshot):
[`benchmarks/2026-08-23.md`](../benchmarks/2026-08-23.md). Snapshot history
and hardware notes: [`benchmarks/README.md`](../benchmarks/README.md).
Machine-readable JSON: `validation/results/bench5k_latest.json`.
