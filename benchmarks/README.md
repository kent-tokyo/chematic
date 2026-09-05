# chematic Benchmarks

Periodic performance snapshots. Each file is a date-stamped record of throughput and accuracy metrics at that version.

## Entries

| Date | Version | Notes |
|------|---------|-------|
| [2026-09-05 hot-path 1.10x gate](2026-09-05-hotpath-110.md) | v1.0.6 local source | Seven alternating A/B pairs: canonical 1.176x, SDF read 1.180x, reused-buffer write 1.419x; parse 1.034x remains below target; exact-output checks and load caveats |
| [2026-09-05 descriptor/streaming](2026-09-05-descriptor-streaming.md) | v1.0.6 release source | Shared descriptor provenance and Rust/Python/Node/WASM fixture contract; 4,999-molecule core parity; 2,000-pass SDF/MOL/XYZ streaming evidence |
| [2026-09-05 prepared index](2026-09-05-prepared-index.md) | v1.0.6 local source | Exact reusable fingerprint index; 7.30x repeated-query speedup on the pinned ten-molecule fixture |
| [2026-09-05 parallel Tanimoto](2026-09-05-tanimoto-parallel.md) | v1.0.6 local source | Row-wise parallel dense matrix with serial parity; 1.21x on the pinned 256x256 lane |
| [2026-09-05 descriptor topology](2026-09-05-descriptor-topology.md) | v1.0.6 local source | Shared Wiener/Kappa/Chi topology context with scalar parity and lazy single-group fallback |
| [2026-09-05 distance descriptors](2026-09-05-distance-descriptors.md) | v1.0.6 local source | Shared AutoCorr2D/Moran/Geary distance matrix with exact scalar parity |
| [2026-09-05 descriptor scaling](2026-09-05-descriptor-scaling.md) | v1.0.6 local source | `descriptors_array` 3/8/all column contract, deterministic digest, and Python-visible allocation record |
| [2026-09-05 MMFF94 prepared nonbonded](2026-09-05-mmff94-prepared-nonbonded.md) | v1.0.6 local source | Prepared vdW parameters and electrostatic charge products with energy parity; analytic gradients remain open |
| [2026-09-05 MMFF94 parallel gradient](2026-09-05-mmff94-gradient-parallel.md) | v1.0.6 local source | Large-molecule finite-difference probes use bounded parallelism with deterministic gradient ordering |
| [2026-09-05 RDKit/Open Babel speed](2026-09-05-rdkit-openbabel-speed.md) | v1.0.6 local wheel | Three-sample in-process chematic/RDKit remeasurement and separate Open Babel CLI boundary lane |
| [2026-09-05 hot-path follow-up](2026-09-05-hot-path-follow-up.md) | v1.0.4-based unreleased source | Additional 1.10x gate: canonical SMILES 1.115x, file-backed SDF read 1.130x, V2000 SDF serialization 3.042x; output and full-workspace gates included |
| [2026-09-04 MMFF94/3D](2026-09-04-mmff94-3d.md) | v1.0.3 | Prepared MMFF94 energy, L-BFGS minimization, ETKDG generation, and 3D minimization microbenchmarks on macOS arm64 |
| [2026-09-04 streaming formats](2026-09-04-streaming-formats.md) | v1.0.3 | File-backed chematic SDF/MOL/XYZ streaming runner and explicitly non-equivalent RDKit block-parser reference |
| [2026-09-04 RDKit/Open Babel](2026-09-04-rdkit-openbabel.md) | v1.0.3 wheel + v1.0.4 source follow-up | Speed comparison with Open Babel 3.2.1 CLI plus source-level canonical and SDF writer A/B checks; CLI conversion remains separate from in-process measurements |
| [2026-09-04 WASM](2026-09-04-wasm-size.md) | v1.0.2 release candidate | Optimized web artifact size, gzip size, SHA-256 digest, tool versions, and exact reproduction commands |
| [2026-09-04 canonical](2026-09-04-canonical-fast-path.md) | v1.0.2 code, measured before metadata bump | Canonical SMILES comparisons on two 5,000-molecule corpora; chematic leads RDKit by 2.5% and 1.47× at the respective medians on the recorded macOS arm64 environment |
| [2026-09-04 SDF](2026-09-04-sdf-fast-path.md) | v1.0.2 code, measured before metadata bump | Graph/property read and serialization-only write improved 1.26× and 1.33× over the preceding chematic implementation; scoped RDKit comparison included |
| [2026-09-03](2026-09-03-competitive.md) | v1.0.1 | Resumable six-operation competitive run for chematic and RDKit; Open Babel not installed |
| [2026-08-23](2026-08-23.md) | v0.18.0 (commit `24a9239`) | v0.19.0 release-prep re-measurement; corpus now committed (`scripts/chembl_accuracy_corpus_4999.smi`); real MW check added; diverse-corpus ECFP4 now reproducible via `benchmark_vs_rdkit.py --corpus`; WASM size rebuilt clean; CIP R/S/E/Z label agreement re-measured (96%→99.7%+) |
| [2026-07-17](2026-07-17.md) | v0.4.29 | Hardware moved to Apple M4; throughput headline (5–14×) did not reproduce even on the same fixture — see file for details; descriptor accuracy holds |
| [2026-06-25](2026-06-25.md) | v0.4.20 | Baseline: ECFP4, descriptor batch, WASM size, RDKit accuracy |

## How to reproduce

### Throughput (Python)

```bash
# Corpus is committed — scripts/chembl_accuracy_corpus_4999.smi
pip install chematic rdkit
python scripts/bench5k.py scripts/chembl_accuracy_corpus_4999.smi
python scripts/bench5k.py scripts/chembl_accuracy_corpus_4999.smi --detail   # show mismatches
```

### Accuracy vs RDKit (4,999-mol ChEMBL-derived corpus)

`scripts/rdkit_benchmark.py` measures RDKit-side timing only (see its own docstring) — it
does not perform an accuracy comparison, despite an earlier version of this file pointing to
it for that purpose. The actual, current accuracy-reproduction path is:

```bash
pip install chematic rdkit
python scripts/bench5k.py scripts/chembl_accuracy_corpus_4999.smi --json /tmp/bench5k.json
python scripts/gen_validation_report.py /tmp/bench5k.json
```

See [`docs/validation.md`](../docs/validation.md) for the canonical, regeneratable report
this produces.

### WASM bundle size

```bash
# Requires wasm-pack
cd crates/chematic-wasm
wasm-pack build --target web --release
ls -lh pkg/chematic_wasm_bg.wasm
gzip -k pkg/chematic_wasm_bg.wasm && ls -lh pkg/chematic_wasm_bg.wasm.gz
```

## Hardware reference

Hardware varies by snapshot — see each dated file's header. 2026-06 was measured on Apple M2
(8-core, 8 GB RAM) / macOS 14; 2026-07 moved to Apple M4 (10-core, 16 GB RAM) / macOS 26.
Rust `cargo bench` numbers are from the same machine as the Python numbers in a given snapshot.

Results vary by CPU, load, and corpus/fixture choice — treat numbers as order-of-magnitude
references, not SLA guarantees. The 2026-07 snapshot found the headline throughput ratio is
fixture-sensitive; read the Notes section of each file before quoting a number.
