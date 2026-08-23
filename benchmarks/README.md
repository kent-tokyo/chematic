# chematic Benchmarks

Periodic performance snapshots. Each file is a date-stamped record of throughput and accuracy metrics at that version.

## Entries

| Date | Version | Notes |
|------|---------|-------|
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
