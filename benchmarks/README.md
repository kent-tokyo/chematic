# chematic Benchmarks

Periodic performance snapshots. Each file is a date-stamped record of throughput and accuracy metrics at that version.

## Entries

| Date | Version | Notes |
|------|---------|-------|
| [2026-06](2026-06.md) | v0.4.20 | Baseline: ECFP4, descriptor batch, WASM size, RDKit accuracy |

## How to reproduce

### Throughput (Python)

```bash
# Requires a SMILES CSV with a 'SMILES' column (e.g. ChEMBL export)
pip install chematic rdkit
python scripts/bench5k.py ~/Downloads/SMILES.csv
python scripts/bench5k.py ~/Downloads/SMILES.csv --detail   # show mismatches
```

### Accuracy vs RDKit (175-mol drug-like corpus, in-repo)

```bash
pip install chematic rdkit pandas
python scripts/rdkit_benchmark.py
```

### WASM bundle size

```bash
# Requires wasm-pack
cd crates/chematic-wasm
wasm-pack build --target web --release
ls -lh pkg/chematic_wasm_bg.wasm
gzip -k pkg/chematic_wasm_bg.wasm && ls -lh pkg/chematic_wasm_bg.wasm.gz
```

## Hardware reference

All Python benchmarks were measured on **Apple M2 (8-core, 8 GB RAM)** running macOS 14, Python 3.13, single process. Rust `cargo bench` numbers are from the same machine.

Results vary by CPU and load — treat numbers as order-of-magnitude references, not SLA guarantees.
