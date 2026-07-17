# chematic Benchmarks

Periodic performance snapshots. Each file is a date-stamped record of throughput and accuracy metrics at that version.

## Entries

| Date | Version | Notes |
|------|---------|-------|
| [2026-07-17](2026-07-17.md) | v0.4.29 | Hardware moved to Apple M4; throughput headline (5–14×) did not reproduce even on the same fixture — see file for details; descriptor accuracy holds |
| [2026-06-25](2026-06-25.md) | v0.4.20 | Baseline: ECFP4, descriptor batch, WASM size, RDKit accuracy |

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

Hardware varies by snapshot — see each dated file's header. 2026-06 was measured on Apple M2
(8-core, 8 GB RAM) / macOS 14; 2026-07 moved to Apple M4 (10-core, 16 GB RAM) / macOS 26.
Rust `cargo bench` numbers are from the same machine as the Python numbers in a given snapshot.

Results vary by CPU, load, and corpus/fixture choice — treat numbers as order-of-magnitude
references, not SLA guarantees. The 2026-07 snapshot found the headline throughput ratio is
fixture-sensitive; read the Notes section of each file before quoting a number.
