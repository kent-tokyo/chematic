# Descriptor contract and streaming benchmark — 2026-09-05

Status: local reproducibility evidence for the current `1.0.6` source tree.
No version, tag, or published-artifact claim is made here.

## Shared descriptor contract

The checked-in `validation/cross_binding_contract.json` now carries descriptor
field semantics, units, source provenance, tolerances, and four fixtures.
Rust, Python, and Node/WASM consume the same descriptor fixture.

The 4,999-molecule core parity lane used the same SMILES input for chematic's
Python binding and RDKit's Python API:

| Field | Exact agreement | Agreement at published tolerance |
| --- | ---: | ---: |
| Molecular weight | 3,827/4,999 | 4,990/4,999 (±0.01 Da) |
| TPSA | 4,999/4,999 | 4,999/4,999 (±0.01 Å²) |
| HBD | 4,999/4,999 | 4,999/4,999 |
| HBA | 4,999/4,999 | 4,999/4,999 |
| Heavy atoms | 4,999/4,999 | 4,999/4,999 |

The molecular-weight residuals are small atomic-mass-table differences, not
parse failures. Full descriptor census remains a separate longer lane and is
not represented as complete by this result.

Reproduction:

```text
TMPDIR=/private/tmp python3 scripts/descriptor_core_parity.py \
  scripts/chembl_accuracy_corpus_4999.smi \
  --json validation/results/descriptor-core-parity-4999.json
```

## Streaming throughput

Release-mode chematic file-backed `BufRead` streaming, 2,000 passes over the
checked-in two-record/two-frame fixtures (4,000 records, zero failures):

| Format | Records/s | Bytes/s |
| --- | ---: | ---: |
| SDF | 60,839.93 | 19,255,838.79 |
| MOL/SDF blocks | 52,859.68 | 16,730,089.03 |
| XYZ | 51,137.19 | 4,116,543.80 |

The matching RDKit reference produced 25,464 SDF, 25,555 MOL, and 196,621
XYZ records/s, but it uses Python block constructors rather than the same
file-backed `BufRead` iterator. These numbers are therefore boundary evidence,
not a fair claim that one engine is faster than the other.

Commands:

```text
cargo run -p chematic-mol --example streaming_benchmark --release --offline -- --format sdf --path benchmarks/fixtures/streaming.sdf --repeats 2000
cargo run -p chematic-mol --example streaming_benchmark --release --offline -- --format mol --path benchmarks/fixtures/streaming.sdf --repeats 2000
cargo run -p chematic-mol --example streaming_benchmark --release --offline -- --format xyz --path benchmarks/fixtures/streaming.xyz --repeats 2000
TMPDIR=/private/tmp python3 scripts/bench_streaming_formats.py --repeats 2000
```
