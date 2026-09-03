# Canonical SMILES fast-path benchmark — 2026-09-04

This follow-up measures canonical SMILES generation only. Both engines receive
already parsed molecules, execute a warm-up pass, and canonicalize the same
5,000-entry corpus. It does not compare parsing, chemical-model agreement, or
canonical string equality between engines.

Environment: chematic 1.0.1 local candidate based on `f0058d2b`, RDKit
2025.09.3, Python 3.13.6, macOS 26.5.2 arm64. Corpus:
`scripts/descriptor_census_corpus.smi`, SHA-256
`d6f2ba3f128296f935007f0b0813aa97b6ebcc2457e014ddca2213ddd655276c`.

| Run | chematic (µs/mol) | RDKit (µs/mol) |
|---:|---:|---:|
| 1 | 26.01 | 25.76 |
| 2 | 25.23 | 25.58 |
| 3 | 25.38 | 25.78 |
| 4 | 24.72 | 25.37 |
| 5 | 24.73 | 25.57 |
| 6 | 24.95 | 25.86 |
| 7 | 24.76 | 25.20 |
| **Median** | **24.95** | **25.58** |

On this machine and corpus, chematic was 2.5% faster by the seven-run median
and faster in six of seven runs. Both engines produced a canonical result for
all 5,000 inputs. The margin is small and environment-specific; this is not a
claim that chematic is faster for every molecule class or platform.

An independent 5,000-entry ChEMBL-derived corpus
(`scripts/chembl_accuracy_corpus_4999.smi`, despite the historical filename)
was also measured five times. Its SHA-256 is
`1c47371dcbe37f4e0a141bf545b72bf238de2761fa3894fa251a552d84728d3e`.

| Run | chematic (µs/mol) | RDKit (µs/mol) |
|---:|---:|---:|
| 1 | 17.79 | 26.49 |
| 2 | 18.30 | 27.45 |
| 3 | 18.25 | 26.82 |
| 4 | 18.27 | 26.29 |
| 5 | 18.54 | 27.08 |
| **Median** | **18.27** | **26.82** |

Here chematic was 1.47× faster by the median and led all five runs. This
second result reduces corpus-specific uncertainty but remains a local
single-platform measurement.

The optimization avoids whole-molecule SSSR work for acyclic/exocyclic E/Z
double bonds, reuses Morgan-refinement storage, reduces canonical-partition
copying, proves simple local-twin orbits without invoking the general
automorphism search, and keeps common degree-four traversal buffers inline.
The canonical idempotency corpus and focused canonical regression suite passed
before measurement.

Reproduce after building and installing the local wheel:

```bash
maturin build --manifest-path crates/chematic-py/Cargo.toml --release --locked \
  --out /private/tmp/chematic-wheels
python3 -m pip install --force-reinstall --no-deps \
  /private/tmp/chematic-wheels/chematic-1.0.1-*.whl
for i in 1 2 3 4 5 6 7; do
  python3 scripts/bench_canonical.py scripts/descriptor_census_corpus.smi \
    --rdkit --n 5000 --json
done
for i in 1 2 3 4 5; do
  python3 scripts/bench_canonical.py scripts/chembl_accuracy_corpus_4999.smi \
    --rdkit --n 5000 --json
done
```
