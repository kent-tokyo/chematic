# Clean install, cold start, throughput, and peak memory — v1.0.10

Status: local current-candidate evidence. This is not a published-wheel or
cross-platform availability claim.

Environment: macOS arm64, CPython 3.13, Rust 1.97.0, maturin 1.13.3,
offline workspace build. The wheel was installed into a fresh virtual
environment with `pip install --no-deps`.

| Check | Result |
|---|---|
| Wheel | `chematic-1.0.10-cp313-cp313-macosx_11_0_arm64.whl` (4,947,923 bytes) |
| Isolated install | passed |
| Import | passed; reported version `1.0.10` |
| First parse smoke | passed; `CCO` -> `C(C)O` |
| Cold import median | 35.6 ms, 10 subprocess samples |
| Import + first parse median | 35.7 ms, 10 subprocess samples |
| SMILES parse throughput | 154,220 molecules/s, 5,000 rows, 6.48 µs/mol |
| Peak RSS | 25,051,136 bytes for 5,000 parses |
| Python binding test suite | 864 passed, 0 failed, 6.27 s |

## Reproduction

```text
maturin build --release --offline -m crates/chematic-py/Cargo.toml -o /tmp/chematic-wheel-next
python3 -m venv /tmp/chematic-clean-venv
/tmp/chematic-clean-venv/bin/pip install --no-deps /tmp/chematic-wheel-next/chematic-1.0.10-cp313-cp313-macosx_11_0_arm64.whl
/tmp/chematic-clean-venv/bin/python scripts/bench_startup.py --json --runs 10
/tmp/chematic-clean-venv/bin/python scripts/bench_smiles_parse.py scripts/descriptor_census_corpus.smi --n 5000 --json
/tmp/chematic-clean-venv/bin/python -c 'import chematic,resource; [chematic.from_smiles("CC(=O)Oc1ccccc1C(=O)O") for _ in range(5000)]; print(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)'
PYTHONPATH=/tmp/chematic-clean-venv/lib/python3.13/site-packages:/Library/Frameworks/Python.framework/Versions/3.13/lib/python3.13/site-packages /tmp/chematic-clean-venv/bin/python -m pytest -p no:asyncio crates/chematic-py/tests --tb=short
```

The startup numbers use a fresh Python process per sample. Peak RSS is a
single-process CPython measurement on this host, not a browser-memory or
million-molecule extrapolation. Competitor measurements remain separate
unless they use the same environment and operation boundary.
