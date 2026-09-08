# Clean install and cold start — v1.0.9

Status: local release-candidate evidence. This is not a published-wheel or
cross-platform availability claim.

Environment: macOS arm64, CPython 3.13.6, Rust 1.97.0, maturin, offline
workspace build. The wheel was built from the current v1.0.9 source tree and
installed into a fresh virtual environment with dependencies disabled.

## Result

| Check | Result |
|---|---|
| Wheel | `chematic-1.0.9-cp313-cp313-macosx_11_0_arm64.whl` |
| Isolated install | passed with `pip install --no-deps` |
| Import | passed; reported version `1.0.9` |
| First parse smoke | passed; `CCO` -> `C(C)O` |
| Cold import median | 13.0 ms, 10 subprocess samples |
| Import + first parse median | 12.9 ms, 10 subprocess samples |

## Reproduction

```text
maturin build --release --offline -m crates/chematic-py/Cargo.toml -o /tmp/chematic-wheel-next
python3 -m venv /tmp/chematic-clean-venv
/tmp/chematic-clean-venv/bin/pip install --no-deps /tmp/chematic-wheel-next/chematic-1.0.9-cp313-cp313-macosx_11_0_arm64.whl
/tmp/chematic-clean-venv/bin/python scripts/bench_startup.py --json --runs 10
```

The startup numbers include a fresh Python process per sample. They are
machine- and filesystem-dependent and do not establish throughput, peak
memory, or competitor parity.
