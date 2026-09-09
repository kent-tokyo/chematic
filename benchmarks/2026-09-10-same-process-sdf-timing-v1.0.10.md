# v1.0.10 same-process SDF timing context

This is a bounded equivalent-operation timing context, not a general
performance ranking. The current source-built schematic Python extension and
RDKit 2025.09.3 were loaded in one CPython 3.13 process and received the same
633-byte SDF fixture. Warm-up was performed before 40 alternating rounds, with
25 parses per timed sample.

| Lane | P50 ns/parse | P95 ns/parse |
| --- | ---: | ---: |
| schematic v1.0.10 | 14,176.7 | 16,825 |
| RDKit 2025.09.3 | 374,594.2 | 501,886 |

The observed RDKit/schematic P50 ratio is 26.4233. Both lanes parsed two
records with zero failures. This result is limited to the checked-in SDF
fixture and these Python/file-backed parser boundaries; it does not establish
cross-format, browser, Open Babel, or production-workload performance.

Reproduce with a clean temporary environment containing the locally built
wheel:

```sh
maturin build --release --offline -m crates/chematic-py/Cargo.toml \
  --out /private/tmp/chematic-wheel-v1.0.10
python3 -m venv --system-site-packages /private/tmp/chematic-timing-venv
/private/tmp/chematic-timing-venv/bin/python -m pip install --no-deps \
  /private/tmp/chematic-wheel-v1.0.10/*.whl
PYTHONPATH=scripts /private/tmp/chematic-timing-venv/bin/python \
  scripts/benchmark_same_process_sdf_timing.py \
  --output benchmarks/2026-09-10-same-process-sdf-timing-v1.0.10.json
```

Machine-readable evidence: [`2026-09-10-same-process-sdf-timing-v1.0.10.json`](2026-09-10-same-process-sdf-timing-v1.0.10.json).
