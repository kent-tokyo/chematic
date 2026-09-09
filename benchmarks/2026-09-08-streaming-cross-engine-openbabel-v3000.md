# 2026-09-08 same-input V3000 streaming comparison

The common cross-engine runner now includes Open Babel for the checked-in
single-record V3000 fixture. This is record accounting, not a same-condition
speed claim: chematic uses its Rust reader, RDKit uses a Python block parser,
and Open Babel starts a CLI process for each repetition.

- Fixture: `benchmarks/fixtures/ethanol.v3000`
- SHA-256: `64ac11b3de1c186deae555d54d6a47e7f3c72607ef0a90b0d93f467496bbd3ad`
- Repetitions: 20; expected records: 20; source bytes: 298 per repetition
- Open Babel: `Open Babel 3.2.1 -- Jul 11 2026`
- RDKit: `2025.09.3`

| Engine | Records | Failures | Seconds | Records/s |
| --- | ---: | ---: | ---: | ---: |
| chematic | 20 | 0 | 0.001083 | 18,472.90 |
| RDKit | 20 | 0 | 0.435837 | 45.89 |
| Open Babel CLI | 20 | 0 | 3.866521 | 5.17 |

Reproduce with:

```text
python3 scripts/check_streaming_cross_engine.py --format v3000 --repeats 20 --openbabel obabel
```
