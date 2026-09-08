# 2026-09-08 same-input MOL2 streaming comparison

The common cross-engine runner now includes Open Babel for the checked-in
single-record MOL2 fixture. This is record accounting, not a same-condition
speed claim: chematic uses its Rust reader, RDKit uses a Python block parser,
and Open Babel starts a CLI process for each repetition.

- Fixture: `benchmarks/fixtures/ethanol.mol2`
- SHA-256: `80144aad97332e2f58eb42bb2da1472587a9dfc2466f9d2fe081eb9834c2224d`
- Repetitions: 20; expected records: 20; source bytes: 361 per repetition
- Open Babel: `Open Babel 3.2.1 -- Jul 11 2026`
- RDKit: `2025.09.3`

| Engine | Records | Failures | Seconds | Records/s |
| --- | ---: | ---: | ---: | ---: |
| chematic | 20 | 0 | 0.001109 | 18,040.37 |
| RDKit | 20 | 0 | 0.480521 | 41.62 |
| Open Babel CLI | 20 | 0 | 3.578724 | 5.59 |

Reproduce with:

```text
python3 scripts/check_streaming_cross_engine.py --format mol2 --repeats 20 --openbabel obabel
```
