# Descriptor selected-column scaling — 2026-09-05

Status: local reproducibility evidence for the current `1.0.6` source tree.
No version, tag, or published-artifact claim is made here.

## Contract checked

`scripts/bench_descriptor_scaling.py` exercises the public
`chematic.bulk.descriptors_array` API with 3, 8, and all 60 supported columns.
For every lane it verifies:

- the returned column set matches the request;
- every column has one value per valid input molecule;
- repeated output has a stable SHA-256 digest; and
- Python-visible peak allocation is recorded with `tracemalloc`.

The mapping insertion order is intentionally not treated as a contract. Native
allocations inside the Rust extension are not visible to `tracemalloc` and are
reported as such.

## Measured lane

Environment: existing chematic `1.0.6` wheel in
`/private/tmp/chematic-speed-venv-2`, Apple arm64 host. Input is the first 250
SMILES from `scripts/chembl_accuracy_corpus_4999.smi`, one measured call per
column set.

| Requested columns | Rows | Seconds | Calls/s | Python peak bytes |
| ---: | ---: | ---: | ---: | ---: |
| 3 | 250 | 0.184132 | 5.4309 | 2,864 |
| 8 | 250 | 0.674350 | 1.4829 | 4,248 |
| 60 (all) | 250 | 3.472336 | 0.2880 | 17,307 |

Stable output digests for this exact fixture are:

| Requested columns | SHA-256 |
| ---: | --- |
| 3 | `2eb3fc6e54f7c94d58e22ad9ebc12615746cd33d60ce81fc6df2b673b1b2642` |
| 8 | `0b44753ec03cca9f78907e1e1b96ea09129cbefa812c83f0920ddb0194614768` |
| 60 | `6e389f612ba538f40a07a1a75061a575d9c2ba275187003538872fb47eab0902` |

Reproduce:

```text
/private/tmp/chematic-speed-venv-2/bin/python \
  scripts/bench_descriptor_scaling.py \
  scripts/chembl_accuracy_corpus_4999.smi \
  --limit 250 --repeats 1
```

The full 4,999-molecule run is a separate long-duration lane; this bounded
release-mode slice is the checked-in scaling and contract gate.
