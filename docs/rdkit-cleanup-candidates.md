# Sparse cleanup candidates

[RDKit #9475](https://github.com/rdkit/rdkit/issues/9475), confirmed open on
October 19, 2026, tracks incremental SMILES parsing and sanitization speedups.
This change targets chematic's native RDKit-compatible `clean_up` stage:
neutral valence-five nitrogen corrections followed by phosphorus/halogen-oxide
corrections. It is independent of the earlier conjugation-stage work.

The former code allocated one candidate flag per molecule atom and scanned the
whole flag array a second time. It also cloned each eligible nitrogen's bond
list in both nitrogen passes. The candidate now keeps eligible atom indices in
a local inline list, visits those indices directly, and borrows unchanged
adjacency lists. Up to eight eligible nitrogens fit inline; larger lists spill
to the heap. Cleanup results and public APIs are unchanged.

## Preserve two-pass chemistry

First-pass discovery order and membership are retained exactly. In particular,
an N=O correction can change a candidate's charge before the later N#N pass;
filtering candidates again using the new charge would change the chemistry.
The implementation keeps the original candidate and preserves bond-neighbour
order, aromatic-flag restoration, property-cache updates and the final
phosphorus/halogen pass.

The previous function is frozen as a test-only reference at base `b7004777`.
The differential check compares the full debug representation of the internal
molecule and the returned result, twice in succession, rather than only its
canonical SMILES. This covers bond order, charges, stereo/source fields,
adjacency, rings and stored atom properties.

The existing 1,131-row canonical fixture supplies 1,130 raw projected internal
models and one explicitly recorded malformed-input refusal. All 1,130 states
agree after each cleanup pass. These are states **before full sanitization**;
accepting a raw model does not imply that its complete sanitization succeeds.

Original controls exercise both nitrogen phases, already-charged and aromatic
atoms, phosphorus/halogens, isotopes, explicit H, chirality, metal/dative bonds,
wildcards, changed aromatic/cache fields, empty input, candidate counts around
the inline capacity and sparse molecules up to 2,051 atoms. Control provenance,
license, minimization, oracle and hashes are retained in the validation record.

## One-host stage measurement

Both functions run in the **same release test binary** on identical fresh raw
internal models. Parsing, model restoration/cloning and destruction are
outside the timer. After warming both paths, nine paired rounds alternate the
measurement order; each round applies cleanup to every model 300 times.

| Packet | Models | Frozen median | Candidate median | Time reduction |
| --- | ---: | ---: | ---: | ---: |
| Existing canonical fixture | 1,130 | 11.138377 ms | 4.814780 ms | 56.77% |
| Original branch controls | 22 | 0.405923 ms | 0.175412 ms | 56.79% |

These are cleanup-stage totals on macOS arm64 with Rust 1.97.0, not whole-parser
or full-sanitization timings. They do not establish an end-to-end gain, RDKit
speed superiority, Enamine REAL representativeness or a released-package result.

Reproduce from the candidate checkout:

```sh
python3 scripts/measure_cleanup_candidates.py \
  --target-dir /private/tmp/chematic-cleanup-target \
  --output /private/tmp/chematic-cleanup-results
```

The script uses offline Cargo, release optimization with debug/LTO/incremental
disabled and 16 codegen units. It records compiler/platform, settings, source
and fixture hashes, every paired observation and refusal counts. It runs timing
alone after differential validation.

Validation: **572 Rust tests passed, 0 failed, 2 ignored** across 28 test
binaries. The three new functional tests and the opt-in timing test are part
of the existing library target; ordinary CI runs the functional tests.
See the [dated evidence](https://github.com/kent-tokyo/chematic/tree/main/validation/results/cleanup-candidates-2026-10-19).
