# Checked chi-square bit ranking

Base `832450783b3faff6ad455a4d514cadf6837b597d` after integrating main. Initial base
`0e7a6d9f8d50fd75a9b1bfada00336f8e320f8ec`. Python, typing, stubs and docs
checks were repeated after integration.
[RDKit #9601](https://github.com/rdkit/rdkit/issues/9601) is OPEN and tracks typed
explicit/sparse inputs to InfoTheory rankers. chematic already accepted those
inputs for ENTROPY, but refused CHISQUARE. This adds the missing uncorrected
Pearson statistic using the existing owned uint64 counts and validated votes.
It does not complete the upstream tracking issue or change RDKit.

## Results

- **225 Python tests passed**: 19 new chi-square cases, 44 existing entropy/input
  cases and 162 RDKit compatibility cases. The prior test that refused CHISQUARE
  is removed. Covers exact tables, empty rows/classes, one observed class,
  complements, reversed/doubled data, input/output ownership, errors, resource
  boundaries and exact counts beyond 65535. Includes the reduced RDKit duplicate-ID
  table as a regression; other unsupported metrics still refuse.
- **768 live RDKit 2026.09.1 rankers attempted**, explicit/sparse binary vectors,
  ENTROPY/CHISQUARE, 192 seeded datasets, 32,896 paired votes. **764** yield a
  complete upstream bit-ID permutation and match class counts exactly. Maximum
  absolute score difference is 1.9539925233402755e-14 (chi-square) and
  1.3322676295501878e-15 (entropy).
- **Four upstream CHISQUARE duplicate-ID outputs**, two datasets each exercised
  with explicit and sparse vectors, are recorded with full inputs/outputs in
  `oracle.json`. They are not counted as parity successes. In the reduced case,
  RDKit outputs bit 0 twice and omits bit 2; source inspection shows its
  cancellation-prone score expression can become negative near zero, after which
  the ranking excludes it and fills missing rows. This explains the observation
  by inference from source, not an instrumented upstream build. No such numerical
  bug is asserted to be the subject of #9601.
- Independent exact-rational Pearson checks pass for **2,064 distinct bit
  tables**, across both binary-vector input representations, including every
  reference anomaly. Maximum difference 1.7763568394002505e-15. **384** vote-order
  reversal, **384** duplication and **18** separate degenerate-ranker checks pass.
- Exact rank order differs in 244 of the 764 comparable rankers. The implementation
  uses descending computed score and ascending bit ID for exact ties. This is not
  exact-order or all-input RDKit parity.
- Strict mypy checks the implementation and a public consumer, with four expected
  invalid-consumer errors. Imported legacy package stubs are followed silently;
  this is not package-wide strict typing. Python 3.13 matches the installed NumPy
  2.5.3 stubs' syntax requirements.
- Source stub, cross-binding, binding, unsafe, refactor and product-name gates;
  strict MkDocs and two executed guide blocks. No Rust source changes or build.

## Contract and boundaries

Score = sum((observed - expected)^2 / expected) for the 2-by-class on/off table.
Expected-zero cells contribute zero. Constant bits/unobserved classes are valid;
zero-observation rankers remain errors. There is no Yates correction, p-value,
biasing, mask, chemistry/label inference, new sparse storage, Rust/WASM API or
performance claim. Default ENTROPY calculation and all existing limits remain:
1..65536 bits, 1..256 classes, <=1048576 count entries and <=1048576 records.
Output is independent float64 data, counts retained in uint64. Invalid votes
commit no partial changes.

A verified existing local native extension from `ca7e14e8` is reused only as the
package import fixture. It is hard-linked into a separate source package and used without modification;
its SHA-256 is compared with that task's recorded provenance. This is not a new
native build, published wheel or older-Python runtime validation. Candidate/source
hashes and exact upstream source revision are in `provenance.json`.

## Reproduction

Set `CANDIDATE_PY` to this source package with the existing extension, and
`RDKIT_PY` to RDKit 2026.09.1. Python 3.13.6 / pytest 9.1.1 / NumPy 2.5.3.

```sh
PYTHONPATH="$CANDIDATE_PY" python3.13 -m pytest \
  crates/chematic-py/tests/test_chi_square_bit_ranker.py \
  crates/chematic-py/tests/test_info_bit_ranker.py \
  crates/chematic-py/tests/test_rdkit_compat.py -q -p no:cacheprovider
PYTHONPATH="$CANDIDATE_PY:$RDKIT_PY" python3.13 \
  validation/results/chi-square-bit-ranker-2026-10-18/oracle.py
MYPYPATH="$CANDIDATE_PY" python3.13 -m mypy --no-incremental \
  --config-file validation/results/chi-square-bit-ranker-2026-10-18/mypy.ini \
  crates/chematic-py/python/chematic/_info_theory.py \
  validation/results/chi-square-bit-ranker-2026-10-18/typing_valid.py
# typing_invalid.py: expect four errors and exit 1.
python3.13 -m mkdocs build --strict --site-dir /tmp/chematic-chi-square-site
```

Independent formula/reference sources: RDKit
[`InfoBitRanker.cpp`](https://github.com/rdkit/rdkit/blob/7f3d37d16b71ff8c8a6b0de82c12a8b371e0a8df/Code/ML/InfoTheory/InfoBitRanker.cpp),
[`InfoGainFuncs.h`](https://github.com/rdkit/rdkit/blob/7f3d37d16b71ff8c8a6b0de82c12a8b371e0a8df/Code/ML/InfoTheory/InfoGainFuncs.h).
