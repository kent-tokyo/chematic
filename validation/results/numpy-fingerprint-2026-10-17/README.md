# NumPy fingerprint output correction — 2026-10-17

Base: `f20304f4699ee803ab993f73b98c3c8be38b58fc` (main).
Related to RDKit #9601's remaining destination-array typing work. This fixes two
measured local adapter defects; installed RDKit already handles ordinary zero
bits correctly. It is not an upstream numeric bug fix or closure of that tracker.

## Before and after

The baseline adapter source is byte-identical to the pinned main file; its hash
is in `source-hashes.json`. With bit 2 set and destination `[7]*8`, baseline
chematic returns `[7,7,1,7,7,7,7,7]`; RDKit 2026.09.1 returns `[0,0,1,0,0,0,0,0]`.
An iterator that yields bit 1 then raises also leaves the baseline destination
partially modified. Both new regressions fail before the change, with exit code
1 (`baseline-tests.log`), and pass afterward.

The corrected converter validates a plain writable numeric NumPy array, extracts
all source IDs into independent bit storage, then overwrites every destination
element. The existing optional destination, array return and `dest` keyword stay
as adapter extensions. Exact size is required; there is no resizing. Source calls
and iterator errors propagate before conversion writes. Width and raw input IDs
are capped separately at 1,048,576. Source callbacks themselves have no timeout.

56 new tests plus 182 existing compatibility/folding tests pass (238 total).
Coverage includes stale zeros, exception atomicity, source/destination aliasing,
noncontiguous/reversed/foreign-endian arrays, bool/integer/float dtypes, read-only
and overlapping views, wrong shape/type/index, exact limits and bounded lazy views.

## Independent corpus

`validation/rdkit_issues/numpy-fingerprint/probe.py` constructs 576 datasets:
eight widths, three densities, four seeds and six destination dtypes. Five source
forms are tested: adapter ExplicitBitVect, RDKit ExplicitBitVect, RDKit
SparseBitVect, lazy IDs with duplicates, and NumPy integer ID views.

11,520 adapter calls (three destination layouts plus optional fresh output) match
an independent binary membership reference. 576 live RDKit ordinary writable
contiguous outputs agree. Invalid, read-only or strided destinations are not sent
to the upstream oracle; adapter views are checked against expected values and
backing padding is verified unchanged. This is binary projection, not molecular
fingerprint-generation accuracy or a performance benchmark.

The installed RDKit conversion does not accept **SparseBitVect** (binary); that
separate diagnostic is retained. Do not confuse this with RDKit's supported sparse
**count** vectors. The adapter's structural source support is a different contract.
No full RDKit API, dtype, resizing, source-type or return parity is claimed.

`corpus-results.json.gz` retains original bits, oracle and adapter outputs.
`summary.json` records versions and exact source/extension hashes.

## Typing and validation boundary

Valid source/destination mypy examples pass. Two invalid calls (missing source
protocol, list destination) produce exactly two argument-type errors after the
change. The same invalid call file passes against the old untyped signature;
logs record both. This is a local API typing check, not whole-package strict typing.

Changed Python sources parse with Python 3.9 grammar and compile. Existing source
stub, binding inventory, cross-binding manifest, product-spelling and strict
MkDocs gates pass. No Rust code/export, dependency or native stub changes are
included; Rust tests/Clippy and a new wheel build were not run for this correction.

Baseline and candidate Python package fixtures use the same prebuilt extension
from `0ed8958cb9a9e6a7140074ceb592e89061eb47cd`, hash recorded explicitly. The changed
conversion is pure Python/NumPy and does not call chematic Rust APIs. Existing
adapter regression tests use that extension fixture; they are not evidence of a
fresh source-main native build or publication. CI rebuilds the exact candidate
and discovers the new Python test file automatically.

## Reproduction

Copy current package Python files (including `_numpy_fingerprint.py`) next to a
compatible `chematic/chematic.so` in a temporary package directory. Put that
parent and an isolated RDKit 2026.09.1 package on PYTHONPATH. Existing installed
chematic packages were not overwritten. This host uses Python 3.13.6 and NumPy
2.5.3. External pytest plugins were disabled because an old host pytest-asyncio
plugin fails collection with pytest 9.

```sh
PYTEST_DISABLE_PLUGIN_AUTOLOAD=1 python3.13 -m pytest \
  crates/chematic-py/tests/test_numpy_fingerprint.py \
  crates/chematic-py/tests/test_rdkit_compat.py \
  crates/chematic-py/tests/test_fold_density.py -q -p no:cacheprovider
python3.13 validation/rdkit_issues/numpy-fingerprint/probe.py <output-directory>
MYPYPATH=crates/chematic-py/python python3.13 -m mypy \
  --follow-imports=silent --ignore-missing-imports --no-incremental \
  validation/rdkit_issues/numpy-fingerprint/typecheck_valid.py
# Using typecheck_invalid.py instead must fail with two arg-type errors.
```

Logs have trailing whitespace normalized. Native packaging, hosted CI and main
integration are separate from these local adapter checks.
