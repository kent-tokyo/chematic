# Checked NumPy fingerprint output

[RDKit #9601](https://github.com/rdkit/rdkit/issues/9601) lists destination-array
typing for `ConvertToNumpyArray` among remaining wrapper work. Investigating that
contract found a local chematic compatibility-layer bug: only on-bits were
written, leaving old values in every off-bit position. Installed RDKit 2026.09.1
already overwrites those positions correctly; this is not an upstream numeric bug.

## Use

```python
import numpy as np
from chematic.rdkit_compat import DataStructs, ExplicitBitVect

fp = ExplicitBitVect(8)
fp.SetBit(2)
dest = np.full(8, 7, dtype=np.int8)
assert DataStructs.ConvertToNumpyArray(fp, dest) is dest
assert dest.tolist() == [0, 0, 1, 0, 0, 0, 0, 0]

# Existing chematic extension: an omitted destination returns a new int8 array.
fresh = DataStructs.ConvertToNumpyArray(fp)
```

Every element, including zero bits, replaces the corresponding old value. The
source is completely read and validated before this converter writes to `dest`;
source-call/iterator failures and validation errors therefore do not cause partial
writes. This also works when a supplied on-bit view aliases the destination.

## Input and output contract

Sources provide callable `GetNumBits()` and `GetOnBits()` methods. Width and bit
IDs use the integer index protocol (including NumPy integer scalars), with Python
and NumPy booleans refused. Indices must be in `[0, width)`. Order and duplicate
IDs do not change binary membership; counts are not represented.

The source protocol supports this adapter's `ExplicitBitVect`, installed RDKit's
explicit and sparse binary vectors, and lazy or NumPy-based on-bit views.
Installed RDKit's own conversion overload does **not** accept `SparseBitVect` on
this host; accepting that structural source is an adapter extension. Count
vectors and real-valued fingerprints are outside this binary contract.

A supplied destination must be a **plain `numpy.ndarray`**, one-dimensional and
exactly the source width, writable, with boolean/signed-integer/unsigned-integer/
floating dtype. Disjoint strided and reversed views and non-native byte order
are supported. Elements whose byte ranges overlap are refused. Array subclasses,
object/string/complex/datetime arrays, lists and byte buffers are refused.
No resizing or implicit container conversion is performed.

| Bound | Maximum |
| --- | ---: |
| Declared width | 1,048,576 bits |
| Raw on-bit IDs, before duplicate removal | 1,048,576 |
| Conversion's temporary bit storage | 1 MiB |

Zero width is valid. At most one excess ID is consumed to detect the raw-input
cap. Bounds exclude caller-owned data and work performed inside arbitrary source
callbacks; those callbacks have no timeout.

Invalid shapes, ranges, writability, overlap or limits raise `ValueError`.
Wrong destination/index/protocol types raise `TypeError`; original source-call
and iterator exceptions propagate. A missing method is a protocol type error.
The existing `dest` keyword, optional destination and returned array are retained;
these differ from RDKit's in-place-only signature/return and resize behavior.
No whole-wrapper parity is claimed.

## Evidence

Two new regressions fail against the exact main Python source at `f20304f4`:
dirty destinations retain 7s, and a source exception leaves a partially written
1. A fresh independent preparation buffer fixes both before the final write.
The new implementation lives in a private Python/NumPy helper; no Rust code,
fingerprint generation or native binding inventory changes are involved.

The original probe at `validation/rdkit_issues/numpy-fingerprint/probe.py` checks
576 datasets (eight widths, three densities, four seeds, six destination dtypes).
11,520 adapter calls cover five source forms, three destination layouts and
optional output. All match an independent bit-set reference; 576 ordinary,
writable contiguous RDKit conversions agree as well. Invalid/read-only/strided
destinations are not passed into the upstream oracle. Supported adapter views
are compared with expected values, not with RDKit's view-handling behavior.

Source annotations are checked with valid and invalid mypy call examples;
this is an API input-type check, not a claim of complete static typing.
Versions, original inputs/results, failing baseline, hashes and logs are retained
at `validation/results/numpy-fingerprint-2026-10-17/`. A prebuilt extension fixture
was reused to import the package; the changed conversion does not call chematic
Rust APIs. No fresh native build, wheel validation or publication is claimed.
