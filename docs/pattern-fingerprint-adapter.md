# Typed Pattern fingerprint adapter

`chematic.rdkit_compat.PatternFingerprint` exposes the existing native
concrete-molecule Pattern calculator through checked RDKit-style arguments.
This covers an input-typing slice of
[RDKit #9601](https://github.com/rdkit/rdkit/issues/9601), not nanobind or the
complete tracking issue.

```python
import numpy as np
from chematic import rdkit_compat as Chem

mol = Chem.MolFromSmiles("CCO")
counts = np.array([10, 20, 30], dtype=np.uint32)
empty_mask = Chem.ExplicitBitVect(128)
fp = Chem.PatternFingerprint(
    mol, fpSize=128, atomCounts=counts, setOnlyBits=empty_mask
)
assert fp.GetNumBits() == 128
assert fp.GetNumOnBits() > 0
assert counts.tolist() == [10, 20, 30]
```

## Counter and mask semantics

In RDKit 2026.09.1, the Pattern implementation checks `atomCounts` and
`setOnlyBits` sizes but does not use them when generating bits. A nonempty seed
counter remains unchanged and an empty mask does not suppress output bits.
This adapter preserves that behavior explicitly. It does not invent atom
contribution counts or turn this argument into a filtering operation.
`LayeredFingerprint` has different, documented counter semantics.

`PatternAtomCounts` is a plain `list[int]` with u32-valued entries or a plain
writable one-dimensional NumPy `uint32` array. Nonempty destinations must cover
all atoms and have at most 65536 entries. Empty destinations are accepted.
Reversed, strided and either-byte-order views are supported; overlapping,
read-only or subclass arrays are refused. Counts are validated but never written.

`PatternBitMask` is a structural protocol requiring `GetNumBits()` and
`GetOnBits()`. Width must equal `fpSize`; on-bit IDs must be unique integers in
range. Explicit, sparse and generator-backed sources are accepted and copied
for validation, then discarded. Source exceptions propagate. User-defined
source side effects are not rolled back.

## Supported calculation

- `mol` must be a concrete compatibility `Mol`, with at most 1024 atoms and
  4096 bonds. SMARTS query objects, native `Mol` objects and bundle overloads
  are rejected by this facade; native molecules retain their existing methods.
- `fpSize` must be one of 1, 2, 4, 8, 16, 32, 64, 128, 256, 512, 1024 or
  2048. Integer-like NumPy values are accepted; boolean/float widths are not.
- `tautomerFingerprints` must be `False`. Unsupported modes fail explicitly.
- Options are validated before the native fingerprint calculation. Results
  are independent mutable `ExplicitBitVect` objects.

The native calculator hashes each feature modulo 2048. For any supported
width `w` dividing 2048, `(hash % 2048) % w == hash % w`, so folding produces
the same width reduction. This identity does not justify arbitrary widths or
widths larger than 2048, which RDKit supports but this adapter rejects.

This reuses the existing native Pattern algorithm and its aromaticity/model
compatibility limits. It introduces no new feature definitions or query/tautomer
algorithm. Graph limits do not create a new matcher visit budget or deadline.

## Validation

76 molecules, 12 widths and three destination/mask profiles give 2736 complete
bit-vector comparisons against RDKit 2026.09.1, all matching. The profiles with
supplied counters additionally check 1824 unchanged counter/mask cases.
There are 54 new tests; 377 tests pass including existing fingerprint,
Layered-counter and compatibility tests. Reproducible source/runtime boundaries
are in `validation/results/typed-pattern-fingerprint-2026-10-19/`.
