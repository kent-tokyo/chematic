# Typed caller-owned Layered counters

`chematic.rdkit_compat.LayeredFingerprint` exposes the existing checked native
Layered calculator through RDKit-style argument names and updates an explicit
caller-owned counter destination. It addresses the `atomCounts` output-argument
boundary in [RDKit #9601](https://github.com/rdkit/rdkit/issues/9601).
The underlying hashing and enumeration implementation is unchanged.

```python
import numpy as np
from chematic import rdkit_compat as Chem

mol = Chem.MolFromSmiles("CCO")
counts = [10, 20, 30, 99]
fp = Chem.LayeredFingerprint(mol, atomCounts=counts)
assert counts == [12, 23, 32, 99]
assert fp.GetNumBits() == 2048

array = np.zeros(mol.GetNumAtoms(), dtype=np.uint32)
Chem.LayeredFingerprint(mol, atomCounts=array)
assert array.tolist() == [2, 3, 2]
```

## Typed destinations and update order

`LayeredAtomCounts` is `list[int] | ndarray[..., dtype[uint32]]`.
Only a plain built-in list or a plain NumPy ndarray is accepted at runtime.
List entries must be nonnegative index integers fitting u32; bool is refused.
Arrays must be one-dimensional, writable and uint32. Native/non-native byte
order, ordinary strides and reversed views are supported; overlapping elements,
array subclasses, masked arrays, float/signed/other-width arrays, tuples and
custom mutation callbacks are refused.

A supplied destination must cover the stored atom count and contain no more than
65,536 entries. Counts accumulate from the supplied seeds with u32 wrapping.
Extra trailing entries remain unchanged. The list/array object and an array's
shape, dtype and backing storage are retained. `None` means no counter update.
A short or empty destination for a nonempty molecule is an error; RDKit's Boost
wrapper treats an empty list as an absent output request, so that case differs
explicitly.

The adapter copies and checks the complete seed state, runs the native calculation,
and constructs the returned fingerprint before writing any caller counters.
Input validation failures, mask iteration errors, unsupported graphs, native
budget exhaustion and fingerprint construction failures leave the destination
unchanged by the adapter. The final complete result is staged before a single
list slice or plain NumPy view assignment. The API provides no synchronization
for concurrent mutations of the same buffer.

The underlying `Mol.rdkit_layered_config(atom_counts=...)` remains a separate
copy-in/copy-out API; it does not mutate its input list.

## Options and graph scope

The positional/keyword parameters follow RDKit's order: `mol`, `layerFlags`,
`minPath`, `maxPath`, `fpSize`, `atomCounts`, `setOnlyBits`, `branchedPaths`,
`fromAtoms`. `maxVisits` and `maxPaths` are additional keyword-only native budgets.
Integral scalars support `__index__`, including NumPy integers, and exclude bool.
`branchedPaths` requires bool. `fromAtoms` accepts a bounded indexed sequence of
integer atom indices; generators and NumPy arrays are outside this sequence
profile. Empty roots use the native unrooted branched calculation. Root duplicates
retain their native multiplicity.

`LayeredBitMask` is a read-only protocol with `GetNumBits()` and `GetOnBits()`;
explicit and sparse binary vectors can provide it. Its width must equal `fpSize`.
Invalid/duplicate set-bit indices and mismatched widths are refused. The mask is
copied without mutation. Masking affects accepted paths and atom counts, in
addition to the returned bits.

See [native Layered options](rdkit-layered-fingerprint-options.md) for the full
hashing, path, aromaticity, explicit-H and ring semantics. This adapter inherits
those limits: 512 atoms, 1,024 bonds, 1–65,536 bits, paths of 1–12 bonds,
4,096 roots, and bounded search. Unrooted `branchedPaths=False` is refused;
rooted linear paths are supported. Graphs needing index-changing preparation,
complex query/dative/zero bonds, loops or parallel edges are outside the native
profile. Native typed calculation errors retain their category/code/format.

This addition does not implement PatternFingerprint atom counts, nanobind
wrappers or other unchecked portions of the upstream tracking issue.

## Evidence

`validation/results/layered-counts-2026-10-17/` records the issue snapshot,
source/native artifact hashes, test/type logs and the live comparison.
Local tests passed 307 cases, including 81 new adapter contracts. Across 16
original molecules and 2,048 option sets, output bits and seeded counters match
RDKit 2026.09.1 exactly. Another 2,048 second calls check accumulation, and 4,096
strided/reversed NumPy checks retain both result counters and untouched storage.
This is measured parity on the stated corpus, with the native experimental
Layered profile's limitations preserved.

The helper passes strict mypy in isolation. A public consumer accepts list/u32
array output and rejects four representative wrong argument types. No
package-wide strict typing or fresh native build is claimed: the documented
extension fixture reuses unchanged relevant native sources.
