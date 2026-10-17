# Cancellable tautomer enumeration

`enumerate_tautomers_checked` runs chematic's Rust RDKit-profile tautomer
algorithm and calls a typed Python callback during enumeration. A callback can
inspect the current keys and stop before another transform is applied. This
addresses the callable-boundary and cancellation use case tracked in
[RDKit #9601](https://github.com/rdkit/rdkit/issues/9601).

```python
import chematic

mol = chematic.from_smiles("CC(=O)CC=O")

def keep_going(progress: chematic.TautomerProgress) -> bool:
    return len(progress.keys) < 3

result = chematic.enumerate_tautomers_checked(
    mol, max_tautomers=1000, max_transforms=1000, callback=keep_going,
)
assert result.status == "Canceled"
print(result.smiles)
```

The callback receives a frozen `TautomerProgress` with `transforms: int` and
`keys: tuple[str, ...]`. The count increases once for each transform rule with
a nonempty match set; multiple checkpoints can share the same count. Keys are
snapshots of the current enumeration map in lexicographic order. Earlier
snapshots remain unchanged as enumeration continues.

Return an actual Python `bool`: `True` continues and `False` cancels. Callable
instances work too. An invalid return type raises `TypeError`; an exception
raised by the callback propagates unchanged. Calls are synchronous on the
calling thread. A molecule with no applicable transform produces no callback.

`TautomerEnumeration` is frozen and contains:

| Field | Meaning |
| --- | --- |
| `smiles: tuple[str, ...]` | Final canonical SMILES, in enumeration-key order |
| `keys: tuple[str, ...]` | Final enumeration map keys |
| `status` | `Completed`, `MaxTautomersReached`, `MaxTransformsReached`, or `Canceled` |

Cancellation stops additional transforms. The algorithm still normalizes stereo
on the partial result before returning it. Progress keys, final keys and final
SMILES can consequently differ. Input molecule graph state is unchanged.

## Limits and scope

Both limits default to 1000 and accept integers from 1 through 10000; bools are
rejected. Python inputs are capped at 1024 graph atoms. Limit checks precede
callbacks: the transform limit wins over the tautomer limit when both are
reached. The cutoff follows the existing RDKit-profile algorithm; for example,
`max_transforms=1` stops at the first applicable transform before applying it.

Callbacks run before individual matched transforms, after SMARTS matching.
They do not preempt a match in progress or impose a wall-clock timeout. This
function accepts synchronous callbacks; asynchronous callbacks do not satisfy
the required bool return contract. Enumeration retains the existing concrete
molecule, sanitization and stereo scope of the Rust RDKit-profile implementation.
This API is distinct from RDKit's `(mol, result)` callback and `SetCallback`
class interface. No RDKit installation is required at runtime.

Rust callers can use
`chematic_smiles::rdkit_enumerate_tautomers_with_callback(mol, max_tautomers,
max_transforms, &mut callback)`. It validates positive limits and supplies
`RdkitTautomerProgress`; the Python-specific upper limits are enforced by the
Python entry point. Existing default enumeration and canonicalization APIs
retain their signatures.

## Validation

A fresh native Python extension was built from this implementation. The checked-in
[oracle runner](https://github.com/kent-tokyo/chematic/blob/main/validation/validate_tautomer_callback.py)
compares 24 concrete molecules under 12 cancellation/limit policies against
RDKit 2026.09.1. All 288 conditions matched final SMILES, keys, status and the
sequence of checkpoint keys (557 chematic callback calls). The cases include
empty graphs, no-match controls, salts, keto/enol, imine, nitro, phosphorus,
aromatic heterocycles and stereo-bearing molecules. These observations do not
establish parity for every molecule or RDKit option.

The Python regression suite adds 62 cases for cancellation, frozen snapshots,
limit precedence, typed inputs/outputs, reentrant calls and exception identity.
The existing 162 compatibility tests also pass. The Rust library's full test
suite covers the refactoring and its added callback tests. Evidence, commands
and native build provenance are retained under
`validation/results/tautomer-callback-2026-10-17/`.
