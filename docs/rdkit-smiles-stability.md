# Checked RDKit-compatible SMILES

`chematic.rdkit_smiles_checked(mol)` adds two bounded checks to the default
RDKit-compatible canonical writer. It first writes, reads and rewrites one
candidate. It then renumbers the same graph in up to 32 deterministic orders
and requires every canonical string to be identical.

```python
import chematic as ch

assert ch.rdkit_smiles_checked(ch.from_smiles("OCC")) == "CCO"

mol = ch.from_smiles("C[N@+]1(F)CCC2(CC[C@H](Cl)CC2)CC1")
try:
    identifier = ch.rdkit_smiles_checked(mol)
except ch.RdkitSmilesStabilityError as error:
    assert error.code == "non_idempotent"
    print(error.first_smiles)
    print(error.roundtrip_smiles)
    print(error.atom_order)
```

[RDKit issue #9391](https://github.com/rdkit/rdkit/issues/9391) reports a spiro
compound whose legacy canonical SMILES alternates between `@` and `@@` on
successive reads and writes. The same oscillation is reproducible in
chematic's compatibility writer. This API detects and refuses the unstable
candidate. It does not repair the underlying stereo perception or select a
member of the cycle. Existing `Mol.rdkit_smiles` behavior is unchanged.

The check uses the default canonical, isomeric writer. It accepts a `Mol`,
leaves that molecule unchanged, and returns a string. Initial writer errors
propagate. Failures after obtaining a candidate raise
`RdkitSmilesStabilityError`, a `ValueError` subclass:

| `code` | Meaning | Diagnostic fields |
|---|---|---|
| `unreadable_output` | The candidate cannot be read | `first_smiles`; `roundtrip_smiles=None` |
| `rewrite_failed` | Reading succeeded, but rewriting failed | `first_smiles`; `roundtrip_smiles=None` |
| `non_idempotent` | The read/write fixed point differs | `first_smiles` and `roundtrip_smiles` |
| `permutation_write_failed` | A renumbered graph cannot be written | `first_smiles`; `atom_order` |
| `permutation_dependent` | A renumbered graph produces another string | both strings; `atom_order` |

Read/rewrite exceptions are retained as `__cause__`. Diagnostic strings are
rejected candidates, not verified identifiers. There is one readback and one
rewrite; later convergence is not accepted.

`max_permutations` includes the identity order, defaults to 32 and accepts
values from 2 through 256. The set includes reversal, rotations, even/odd
partitions and deterministic shuffles. It never grows factorially.

Fixed-point plus bounded atom-order invariance does **not** prove chemical
identity preservation, stereochemical correctness, or canonical uniqueness
over every equivalent input. A writer can lose information and still pass
both checks, and larger graphs have unsampled orders. In particular, this
check does not establish that legacy stereo perception retained every center
specified by the caller. It is a fail-closed guard against two observable
failure modes, not a complete resolution of #9391.

The implementation is Python orchestration over the native parser and writer;
it does not import RDKit. Regression tests include the reported ammonium
example, an explicitly tagged spiro junction, a carbon analogue, stable
spiro/cage/ring controls, distinct tetrahedral and E/Z isomers, deterministic
order generation and injected permutation failures.
