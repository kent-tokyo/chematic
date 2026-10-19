# Checked RDKit-compatible SMILES

`chematic.rdkit_smiles_checked(mol)` adds an opt-in idempotence check to the
default RDKit-compatible canonical writer. It writes a candidate, reads that
candidate, and writes it once more. Only identical strings are returned.

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
| `non_idempotent` | The two strings differ | `first_smiles` and `roundtrip_smiles` |

Read/rewrite exceptions are retained as `__cause__`. Diagnostic strings are
rejected candidates, not verified identifiers. There is one readback and one
rewrite; later convergence is not accepted.

Textual idempotence does **not** prove chemical identity preservation,
stereochemical correctness, or canonical uniqueness across equivalent inputs.
A writer can lose information and still produce a stable string. In particular,
this check does not establish that legacy stereo perception retained every
stereo center specified by the caller. It is a guard against one observable
failure mode, not a complete resolution of #9391.

The implementation is Python orchestration over the native parser and writer;
it does not import RDKit. Regression tests include the reported ammonium
example, an explicitly tagged spiro junction, a carbon analogue, stable
spiro/cage/ring controls, and distinct tetrahedral and E/Z isomers.
