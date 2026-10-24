# MolJSON non-node hydrogen counts

The existing chematic atom-record MolJSON reader now retains `hydrogens` as
chemical data. For example, an atom record describing carbon with three
hydrogens yields `[CH3]` after reading, instead of becoming methane (`C`).
Aromatic `[nH]` and a declared zero H count also survive a read/write cycle.

```python
import chematic

source = '{"atoms":[{"id":"a1","element":"C","hydrogens":3}],"bonds":[]}'
mol = chematic.from_moljson(source)
assert mol.smiles == "[CH3]"
assert chematic.from_moljson(mol.to_moljson()).smiles == "[CH3]"
```

A present `hydrogens` value must be a JSON integer from 0 through 255, the
core's storage range. Zero is an authoritative zero, not a request to infer
hydrogens. The field counts H atoms that are **not** separate graph nodes;
explicit H nodes are additional. For instance, a carbon with `hydrogens: 3`
and one bonded H node represents four hydrogens in total. When the field is
absent, existing valence-based inference remains in effect. Present null,
booleans, strings, fractional values and out-of-range counts fail with a
`hydrogens` error; they are not silently defaulted or narrowed.

The shared native reader supplies Rust, Python `from_moljson`, WASM
`mol_from_moljson`, common MolJSON conversion routes, and CLI `convert`.
No new public export or error-enum variant is added. This changes acceptance
of invalid H fields and interpretation of counts that were previously ignored.

## Dialect boundary

This is the existing chematic **atom-record dialect**, with per-atom `id`,
`element`, `charge`, `isotope`, `hydrogens` and `aromatic` fields. It differs
from the [current oxpig MolJSON schema](https://github.com/oxpig/MolJSON#format-notes),
which uses separate `charges` and `aromatic_n_h` lists. This change does not
implement that upstream schema or silently translate between the two.

This is an H-authority correction, not a strict validator of every MolJSON
field or a general lossless writer. The storage range does not establish
chemical valence validity. Spin multiplicity and stereochemistry are outside
this contract. The later [graph input integrity](moljson-graph-input.md) change
validates other known graph fields and corrects quadruple-bond reading; the
H-specific measurements below retain their original dated scope.

## Local validation — 2026-10-24

Thirteen original achiral graph inputs cover methane and ammonium controls,
methyl, zero-H carbon/oxygen, hydroxyl, amino, isotopic methylene, an anion,
pyrrole, one/four explicit H nodes and an omitted-field compatibility control.
The expected non-node H arrays and reference SMILES are authored in the packet.

| Python and actual Node/WASM, each | Previous local runtime | Candidate |
|---|---:|---:|
| Exact H array and both direct/common SMILES graphs | 5/13 | 13/13 |
| Invalid present H values refused | 0/9 | 9/9 |

RDKit 2026.09.1 observes each resulting SMILES with sanitization enabled and H
removal disabled. Its canonical graph string, formula, atom count and radical
electron count must match the authored reference SMILES. It is an independent
SMILES observer, not a MolJSON parser. Positive graph preservation and invalid
input refusals are separate results. Native/Python/Node tests also exercise
three successive read/write cycles; CLI tests cover direct/common conversion.

Open Babel 3.2.1's installed format inventory advertises no MolJSON format, so
it receives no same-format preservation score here. This packet closes a
measured local semantic gap; it does not establish Open Babel-wide superiority,
upstream MolJSON compatibility, an independent holdout or a published release.

Reproduction and compressed raw records are in
`validation/moljson_hydrogens/README.md` and
`validation/results/moljson-hydrogens-{baseline,candidate}-2026-10-24.json.gz`.
