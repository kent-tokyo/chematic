# Checked CJSON output: refuse unrepresented graph meaning

Rust `write_cjson_checked` and Python `Mol.to_cjson_checked` serialize the local
ChemicalJSON graph/Cartesian profile only when its in-memory properties can be
retained. Python `convert_format` and WASM `convert_common_format` use this guard
for graph-projected CJSON output; the CLI's topology writer uses the same guard.
This does not add coordinate retention to the CLI's topology bridge.
Exact CJSON retention and the separate MOL2
selected-field projection keep their own contracts.

The writer retains atom/bond indices, elements, formal charges, positive isotope
masses, single/double/triple/aromatic bond orders and optional finite Cartesian
coordinates in Angstrom. Coordinates must be absent or have one complete row
per stored atom, including explicit H nodes. Supported finite binary64 values,
including signed zero, roundtrip through the local reader.

The guard refuses zero, dative, query and directional bonds. It also refuses
quadruple bonds: this reader interprets numeric order `4` as an aromatic alias,
so writing a quadruple bond as `4` changes it on reparse. Atom chirality, explicit
H counts (including zero), aromatic atom flags, wildcards, maps, tags, R-group
labels, CIP descriptors, isotope mass zero, enhanced stereo groups and stashed
bond directions are refused. Self bonds outside the checked reader's profile
are also refused. No stereo inference, H expansion or kekulization
is performed. Errors identify the property and atom/bond index when applicable;
Rust provides typed errors; Python/WASM raise an exception with no candidate output.

```python
import chematic

mol = chematic.from_smiles("CCO")
text = mol.to_cjson_checked()  # topology; 3D remains absent

# Raises ValueError: cannot preserve bond order Up at bond 0.
chematic.convert_format("F/C=C/F", "smiles", "cjson")
```

This is intentionally conservative. Aromatic atom flags are separate from
aromatic bond orders, and explicit H declarations are separate from inferred
hydrogens. A molecule with either requires refusal even when a downstream
reader could infer an equivalent chemical structure. For example, chematic's
CML writer declares total H counts; CML → CJSON now refuses to omit them.
Use the legacy `Mol.to_cjson`/`write_cjson` only when accepting their lossy
projection. Same-format CJSON source retention remains exact.

On eight original CML graph records, Python/WASM each retained all three
supported graphs/XYZ and explicitly refused all five zero/query records.
Open Babel 3.2.1, in both default and `-xs` modes, retained the three positives
and emitted ordinary bond orders for the five unsupported records. These are
synthetic tests of a local reader profile; they do not establish CML specification
violations or toolkit-wide superiority. See the
[fixtures, complete results and reproduction instructions](https://github.com/kent-tokyo/chematic/tree/main/validation/cjson_output).

The guard starts after parsing. Source fields already ignored by a reader are
outside this check, as are metadata and cross-format stereo reconstruction.
The [ChemicalJSON format description](https://github.com/OpenChemistry/chemicaljson/blob/main/chemicaljson.md)
describes parallel atom, coordinate and bond arrays; the profile and refusals
above are chematic's explicit contract rather than a universal CJSON definition.
