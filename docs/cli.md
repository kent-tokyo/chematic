# Command-line interface

`chematic-cli` publishes the `chematic` binary for small, composable molecule
conversion workflows. It reads stdin and writes stdout by default, so it can be
used in shell pipelines without a temporary file.

```sh
printf 'CCO\n' | chematic convert \
  --input-format smiles --output-format mol2

chematic convert --input-format mol2 --output-format smiles \
  --input ligand.mol2 --output ligand.smi
```

```sh
chematic descriptors 'CC(=O)Oc1ccccc1C(=O)O'
```

The initial CLI contract covers topology-bearing SMILES, MOL/SDF (V2000), MOL
V3000, MOL2, CML, ChemicalJSON, MolJSON, and CDXML. Format aliases such as
`smi`, `.mol2`, and `sdf` are accepted. Coordinates and format-specific
metadata are intentionally outside this first command; use the format-specific
Rust, Python, or WASM APIs when they must be preserved.

Errors are written to stderr and return exit status 2. Unsupported formats and
parse failures are explicit errors; the command never silently forwards the
input unchanged.

`descriptors` emits compact JSON containing canonical SMILES, formula, atom
counts, molecular weight, exact mass, Crippen logP, TPSA, HBD, HBA, and
rotatable-bond count.

`fingerprint <SMILES>` emits JSON with the selected algorithm (`ecfp4`,
`ecfp6`, or `maccs`), width, set-bit indices, popcount, and canonical SMILES.
Set-bit indices
are used instead of an opaque binary string so the output is easy to inspect
and remains language-neutral.
