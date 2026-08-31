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

`similarity <SMILES_A> <SMILES_B>` emits the selected algorithm and Tanimoto
similarity. It accepts the same `ecfp4`, `ecfp6`, and `maccs` algorithms as
`fingerprint`.

`substructure <SMILES> <SMARTS>` emits canonical SMILES, the query, match
count, and query-order atom-index mappings as JSON. Invalid SMARTS is reported
as a command error.

`standardize <SMILES>` applies the default standardization pipeline and emits
the canonical input and output SMILES plus an audit report. The report includes
the overall status, whether the structure changed, each pipeline stage's
enabled/changed flags and before/after snapshots, and machine-readable
warnings. This makes normalization decisions inspectable in shell pipelines.

```sh
chematic standardize 'C[NH3+]'
```

`report <SMILES>` emits the complete single-molecule analysis record used by
the workflow API: canonical identity, formula, Murcko scaffold, descriptor
summary, drug-likeness filters and alert names, functional groups, and named
groups. Invalid input is rejected with a non-zero exit status.

```sh
chematic report 'CC(=O)Oc1ccccc1C(=O)O'
```
