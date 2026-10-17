# Open Babel 3.2.1: two-vacancy non-tetrahedral SMILES

## Scope

This source-built development comparison reuses the synthetic two-vacancy
corpora and CheMatic fix described in the [RDKit/COSMolKit report](../cosmolkit_comparison/two-vacancy-nontetra-2026-10-17.md).
RDKit 2025.09.3 accepts every source and supplies the canonical isomeric
SMILES identity. It then rereads each native writer's output. This measures
recognition and stereo preservation within these inputs, with RDKit as the
oracle; it is not an independent chemical truth set.

The identical source rows were passed to CheMatic and Open Babel. Each Open
Babel row used a separate `obabel -:SMILES -osmi` process: Open Babel stops
at the first rejected row of a multiline SMILES stream. A tetrahedral
control, `F[C@](Cl)(Br)I`, passed through the same command with stereo intact.

| Input | Rows | CheMatic exact | Open Babel exact | Open Babel outcome |
| --- | ---: | ---: | ---: | --- |
| TB3 and OH4, two vacant sites | 1,680 | **1,680** | 0 | 1,680 rejected |
| SP2, two vacant sites | 12 | **12** | 0 | 12 read with stereo lost |
| **Total** | **1,692** | **1,692** | **0** | |

The TB/OH rows have three or four distinct halogen ligands with all supported
stereo tags, ligand orders, and root/incoming positions. Open Babel returned
no molecule for these rows. On all SP rows it returned a molecule with the
same graph after removing stereo, but no SP assignment. These are support
and identity-preservation outcomes on this grammar, not a ranking of general
chemical accuracy. The corpus was built after inspecting an earlier
one-vacancy fixture, so it is development data rather than a sealed holdout.

## Reproduce

The [row-level packet](results/two-vacancy-nontetra-2026-10-17.json.gz)
contains every source, expected identity, native output, classification,
tool version, and input/binary/script hash. Its compressed SHA-256 is
`745621145f238a6b886342f3cb67337e42c19f561a72761d8c7a8f7b7e5471fe`.
The CheMatic example binary SHA-256 is
`36be62109c362340a678bacd60894eba54a83af8619a3982cf1e3d230710eb4e`;
the Open Babel binary SHA-256 is
`733b2c925d60b17a5c9d8ce9ec841e918cc20c7634301e23f9008a187d1d9bc5`.
The script SHA-256 is
`7f53f44a4bddb4289e31e84e6e52970d97680fee8fcb5515f52144a931a74d9a`.
The TB/OH and SP corpus SHA-256 values are respectively
`cc357f19953c026eba2f0f9811c6c446e96507395333f5a82f0041b46893f33e`
and `ddb9bb77793b7ce6e68ab36c2beb2d83adf4b037c6b6296f6c16a0ac4195d309`.

From the repository root with RDKit 2025.09.3 and Open Babel 3.2.1:

```sh
cargo build --offline -p chematic-smiles --example canonical_roundtrip_dump
python3 validation/openbabel_comparison/check_two_vacancy_smiles.py \
  --corpus validation/cosmolkit_comparison/corpora/two-vacancy-nontetra-v1.smi \
  --corpus validation/cosmolkit_comparison/corpora/two-vacancy-square-planar-v1.smi \
  --chematic target/debug/examples/canonical_roundtrip_dump \
  --openbabel /opt/homebrew/bin/obabel \
  --output /tmp/two-vacancy-openbabel.json.gz
```

The source-built result establishes a measured lead over Open Babel 3.2.1
on this narrow SMILES stereo task. It does not establish published-package
behavior or toolkit-wide superiority.
