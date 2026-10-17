# Two-vacancy non-tetrahedral SMILES recognition

## Scope

This is a source-built, synthetic development comparison for SMILES with two
unoccupied coordination sites. It does not measure general toolkit accuracy.
RDKit supplies the expected isomeric SMILES identity and independently reads
each toolkit's native output. RDKit is therefore the oracle, not a competing
native writer in this experiment. [RDKit's documented convention](https://rdkit.org/docs/RDKit_Book.html#missing-ligands)
places missing ligands after the incoming ligand in SMILES.

`corpora/two-vacancy-nontetra-v1.smi` contains all permutations of three
distinct TB ligands (F, Cl, Br) with every `@TB1`–`@TB20` tag, and all
permutations of four distinct OH ligands (F, Cl, Br, I) with every
`@OH1`–`@OH30` tag. Each appears with the center both first and after an
incoming ligand: 240 TB and 1,440 OH rows. RDKit 2025.09.3 accepts all 1,680.
The corpus SHA-256 is
`cc357f19953c026eba2f0f9811c6c446e96507395333f5a82f0041b46893f33e`.
The separate 12-row `corpora/two-vacancy-square-planar-v1.smi` covers both
orders of F/Cl and every `@SP1`–`@SP3` tag, again with root and incoming
positions. Its SHA-256 is
`ddb9bb77793b7ce6e68ab36c2beb2d83adf4b037c6b6296f6c16a0ac4195d309`.

## Result

| Native writer | Before | After | Strict RDKit read-back |
| --- | ---: | ---: | ---: |
| CheMatic | 92 / 1,680 | **1,680 / 1,680** | 1,680 / 1,680 |
| COSMolKit 0.5.0 | 240 / 1,680 | 240 / 1,680 | 240 / 1,680 |

On the separate square-planar set, CheMatic matches **12 / 12** and
COSMolKit **4 / 12**. The [SP row-level packet](results/two-vacancy-square-planar-candidate-2026-10-17.json.gz)
has compressed SHA-256
`16b295e453dd8126e257fb9f4623781e64fb7bf9e5c136f69a92e608f0161b5e`.

The before CheMatic binary was built from `origin/main` at `35b932f9`; its
SHA-256 is `2df718671c765f2d4c5f05280ae3a14d19cc48b567b297e1cc4689be9be614f9`.
The after binary SHA-256 is
`36be62109c362340a678bacd60894eba54a83af8619a3982cf1e3d230710eb4e`.
The [before packet](results/two-vacancy-nontetra-main-2026-10-17.json.gz)
has compressed SHA-256
`97deb22dd1e1c00a6be6b950847a2f6e6055e1a2faae01d70e23feefa4292dd5`;
the [after packet](results/two-vacancy-nontetra-candidate-2026-10-17.json.gz)
has compressed SHA-256
`77b47f22e959ae00ae9bdaa572087e80ebdb8ad2cc68e3f51ca4c78d8797a950`.
Both packets contain row-level classifications and all unresolved indices.

The fix stores both empty ligand positions and exchanges their order when the
center moves between root and incoming-ligand SMILES positions. Full and
one-vacancy centers retain their previous representation; centers missing
more than two sites remain outside this change.

For an independent traversal check, `random_roundtrip_dump` emitted 16
randomized SMILES per source. RDKit 2025.09.3 reread **26,880 / 26,880** TB/OH
and **192 / 192** SP outputs with the same canonical isomeric identity. The
[TB/OH row-level dump](results/two-vacancy-nontetra-random-16-2026-10-17.tsv.gz)
has compressed SHA-256
`b2c4ccc410ef16ba5351e3f48201977963177d282a8b49b3381ed3aa6dda4978`;
its uncompressed SHA-256 is
`c8a1ed58c23c6aad6f4977fc048c93f502b664078e0a7d0d22a65d05ab1efc5a`.
The source-built random example binary SHA-256 is
`1300bc4070c89ba284e236cf07b9fe3e0e8e5b8331704caf0fea78358453d24b`.
The [SP random dump](results/two-vacancy-square-planar-random-16-2026-10-17.tsv.gz)
has compressed SHA-256
`c348074fdc90491e2202290d62c20a7b3a2089c5f1db511d39cf8edf6d3ca42e`.

RDKit 2026.09.1 separately accepted all 1,680 sources and confirmed the
same CheMatic 1,680 and COSMolKit 240 strict identities. The prior one-vacancy
4,080-row corpus remained 4,080 / 4,080 for both toolkits. On the separate
4,072-row stress corpus, CheMatic retained 4,072 verified-equivalent
read-backs (4,050 strict plus 22 classified equivalent); COSMolKit retained
4,064 strict read-backs. The added Rust stereo regression, all
`chematic-smiles` tests, strict Clippy, format, and refactor-budget checks
passed locally.

## Reproduce

From the repository root, with RDKit 2025.09.3 and COSMolKit 0.5.0 installed
in the same Python environment:

```sh
python3 validation/cosmolkit_comparison/corpora/make_two_vacancy_nontetra_corpus.py \
  --output /tmp/two-vacancy.smi \
  --square-planar-output /tmp/two-vacancy-sp.smi
cargo build --offline -p chematic-smiles \
  --example canonical_roundtrip_dump --example random_roundtrip_dump
python3 validation/cosmolkit_comparison/check_native_smiles_roundtrip.py \
  --corpus /tmp/two-vacancy.smi \
  --chematic target/debug/examples/canonical_roundtrip_dump
target/debug/examples/random_roundtrip_dump /tmp/two-vacancy.smi 16 \
  > /tmp/two-vacancy-random.tsv
python3 validation/cosmolkit_comparison/check_random_smiles_stereo.py \
  --corpus /tmp/two-vacancy.smi --dump /tmp/two-vacancy-random.tsv
# Repeat the two checks with /tmp/two-vacancy-sp.smi.
```

The corpus was designed after inspecting a one-vacancy development fixture.
The 1,680/1,680 result is a measured lead over COSMolKit 0.5.0 on this
specific input and oracle. It is not a held-out estimate or proof of complete
superiority over RDKit or COSMolKit.
