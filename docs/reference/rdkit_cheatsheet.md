# RDKit to chematic cheat sheet

Use this page to translate common Python workflows. chematic is not a drop-in
replacement for every RDKit API: choose the named compatibility method when
exact RDKit behavior matters, and review the [migration matrix](rdkit-migration.md)
for partial or unsupported areas.

```bash
pip install chematic
```

## Molecules and identifiers

| RDKit | chematic | Notes |
|---|---|---|
| `Chem.MolFromSmiles(s)` | `chematic.from_smiles(s)` | Invalid input raises a typed Python error |
| `Chem.MolToSmiles(m)` | `m.smiles` | Canonical chematic output |
| `Chem.MolFromMolBlock(b)` | `chematic.from_mol_block(b)` | Use checked/source-preserving APIs when loss reporting matters |
| `Chem.inchi.MolFromInchi(s)` | `chematic.from_inchi(s)` | Standard InChI support |
| `Chem.inchi.MolToInchi(m)` | `m.inchi` |  |
| `Chem.inchi.MolToInchiKey(m)` | `m.inchikey` |  |
| `SDMolSupplier(path)` | `chematic.iter_sdf(path)` | Streaming records |

## Descriptors and filters

| RDKit | chematic |
|---|---|
| `Descriptors.MolWt(m)` | `m.mw` |
| `Descriptors.ExactMolWt(m)` | `m.exact_mass` |
| `Descriptors.MolLogP(m)` | `m.logp` |
| `rdMolDescriptors.CalcTPSA(m)` | `m.tpsa` |
| `Descriptors.NumHDonors(m)` | `m.hbd` |
| `Descriptors.NumHAcceptors(m)` | `m.hba` |
| `Descriptors.NumRotatableBonds(m)` | `m.rotatable_bonds` |
| `Descriptors.FractionCSP3(m)` | `m.fsp3` |
| `Descriptors.RingCount(m)` | `m.ring_count` |
| `Descriptors.qed(m)` | `m.qed` |
| `rdMolDescriptors.CalcMolFormula(m)` | `m.formula` |

Common filters are exposed as properties such as `m.lipinski_passes`,
`m.veber_passes`, `m.pains_passes`, and `m.brenk_passes`. `pka()` and
`admet()` are rule-based screening helpers; they are not experimental or
clinical prediction models.

## Fingerprints and similarity

| RDKit workflow | chematic |
|---|---|
| Native circular fingerprint | `m.ecfp4()`, `m.ecfp6()`, `m.fcfp4()` |
| RDKit-compatible Morgan bits | `m.rdkit_ecfp4_numpy()` or the named packed-byte API |
| MACCS | `m.maccs()` |
| Atom pair | `m.atom_pair_fp()` |
| Topological torsion | `m.torsion_fp()` |
| Layered | `m.layered_fp()` |
| Tanimoto | `chematic.tanimoto(fp1, fp2)` |

Do not infer bit identity from the generic `ecfp4()` name. The native
fingerprint has its own documented definition; use an `rdkit_*` method when
RDKit-compatible bits are required.

```python
smiles = ["CCO", "c1ccccc1", "CC(=O)O"]

X = chematic.bulk.ecfp4(smiles)            # uint8, shape (N, 2048)
S = chematic.bulk.tanimoto(smiles, smiles) # float32, shape (N, N)

index = chematic.SimilarityIndex.from_smiles(smiles)
hits = index.search("CCO", threshold=0.3, k=10)
```

## SMARTS, MCS, and R-groups

| RDKit | chematic |
|---|---|
| `m.GetSubstructMatches(q)` | `chematic.smarts_find(pattern, m)` |
| `m.HasSubstructMatch(q)` | `chematic.smarts_match(pattern, m)` |
| `rdFMCS.FindMCS(mols)` | `chematic.find_mcs(mols)` |
| `RGroupDecompose(...)` | `chematic.rgroup_decompose_labeled(...)` |
| R-group table | `chematic.rgroup_decompose_columns(...)` |

```python
matches = chematic.smarts_find("[OH]", mol)
mcs = chematic.find_mcs([mol1, mol2])

columns = chematic.rgroup_decompose_columns(
    "c1cc([*:1])ccc1[*:2]", mols
)
```

R-group support is a bounded mapped-label subset. Unsupported or ambiguous
contracts fail explicitly; RDKit's full core enumeration and GA scoring are
not implied.

## Reactions and standardization

```python
# RDKit: ReactionFromSmarts(...).RunReactants((mol,))
products = chematic.run_smirks("[OH:1]>>[O-:1]", [mol])

# RDKit: Cleanup + LargestFragmentChooser + Uncharger
clean = mol.standardize()

# Individual operations
largest = mol.largest_fragment()
neutral = mol.neutralize()
canonical = mol.canonical_tautomer()
tautomers = mol.enumerate_tautomers()
```

The chematic standardization pipeline is explicit and tested, but not every
RDKit cleanup policy is identical. Treat [compatibility scope](compatibility-scope.md)
as the contract.

## Scaffolds, fragments, and depiction

| RDKit | chematic |
|---|---|
| `MurckoScaffold.GetScaffoldForMol(m)` | `m.scaffold()` |
| `MurckoScaffold.MakeScaffoldGeneric(m)` | `m.generic_scaffold()` |
| BRICS decomposition | `m.brics_fragments()` |
| `Draw.MolToImage(m)` | `m.svg()` |
| highlighted drawing | `m.svg_highlighted(atoms, color=...)` |
| grid drawing | `chematic.depict_grid(mols, cols=...)` |

## 3D coordinates and force fields

chematic coordinate transforms return a new coordinate array instead of
mutating a conformer:

```python
coords = mol.set_bond_length(coords, 0, 1, 1.52)
coords = mol.set_bond_angle(coords, 0, 1, 2, 109.5)
coords = mol.set_dihedral(coords, 0, 1, 2, 3, 180.0)
transform = mol.compute_canonical_transform(coords)
coords = mol.canonicalize_conformer(coords)
```

Ring edits that cannot be represented safely are rejected. Platform-dependent
eigenvector signs are not a stable comparison target; compare geometry and
handedness instead. For 3D generation, MMFF/UFF, constraints, and failure
diagnostics, use the checked APIs described in [errors and limits](error-and-limits.md).

## Bulk DataFrame workflow

```python
import pandas as pd

df = pd.DataFrame(chematic.bulk.descriptors(smiles_list))
```

## Where chematic differs

| Area | Contract |
|---|---|
| Browser | Typed WASM/Node package and local execution; see [bindings](language-bindings.md) |
| AI tools | Optional local MCP server; see [MCP protocol](mcp-protocol-reference.md) |
| Rust | Rust-native core; native InChI is the documented optional C FFI exception |
| Screening | Built-in rule-based pKa/ADMET helpers, with explicit non-clinical limits |
| Failures | Checked APIs return typed refusal or diagnostics instead of inventing unsupported chemistry |

For API details, see the [cookbook](cookbook.md), [language bindings](language-bindings.md),
and [compatibility scope](compatibility-scope.md).
