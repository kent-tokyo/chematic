# RDKit → chematic migration cheatsheet

Side-by-side API reference for users coming from RDKit. chematic is pure Rust with zero C/C++ dependencies — install with `pip install chematic`, no conda required.

See also: [`rdkit-migration.md`](rdkit-migration.md) for a feature-by-feature Supported/Partial/Not-supported breakdown.

## Installation

```bash
# RDKit
conda install -c conda-forge rdkit   # or: pip install rdkit

# chematic
pip install chematic                  # no extra dependencies
```

---

## Basic I/O

| RDKit | chematic | Notes |
|-------|----------|-------|
| `Chem.MolFromSmiles(smi)` | `chematic.from_smiles(smi)` | SMILES → Mol |
| `Chem.MolToSmiles(mol)` | `mol.smiles` | Mol → canonical SMILES |
| `Chem.MolFromMolBlock(block)` | `chematic.from_mol_block(block)` | MOL block → Mol |
| `Chem.inchi.MolFromInchi(inchi)` | `chematic.from_inchi(inchi)` | InChI → Mol |
| `Chem.inchi.MolToInchi(mol)` | `mol.inchi` | Mol → InChI |
| `Chem.inchi.InchiToInchiKey(inchi)` | `mol.inchikey` | InChIKey |
| `Chem.MolToSmiles(mol)` | `str(mol)` | string conversion |

## Descriptors

| RDKit | chematic |
|-------|----------|
| `Descriptors.MolWt(mol)` | `mol.mw` |
| `Descriptors.ExactMolWt(mol)` | `mol.exact_mass` |
| `Descriptors.MolLogP(mol)` | `mol.logp` |
| `Descriptors.TPSA(mol)` | `mol.tpsa` |
| `Descriptors.NumHDonors(mol)` | `mol.hbd` |
| `Descriptors.NumHAcceptors(mol)` | `mol.hba` |
| `Descriptors.NumRotatableBonds(mol)` | `mol.rotatable_bonds` |
| `Descriptors.FractionCSP3(mol)` | `mol.fsp3` |
| `Descriptors.HeavyAtomCount(mol)` | `mol.heavy_atoms` |
| `Descriptors.RingCount(mol)` | `mol.ring_count` |
| `Descriptors.NumAromaticRings(mol)` | `mol.aromatic_ring_count` |
| `Descriptors.MolMR(mol)` | `mol.molar_refractivity` |
| `Descriptors.qed(mol)` | `mol.qed` |
| `rdMolDescriptors.CalcTPSA(mol)` | `mol.tpsa` |
| `rdMolDescriptors.CalcNumStereocenters(mol)` | `mol.num_stereocenters` |
| `rdMolDescriptors.CalcMolFormula(mol)` | `mol.formula` |
| `sascorer.calculateScore(mol)` | `mol.sa_score` |
| `Descriptors.ESOL(mol)` *(rdkit-contrib)* | `mol.esol` |

## Drug-likeness filters

| RDKit | chematic |
|-------|----------|
| `FilterCatalog.FilterCatalogParams(PAINS)` | `mol.pains_passes` |
| *(manual implementation)* | `mol.lipinski_passes` |
| *(manual implementation)* | `mol.veber_passes` |
| *(manual implementation)* | `mol.ghose_passes` |
| *(manual implementation)* | `mol.egan_passes` |
| *(manual implementation)* | `mol.reos_passes` |
| *(manual implementation)* | `mol.brenk_passes` |

## Fingerprints

| RDKit | chematic |
|-------|----------|
| `AllChem.GetMorganFingerprintAsBitVect(mol, 2, 2048)` | `mol.ecfp4()` |
| `AllChem.GetMorganFingerprintAsBitVect(mol, 3, 2048)` | `mol.ecfp6()` |
| `AllChem.GetMorganFingerprintAsBitVect(mol, 2, 2048, useFeatures=True)` | `mol.fcfp4()` |
| `AllChem.GetMACCSKeysFingerprint(mol)` | `mol.maccs()` |
| `rdMolDescriptors.GetAtomPairFingerprintAsBitVect(mol)` | `mol.atom_pair_fp()` |
| `rdMolDescriptors.GetTopologicalTorsionFingerprintAsBitVect(mol)` | `mol.torsion_fp()` |
| `rdMolDescriptors.GetHashedTopologicalTorsionFingerprintAsBitVect(mol)` | `mol.layered_fp()` |

Return values are `bytes` (chematic) or `DataStructs.ExplicitBitVect` (RDKit). Pass chematic bytes directly to `chematic.tanimoto()`.

## Similarity

```python
# RDKit
from rdkit.DataStructs import TanimotoSimilarity
sim = TanimotoSimilarity(fp1, fp2)

# chematic
sim = chematic.tanimoto(mol1.ecfp4(), mol2.ecfp4())
```

## Bulk processing

```python
# RDKit (single-threaded)
mols = [Chem.MolFromSmiles(s) for s in smiles_list]
fps  = [AllChem.GetMorganFingerprintAsBitVect(m, 2, 2048) for m in mols]

# chematic (automatic multi-core parallelism via Rayon)
fps = chematic.bulk.ecfp4(smiles_list)  # numpy array (N, 2048)
```

## Similarity matrix

```python
# RDKit
from rdkit.DataStructs import BulkTanimotoSimilarity
matrix = [[TanimotoSimilarity(fp, fp2) for fp2 in fps] for fp in fps]

# chematic
matrix = chematic.bulk.tanimoto(smiles_a, smiles_b)  # numpy (M, N) float32
```

## SMARTS substructure search

```python
# RDKit
query = Chem.MolFromSmarts("[OH]")
matches = mol.GetSubstructMatches(query)

# chematic — bool
chematic.smarts_match("[OH]", mol)

# chematic — atom indices
chematic.smarts_find("[OH]", mol)   # [[3], [7], ...]
```

## Reaction SMARTS

```python
# RDKit
from rdkit.Chem import rdChemReactions
rxn = rdChemReactions.ReactionFromSmarts("[OH]>>[O-]")
matches = rxn.RunReactants((mol,))

# chematic
chematic.reaction_smarts_match("[OH]>>[O-]", "CCO>>CC[O-]")  # bool
```

## Maximum Common Substructure (MCS)

```python
# RDKit
from rdkit.Chem import rdFMCS
result = rdFMCS.FindMCS([mol1, mol2])
mcs_mol = Chem.MolFromSmarts(result.smartsString)

# chematic
mcs = chematic.find_mcs([mol1, mol2])
print(mcs.smiles)
```

## SMIRKS reactions

```python
# RDKit
from rdkit.Chem import AllChem
rxn = AllChem.ReactionFromSmarts("[OH:1]>>[O-:1]")
products = rxn.RunReactants((mol,))

# chematic
products = chematic.run_smirks("[OH:1]>>[O-:1]", [mol])
```

## Standardisation

```python
# RDKit
from rdkit.Chem.MolStandardize import rdMolStandardize
clean = rdMolStandardize.Cleanup(mol)
frags = rdMolStandardize.LargestFragmentChooser().choose(clean)
uncharge = rdMolStandardize.Uncharger().uncharge(frags)

# chematic (one-liner)
clean = mol.standardize()

# individual steps
mol.largest_fragment()
mol.neutralize()
mol.remove_stereo()
mol.remove_isotopes()
```

## Tautomers

```python
# RDKit
from rdkit.Chem.MolStandardize import rdMolStandardize
enumerator = rdMolStandardize.TautomerEnumerator()
canonical = enumerator.Canonicalize(mol)
all_tautomers = enumerator.Enumerate(mol)

# chematic
canonical = mol.canonical_tautomer()
all_tautomers = mol.enumerate_tautomers()
```

## Murcko scaffold

```python
# RDKit
from rdkit.Chem.Scaffolds import MurckoScaffold
scaffold = MurckoScaffold.GetScaffoldForMol(mol)
generic = MurckoScaffold.MakeScaffoldGeneric(scaffold)

# chematic
scaffold = mol.scaffold()
generic  = mol.generic_scaffold()
```

## SDF file reading

```python
# RDKit
from rdkit.Chem import SDMolSupplier
for mol in SDMolSupplier("library.sdf"):
    if mol: print(Descriptors.MolWt(mol))

# chematic
for record in chematic.iter_sdf("library.sdf"):
    print(record.mol.mw)
```

## 2D depiction

```python
# RDKit (Jupyter, returns PIL image)
from rdkit.Chem import Draw
Draw.MolToImage(mol)

# chematic (Jupyter, SVG)
from IPython.display import SVG
SVG(mol.svg())

# with atom highlighting
SVG(mol.svg_highlighted([0, 1, 2], color="#FF6B6B"))

# grid
SVG(chematic.depict_grid([mol1, mol2, mol3], cols=3))
```

## SASA (3D solvent-accessible surface area)

```python
# RDKit
from rdkit.Chem import rdFreeSASA
rdFreeSASA.CalcSASA(mol_with_3d_coords)

# chematic (generates 3D coords internally)
mol.sasa()           # total SASA in A^2
mol.sasa_per_atom()  # per-atom list
```

## pKa prediction

```python
# RDKit — not in the standard library (external tools required)

# chematic — built-in
pka = mol.pka()
print(pka["most_acidic"])   # 3.49
print(pka["most_basic"])    # None
```

## ADMET profile

```python
# RDKit — not in the standard library (pkasolver, SwissADME, etc. required)

# chematic — built-in
profile = mol.admet()
# {"bbb": False, "bbb_score": -1.2, "caco2": -5.1, "herg_risk": 0.3, "cyp3a4_risk": 0.4}
```

## All descriptors to a DataFrame

```python
# RDKit
from rdkit.Chem import Descriptors
import pandas as pd
desc_names = [x[0] for x in Descriptors._descList]
data = [{name: Descriptors.__dict__[name](mol) for name in desc_names} for mol in mols]
df = pd.DataFrame(data)

# chematic
df = pd.DataFrame(chematic.bulk.descriptors(smiles_list))
```

---

## Features only in chematic

- **pKa prediction** — built-in, no external tools
- **ADMET profile** — BBB, Caco-2, hERG, CYP3A4 in a single call
- **WASM support** — runs in the browser (~1.1 MB gzip bundle)
- **MCP server** — direct integration with AI agents
- **Pure Rust** — no conda, works in Docker / serverless / CI without extra setup
- **Atropisomer detection** — `mol.atropisomers()` detects biaryl and allene axes
