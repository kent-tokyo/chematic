# chematic cookbook

Copy-paste examples for common Python workflows. Install with `pip install chematic`, then:

```python
import chematic

aspirin = chematic.from_smiles("CC(=O)Oc1ccccc1C(=O)O")
```

## Inspect a molecule

```python
print(aspirin.smiles)
print(aspirin.formula, aspirin.mw)
print(aspirin.logp, aspirin.tpsa)
print(aspirin.hbd, aspirin.hba, aspirin.qed)

print(aspirin.lipinski_passes)
print(aspirin.veber_passes)
print(aspirin.pains_passes)
print(aspirin.brenk_passes)
```

`pka()` and `admet()` are rule-based screening helpers, not experimental or clinical predictions:

```python
print(aspirin.pka())
print(aspirin.admet())
print(aspirin.describe())
```

## Batch descriptors and fingerprints

```python
import pandas as pd

smiles = ["CCO", "c1ccccc1", "CC(=O)O", "c1cccnc1"]
df = pd.DataFrame(chematic.bulk.descriptors(smiles))
print(df[["mw", "logp", "tpsa", "qed"]])

# Dense uint8 matrix, shape (N, 2048)
X = chematic.bulk.ecfp4(smiles)

# Use named compatibility APIs when RDKit bit identity is required.
morgan = aspirin.rdkit_ecfp4_numpy()
avalon = aspirin.rdkit_avalon_fp_numpy()
```

## Similarity search

```python
ibuprofen = chematic.from_smiles("CC(C)Cc1ccc(CC(C)C(=O)O)cc1")
score = chematic.tanimoto(aspirin.ecfp4(), ibuprofen.ecfp4())

matrix = chematic.bulk.tanimoto(smiles, smiles)
scores = chematic.bulk.tanimoto_search(aspirin.smiles, smiles)

index = chematic.SimilarityIndex.from_smiles(smiles)
hits = index.search(aspirin.smiles, threshold=0.3, k=10)
```

## SMARTS and structure operations

```python
pattern = "[CX3](=O)[OX2H1]"
print(chematic.smarts_match(pattern, aspirin))
print(chematic.smarts_find(pattern, aspirin))

salt = chematic.from_smiles("[Na+].[O-]c1ccccc1")
print(salt.standardize().smiles)
print(aspirin.scaffold().smiles)
print(aspirin.generic_scaffold().smiles)

tautomers = chematic.from_smiles("OC1=CC=CC=N1").enumerate_tautomers()
canonical = chematic.from_smiles("OC1=CC=CC=N1").canonical_tautomer()

for fragment in chematic.from_smiles("CC(=O)Nc1ccc(O)cc1").brics_fragments():
    print(fragment.smiles)
```

`standardize()` applies the documented chematic pipeline. If RDKit replacement compatibility matters, review the [compatibility scope](compatibility-scope.md) instead of assuming every cleanup policy is identical.

## MCS, R-group decomposition, and reactions

```python
paracetamol = chematic.from_smiles("CC(=O)Nc1ccc(O)cc1")
mcs = chematic.find_mcs([aspirin, paracetamol])
print(mcs.smiles if mcs else None)

mols = [
    chematic.from_smiles("Cc1ccc(CC)cc1"),
    chematic.from_smiles("CCc1ccc(N)cc1"),
]
columns = chematic.rgroup_decompose_columns(
    "c1cc([*:1])ccc1[*:2]", mols
)
print(columns["R1"], columns["R2"])

phenol = chematic.from_smiles("c1ccccc1O")
products = chematic.run_smirks("[OH:1]>>[O-:1]", [phenol])
print([[mol.smiles for mol in product_set] for product_set in products])
```

R-group decomposition intentionally supports a bounded contract. Ambiguous or unsupported labels fail explicitly rather than being guessed.

## Stream SDF and build a report

```python
records = chematic.iter_sdf("library.sdf")
rows = [
    {"name": r.name, "smiles": r.mol.smiles, "mw": r.mol.mw}
    for r in records
]
df = pd.DataFrame(rows)

report = chematic.report(
    [aspirin, ibuprofen],
    names=["aspirin", "ibuprofen"],
    title="Compound review",
)
report.save("report.html")
```

`iter_sdf()` avoids loading the entire file at once. Checked/source-preserving I/O APIs are preferable when loss reporting matters; see [errors and limits](error-and-limits.md).

## Draw in Jupyter

```python
from IPython.display import SVG, display

display(SVG(aspirin.svg()))

matches = chematic.smarts_find("[CX3](=O)[OX2H1]", aspirin)
atoms = [atom for match in matches for atom in match]
display(SVG(aspirin.svg_highlighted(atoms, color="#FF6B6B")))

display(SVG(chematic.depict_grid([aspirin, ibuprofen], cols=2)))
```

## Convert InChI

```python
print(aspirin.inchi)
print(aspirin.inchikey)

roundtrip = chematic.from_inchi(aspirin.inchi)
assert roundtrip.smiles == aspirin.smiles
```

For the complete API, see the [Python reference](language-bindings.md). For RDKit migration, use the [RDKit cheat sheet](rdkit_cheatsheet.md).
