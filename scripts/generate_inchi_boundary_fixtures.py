"""Record native InChI reconstruction references with RDKit 2026.03.1."""
import json
from pathlib import Path
from rdkit import Chem, rdBase
assert rdBase.rdkitVersion == "2026.03.1",rdBase.rdkitVersion
smiles=[
    "C[N+](=O)[O-]","O=[N+]([O-])c1ccccc1","[N-]=[N+]=N","CN=[N+]=[N-]",
    "N#N","CN=O","CN(O)O","C[N+](C)([O-])C","C[NH2+]C","C[NH3+]",
    "CS(=O)C","CS(=O)(=O)C","CS(=O)(=O)O","CS(=O)[O-]","O=S=O","O=S(O)O",
    "OS(=O)(=O)O","O=S(=O)([O-])[O-]","CSSC","C=S=O","N=S(=O)(C)C",
    "O=Cl(=O)(=O)[O-]","O=Cl(=O)[O-]","O=Cl[O-]","[Cl+3]([O-])([O-])[O-]",
    "O=Br(=O)[O-]","C[Se](=O)C","C[Se](=O)(=O)C","P(=O)(O)(O)O","CP(=O)(O)O",
    "F[B-](F)(F)F","[Na+].[O-]C=O","[Li+].[Cl-]","[CH3]","[O][O]","[NH2]",
    "[13CH3][15NH2]","[2H]O[2H]","[18OH2]","[2H][C@](F)(Cl)Br","[13CH4]",
    "F/C=C/Cl","F/C=C\\Cl","C/C=C/C=C/C","F[C@H](Cl)Br","F[C@@H](Cl)Br",
    "N[C@@H](C)C(=O)O","N[C@H](C)C(=O)O","C1CC2CCC1C2","c1ccc2[nH]ccc2c1",
    "O=c1cccc[nH]1","c1cc[n+]([O-])cc1","[nH]1nnnc1","N#CC#N","N=C(N)N",
]
rows=[]
for text in smiles:
    mol=Chem.MolFromSmiles(text)
    assert mol is not None,text
    inchi=Chem.MolToInchi(mol)
    assert inchi,text
    rebuilt=Chem.MolFromInchi(inchi)
    rows.append(dict(source_smiles=text,inchi=inchi,smiles=Chem.MolToSmiles(rebuilt) if rebuilt is not None else None))
path=Path(__file__).resolve().parents[1]/"validation/rdkit-2026.03.1-inchi-boundary.json"
path.write_text(json.dumps(dict(rdkit_version=rdBase.rdkitVersion,rows=rows),indent=2)+"\n")
print(f"wrote {len(rows)} cases")
