"""Freeze native InChI reading of unusual-valence source encodings.

Sources exercise the cleanup families documented in RDKit Release_2026_03_1
External/INCHI-API/inchi.cpp. The unsanitized sources are generator inputs;
Rust receives only the resulting identifiers and compares native read/reject.
"""
import json
from pathlib import Path
from rdkit import Chem,rdBase
assert rdBase.rdkitVersion=='2026.03.1'
sources=['C1=NN=[N-]=N1','CN1=CCN=CC=1','CN1=NCOCC=1',
         'CN1=NC(=O)OCC=1','N=C1N=CN=N1','N=C1C=CN=N1',
         'CN(=CC=NNC)=C','CN(=C)=C','CN(=C)=[NH2+]',
         'CN(=C)#C[NH-]','CN(=[Si-](C)(C)C)=[Si-](C)(C)C',
         'CN(=C)=O','CN(=C)=S','CN(=C)=F','CN(=C)=Cl',
         'C[S-](=O)(=O)=O','[S-](=N)(=O)(=O)C',
         '[S-](=O)(=O)=CC=N','[S-](=N)(=C)#N',
         '[Cl-](=O)(=O)(=O)=O','[Cl+](=O)(=O)[O-]',
         'Cl#S','Br#[Se]','C[N-](=N)=C','CN(=C)=[NH+]C']
rows=[]
for source in sources:
 raw=Chem.MolFromSmiles(source,sanitize=False);assert raw is not None,source
 raw.UpdatePropertyCache(strict=False)
 inchi=Chem.MolToInchi(raw)
 assert inchi,source
 result=Chem.MolFromInchi(inchi)
 expected=Chem.MolToSmiles(result) if result is not None else None
 rows.append(dict(raw_source=source,inchi=inchi,canonical=expected))
 print(source,expected)
(Path(__file__).resolve().parents[1]/'validation/rdkit-2026.03.1-inchi-unusual-boundary.json').write_text(json.dumps(dict(rdkit_version=rdBase.rdkitVersion,rows=rows),indent=2)+'\n')
print(len(rows),'InChI input references')
