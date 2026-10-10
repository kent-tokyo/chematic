"""Record bounded SMARTS serialization and 2D references with RDKit 2026.03.1."""
import itertools
import json
from pathlib import Path
from rdkit import Chem, rdBase
from rdkit.Chem import rdDepictor
assert rdBase.rdkitVersion == '2026.03.1', rdBase.rdkitVersion
root=Path(__file__).resolve().parents[1]
primitives=['*','#6','#7','n','H0','D2','X3','v4','h1','x2','z1','Z1','r5','R2','^2','+2','13C','!C','$([C]=O)','!$(*#*)','D{2-4}','r{5-}','R{-2}']
queries=['['+a+op+b+']' for a,b in itertools.product(primitives,repeat=2) for op in [';',',','&']]
bonds=['-','=','#',':','~','@','!@','!~','!-','!=','!#','!:']
queries += ['*'+a+op+b+'*' for a,b in itertools.product(bonds,repeat=2) for op in [';',',','&']]
queries+=['[C@H](F)(Cl)Br','[C@@]1(F)CCCCC1','[nH+]1ccccc1','[#6:12]~[#7:3]','C1~C~C~C~C~C1','[$(C=O),$(N#C);!H0]']
query_rows=[]
for query in queries:
    mol=Chem.MolFromSmarts(query)
    assert mol is not None,query
    query_rows.append(dict(query=query,canonical=Chem.MolToSmarts(mol)))
sources=['Cl[Pt@SP1](Cl)(N)N','N[Co@OH1](N)(N)(N)(Cl)Cl','F[P@TB1](F)(F)(Cl)Br','[H]/C(=C(/[H])C)C','C1C2CC3CC1CC(C2)C3','C1CC2CCC3CCCC4CCCC1C2C34','CC.CCC.O','[Na+].[Cl-]']
for i,line in enumerate((root/'validation/results/descriptor_census_unbound.jsonl').read_text().splitlines()):
    if i%17==0 and i<10000:
        row=json.loads(line)
        if row.get('parse_ok') and row.get('smiles'):
            sources.append(row['smiles'])
geometry_rows=[]
for source in dict.fromkeys(sources):
    mol=Chem.MolFromSmiles(source)
    assert mol is not None,source
    rdDepictor.Compute2DCoords(mol,forceRDKit=True)
    xyz=mol.GetConformer()
    geometry_rows.append(dict(smiles=source,coords=[[xyz.GetAtomPosition(i).x,xyz.GetAtomPosition(i).y] for i in range(mol.GetNumAtoms())]))
path=root/'validation/rdkit-2026.03.1-serialization-geometry-boundary.json'
path.write_text(json.dumps(dict(rdkit_version=rdBase.rdkitVersion,queries=query_rows,geometries=geometry_rows),indent=2)+'\n')
print(f'wrote {len(query_rows)} queries and {len(geometry_rows)} geometries')
