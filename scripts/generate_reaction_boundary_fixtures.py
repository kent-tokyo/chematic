"""Record tetrahedral reaction-template permutations with RDKit 2026.03.6."""
import itertools,json
from pathlib import Path
from rdkit import Chem,rdBase
from rdkit.Chem import AllChem
assert rdBase.rdkitVersion=='2026.03.6',rdBase.rdkitVersion
rows=[]
def record(template,source):
    rxn=AllChem.ReactionFromSmarts(template)
    products=[]
    for result in rxn.RunReactants((Chem.MolFromSmiles(source),)):
        values=[]
        for mol in result:
            Chem.SanitizeMol(mol)
            values.append(Chem.MolToSmiles(mol))
        products.append(values)
    rows.append(dict(template=template,reactant=source,products=products))
for atoms,source in [(['[F:2]','[Cl:3]','[Br:4]','[I:5]'],'F[C@](Cl)(Br)I'),
                     (['[F:2]','[Cl:3]','[Br:4]'],'F[C@H](Cl)Br')]:
    hydrogen='H' if len(atoms)==3 else ''
    for react_tag,prod_tag in itertools.product(['','@','@@'],repeat=2):
        react='[C'+react_tag+hydrogen+':1]'+''.join('('+a+')' for a in atoms[:-1])+atoms[-1]
        for perm in itertools.permutations(atoms):
            product='[C'+prod_tag+hydrogen+':1]'+''.join('('+a+')' for a in perm[:-1])+perm[-1]
            record(react+'>>'+product,source)
for source in ['F[C@](Cl)(Br)I','F[C@@](Cl)(Br)I','F[C@H](Cl)Br','F[C@@H](Cl)Br','[2H][C@](F)(Cl)Br','F[C@H]1CCC(C)C1']:
    for template in ['[F:1]>>[F:1]','[F:1][C:2]>>[F:1][C:2]','[F:1][C:2]>>[C:2][F:1]']:
        record(template,source)
root=Path(__file__).resolve().parents[1]
(root/'validation/rdkit-2026.03.6-reaction-boundary.json').write_text('{\n  "rdkit_version": "'+rdBase.rdkitVersion+'",\n  "rows": [\n'+',\n'.join('    '+json.dumps(row) for row in rows)+'\n  ]\n}\n')
print(f'wrote {len(rows)} cases')
