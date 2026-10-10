"""Freeze RDKit 2026.03.1 query atom sets for bounded SMARTS primitives."""
import json
from pathlib import Path
from rdkit import Chem,rdBase
assert rdBase.rdkitVersion=='2026.03.1',rdBase.rdkitVersion
queries=['[C;H3]','[C;v4]','[^1]','[^2]','[^3]','[D0]','[D1]','[X4]','[H0]','[h1]','[13C]','[!C]','[C,N]','[C;!H0]','[R]','[R0]','[R1]','[R2]','[r3]','[r4]','[r5]','[r6]','*@*','*!@*','*~*','[N+]','[O-]','[$(C=O)]','[$([C;$([C]=O)])]','[!$(C=O)]','[$(N-C=O),$(O-C=O)]','[$(C1CC1)]','[$(c1ccccc1)]','[$([NX3;H2,H1;!$(NC=O)])]','[C@H]','[C@@H]','C.C','[#6]-;!@[#7]']
sources=['','C','CCO','CNC=O','CC(=O)O','[NH4+]','C[N+](C)(C)C','CC(=O)[O-]','[13CH3]C','[2H]O','C#N','C=C','N=N','c1ccccc1','c1ccncc1','c1cc[nH]c1','C1CC1','C1CCC1','C1CCCC1','C1CCCCC1','C1CC2CCC1C2','F[C@H](Cl)Br','F[C@@H](Cl)Br','NC(=O)N','[He]','CC.CN']
rows=[]
queries += ['[z1]', '[Z1]', '[D2]', '[X2]', '[x0]', '[x2]', '[v0]', '[^0]', '[^4]', '[^5]',
            'F[C@H](Cl)Br', 'F[C@@H](Cl)Br', '[$(F[C@H](Cl)Br)]', 'F/C=C/Cl', 'F/C=C\\Cl']
sources += ['P(F)(F)(F)(F)F','S(F)(F)(F)(F)(F)F','[O]','[C]','[H][H]','[Na+]','[NH2]',
            'F/C=C/Cl','F/C=C\\Cl']
for source in sources:
    mol=Chem.MolFromSmiles(source);assert mol is not None
    for query in queries:
        q=Chem.MolFromSmarts(query);assert q is not None,query
        # The default matcher contract is chirality-agnostic. Full relative
        # stereochemistry matching is a separate compatibility workstream.
        for use_chirality in [False]:
            matches=sorted({tuple(sorted(m)) for m in mol.GetSubstructMatches(q,maxMatches=10000,useChirality=use_chirality)})
            rows.append(dict(smiles=source,smarts=query,use_chirality=use_chirality,sets=matches))
root=Path(__file__).resolve().parents[1]
(root/'validation/rdkit-2026.03.1-smarts-boundary.json').write_text('{\n  "rdkit_version": "'+rdBase.rdkitVersion+'",\n  "rows": [\n'+',\n'.join('    '+json.dumps(row) for row in rows)+'\n  ]\n}\n')
print('wrote',len(rows),'cases')
