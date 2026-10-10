"""Generate descriptor boundary references; requires RDKit 2026.03.1."""
import json
from pathlib import Path
from rdkit import Chem, rdBase
from rdkit.Chem import rdMolDescriptors as d, GraphDescriptors as g
assert rdBase.rdkitVersion == "2026.03.1", rdBase.rdkitVersion
smiles = [
    "", "[H][H]", "[He]", "B", "B(F)(F)F", "BF", "C#N", "CC#CC", "C=C", "CC(=C)C",
    "CC(C)C", "CC(C)(C)C(C)(C)C", "OC(F)(Cl)Br", "CCI", "CP(=O)(O)O", "P", "[PH4+]",
    "CS", "CS(=O)C", "CS(=O)(=O)C", "CS(=O)(=O)[O-]", "[S-]", "O=S=O", "CS(=S)S",
    "N", "N=N", "N#N", "CN=O", "[NH4+]", "[N-]=[N+]=N", "C[N+](C)(C)C",
    "[OH-]", "[OH3+]", "[O-]C=O", "[Cl-]", "[Br-]", "[I-]", "[Na+]", "[Mg+2]",
    "[SiH4]", "[AsH3]", "[SeH2]", "[Zn+2]", "O.O", "CC.[Na+].[Cl-]",
    "c1ccoc1", "c1ccsc1", "c1ncc[nH]1", "c1cc[nH+]cc1", "c1cc[n+]([O-])cc1",
    "C1=CC=CC=C1", "c1ccc2ccccc2c1", "C1CC2CCC1C2", "C12C3C4C1C5C2C3C45",
    "C1CC2CC3CC1CC(C2)C3", "C1CC2(C1)CCC1(CC2)CCC1", "C1CCC2(CC1)CCCC2",
    "C1=CC2CCC1C2", "C1CCC2C3CCC4CCCCC4C3CCC12", "F[C@H](Cl)Br", "CC(O)C(Cl)Br",
]
rows = []
for text in smiles:
    mol = Chem.MolFromSmiles(text)
    assert mol is not None, text
    for explicit_h in [False, True]:
        source = Chem.AddHs(mol) if explicit_h else mol
        encoded = Chem.MolToSmiles(source, allHsExplicit=explicit_h)
        # These descriptors use the hydrogen-suppressed graph.
        heavy = Chem.RemoveHs(source)
        Chem.GetSymmSSSR(heavy)
        c = [0]*8
        keys = [(Chem.HybridizationType.SP,1),(Chem.HybridizationType.SP,2),
                (Chem.HybridizationType.SP2,1),(Chem.HybridizationType.SP2,2),(Chem.HybridizationType.SP2,3),
                (Chem.HybridizationType.SP3,1),(Chem.HybridizationType.SP3,2),(Chem.HybridizationType.SP3,3)]
        for atom in heavy.GetAtoms():
            key = (atom.GetHybridization(), sum(a.GetAtomicNum()!=1 for a in atom.GetNeighbors()))
            if atom.GetAtomicNum()==6 and key in keys: c[keys.index(key)]+=1
        rows.append(dict(smiles=encoded,mqn=list(d.MQNs_(heavy)),alpha=d.CalcHallKierAlpha(heavy),
                         crippen_no_hs=list(d._CalcCrippenContribs(Chem.MolFromSmiles(encoded))),
                         chi_v=[getattr(d,f"CalcChi{i}v")(heavy) for i in range(5)],
                         rings=d.CalcNumRings(heavy),carbon_types=c,
                         counts=[sum(a.GetAtomicNum()==z for a in heavy.GetAtoms()) for z in [6,7,8,9,17,35,53,16,15]]))
path=Path(__file__).resolve().parents[1]/"validation/rdkit-2026.03.1-descriptor-boundary.json"
path.write_text('{\n  "rdkit_version": "' + rdBase.rdkitVersion + '",\n  "rows": [\n' +
                ',\n'.join('    ' + json.dumps(row) for row in rows) + '\n  ]\n}\n')
print(f"wrote {len(rows)} cases")
