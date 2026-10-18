"""Independent RDKit SD generation and Maestro cross-read; no RDKit at runtime."""
import hashlib, io, json, math
from pathlib import Path
import chematic as ch
import rdkit
from rdkit import Chem

ROOT = Path(__file__).resolve().parent
SMILES = ['O','N','C','CC','CCC','CCCC','CO','CCO','COC','CCOC','C(CO)O',
 'C=O','CC=O','CC(C)=O','O=C=O','CC(=O)O','CC(=O)[O-]','C#N','CC#N',
 'C=C','CC=C','C#C','CC#C','C1CC1','C1CCC1','C1CCCCC1','CC(C)O',
 '[NH4+]','C[NH3+]','[H][H]','[H+]','[Na+].[Cl-]','[Li+].[F-]',
 'P','S','CF','CCl','CBr','CI']
report = dict(rdkit=rdkit.__version__, molecules=len(SMILES), v2000_cases=0,
              graph_cases=0, checked_atoms=0, checked_bonds=0, max_coordinate_error=0.)


def read(text):
    mol = next(Chem.MaeMolSupplier(io.BytesIO(text.encode('ascii')), sanitize=True, removeHs=False))
    assert mol is not None
    return mol


def verify(expected, actual, title, xyz):
    assert expected.GetNumAtoms() == actual.GetNumAtoms()
    assert expected.GetNumBonds() == actual.GetNumBonds()
    assert actual.GetProp('_Name') == title
    for a, b in zip(expected.GetAtoms(), actual.GetAtoms()):
        assert (a.GetAtomicNum(),a.GetFormalCharge(),a.GetIsotope()) == (b.GetAtomicNum(),b.GetFormalCharge(),b.GetIsotope())
        assert a.GetNumImplicitHs() == b.GetNumImplicitHs() == 0
    bonds = lambda mol: sorted((min(b.GetBeginAtomIdx(),b.GetEndAtomIdx()),max(b.GetBeginAtomIdx(),b.GetEndAtomIdx()),str(b.GetBondType())) for b in mol.GetBonds())
    assert bonds(expected) == bonds(actual)
    for i, point in enumerate(xyz):
        got = tuple(actual.GetConformer().GetAtomPosition(i))
        error = max(abs(a-b) for a,b in zip(point,got))
        assert error <= max(1e-12, max(abs(v) for v in point)*1e-15)
        report['max_coordinate_error'] = max(report['max_coordinate_error'],error)
    report['checked_atoms'] += expected.GetNumAtoms()
    report['checked_bonds'] += expected.GetNumBonds()


for index, smiles in enumerate(SMILES):
    base = Chem.AddHs(Chem.MolFromSmiles(smiles))
    # Nonaromatic, stereo-free authored graphs, ordinary closed-shell valences.
    for scale in (.0001, .125, 1., 25.):
        for offset in (-7., 0., 37.):
            for title in ('', 'ASCII title', 'braces {} ::: # "quote" \\ slash'):
                mol = Chem.Mol(base)
                mol.SetProp('_Name', title)
                conf = Chem.Conformer(mol.GetNumAtoms()); conf.Set3D(True)
                for i in range(mol.GetNumAtoms()):
                    conf.SetAtomPosition(i, (round(offset+scale*math.sin(i+index),4),
                                             round(scale*math.cos(i*.7),4),
                                             round(scale*(i%5-.3),4)))
                mol.AddConformer(conf,assignId=True)
                block = Chem.MolToMolBlock(mol,forceV3000=False)
                sdf = block + '>  <Cat. No.>\nfixture-001\n\n$$$$\n'
                expected = Chem.MolFromMolBlock(block,sanitize=True,removeHs=False)
                assert expected is not None
                xyz = [tuple(expected.GetConformer().GetAtomPosition(i)) for i in range(expected.GetNumAtoms())]
                result = ch.v2000_to_mae_checked(sdf,allow_metadata_loss=True)
                assert result.source_text == sdf
                verify(expected,read(result.text),title,xyz)
                report['v2000_cases'] += 1
                native = ch.from_mol_block(block)
                graph = ch.write_mae_checked(native,xyz,title=title,allow_metadata_loss=True)
                assert result.text == graph.text
                verify(expected,read(graph.text),title,xyz)
                report['graph_cases'] += 1
                if index == 0 and scale == 1. and offset == 0. and title == 'ASCII title':
                    (ROOT/'water.sdf').write_text(sdf)
                    (ROOT/'water.mae').write_text(result.text)

# Direct geometry edge controls, independently parsed by the consumer.
for value in (0.,1e-99,-1e-99,1e8,-1e8):
    mol=ch.from_smiles('[H][H]')
    xyz=[[value,0.,0.],[0.,value,0.]]
    out=ch.write_mae_checked(mol,xyz,allow_metadata_loss=True)
    ref=Chem.MolFromSmiles('[H][H]')
    verify(ref,read(out.text),'',xyz)
report['coordinate_edge_cases']=5
for smiles in ('FP(F)(F)(F)F', 'FS(F)(F)(F)(F)F'):
    reference=Chem.MolFromSmiles(smiles)
    assert reference is not None
    native=ch.from_smiles(smiles)
    xyz=[[i*.1,0.,.25] for i in range(native.atom_count)]
    out=ch.write_mae_checked(native,xyz,allow_metadata_loss=True)
    verify(reference,read(out.text),'',xyz)
report['hypervalence_graph_controls']=2
report['reader_binary_sha256']=hashlib.sha256(Path(Chem.rdmolfiles.__file__).read_bytes()).hexdigest()
report['oracle_sha256']=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
(ROOT/'oracle.json').write_text(json.dumps(report,indent=2)+'\n')
print(json.dumps(report,sort_keys=True))
