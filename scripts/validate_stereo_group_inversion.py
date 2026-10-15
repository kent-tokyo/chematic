#!/usr/bin/env python3
"""Independent drawn-configuration oracle for relative stereo-group alternatives."""
import argparse
import hashlib
import json
from pathlib import Path
import random
import sys
import chematic
from rdkit import Chem, RDLogger, rdBase

ROOT=Path(__file__).resolve().parents[1]
FIXTURE=ROOT/'validation/rdkit_issues/stereo/group-inversion.json'


def digest(data):
    return hashlib.sha256(data).hexdigest()


def parse(smiles):
    params=Chem.SmilesParserParams();params.removeHs=False
    mol=Chem.MolFromSmiles(smiles,params)
    assert mol is not None
    return mol


def inchi(smiles):
    return Chem.MolToInchi(parse(smiles))


def reference(case, native):
    mol=parse(case['smiles'])
    for index in case['atoms']:mol.GetAtomWithIdx(index).InvertChirality()
    for index in case['bonds']:
        a,b,*_=native.bond_table[index]
        bond=mol.GetBondBetweenAtoms(a,b)
        assert bond is not None
        stereo=bond.GetStereo()
        assert stereo in (Chem.BondStereo.STEREOE,Chem.BondStereo.STEREOZ)
        bond.SetStereo(Chem.BondStereo.STEREOZ if stereo==Chem.BondStereo.STEREOE else Chem.BondStereo.STEREOE)
    for bond in mol.GetBonds():bond.SetBondDir(Chem.BondDir.NONE)
    Chem.SetDoubleBondNeighborDirections(mol)
    smiles=Chem.MolToSmiles(mol)
    expected=Chem.MolToInchi(mol)
    assert inchi(smiles)==expected, 'RDKit reference SMILES lost its declared configuration'
    return smiles,expected


def orders(n):
    rng=random.Random(9705)
    result=[list(range(n)),list(reversed(range(n)))]
    while len(result)<24:
        order=list(range(n));rng.shuffle(order)
        if order not in result:result.append(order)
    return result


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,default=ROOT/'validation/results/rdkit-stereo-group-inversion-2026-10-15.json')
    args=parser.parse_args()
    if rdBase.rdkitVersion!='2026.09.1':raise RuntimeError(f'expected RDKit 2026.09.1, got {rdBase.rdkitVersion}')
    Chem.SetUseLegacyStereoPerception(True);RDLogger.DisableLog('rdApp.warning')
    fixture=json.loads(FIXTURE.read_text());checked=[];references=[]
    for case in fixture['cases']:
        base=chematic.from_smiles(case['smiles']);external=parse(case['smiles'])
        assert base.atom_count==external.GetNumAtoms()
        assert [(a[1],a[2],a[3]) for a in base.atom_table]==[(a.GetAtomicNum(),a.GetFormalCharge(),a.GetIsAromatic()) for a in external.GetAtoms()]
        expected_bonds=[(b.GetBeginAtomIdx(),b.GetEndAtomIdx(),str(b.GetBondType()),b.GetIsAromatic()) for b in external.GetBonds()]
        normalize=lambda bonds:sorted((min(a,b),max(a,b),kind,aromatic) for a,b,kind,aromatic in bonds)
        assert normalize(base.bond_table)==normalize(expected_bonds), 'oracle/native indexed graph mismatch'
        original=inchi(case['smiles']);expected_smiles,expected=reference(case,base)
        assert inchi(base.smiles)==original
        references.append({'id':case['id'],'source_inchi':original,'alternative_smiles':expected_smiles,'alternative_inchi':expected,'same_drawn_configuration':expected==original})
        for kind in ('and','or'):
            source=base.with_stereo_groups([(kind,1,case['atoms'],case['bonds'])])
            for order in orders(base.atom_count):
                inverse={old:new for new,old in enumerate(order)}
                numbered=source.renumber_atoms(order)
                before=numbered.smiles
                out=numbered.invert_stereo_group(0)
                assert out.stereo_groups==[(kind,1,[inverse[a] for a in case['atoms']],case['bonds'])]
                assert out.bond_table==numbered.bond_table
                assert out.atom_table==numbered.atom_table
                assert out.atom_count==numbered.atom_count
                assert inchi(out.smiles)==expected
                assert inchi(out.invert_stereo_group(0).smiles)==original
                assert numbered.smiles==before
                checked.append({'id':case['id'],'kind':kind,'order':order,'alternative_inchi':expected,'output_smiles':out.smiles})
    # Reproduce the upstream mixed canonicalization defect separately; it is
    # deliberately not the alternative-configuration oracle used above.
    mol=Chem.RWMol(parse('C/C=C/[C@H](O)CC'))
    mol.SetStereoGroups([Chem.CreateStereoGroup(Chem.StereoGroupType.STEREO_AND,mol,[3],[1])])
    canonical=Chem.CanonicalizeStereoGroups(mol.GetMol())
    after=[{'atoms':[a.GetIdx() for a in g.GetAtoms()],'bonds':[b.GetIdx() for b in g.GetBonds()]} for g in canonical.GetStereoGroups()]
    assert sum(len(g['bonds']) for g in after)==0, 'upstream behavior changed; inspect before rebaseline'
    paths=['crates/chematic-chem/src/stereo_group_inversion.rs','crates/chematic-py/src/mol_methods/stereo_groups.rs','crates/chematic-chem/tests/stereo_group_inversion.rs','crates/chematic-py/tests/test_stereo_group_inversion.py','crates/chematic-py/python/chematic/__init__.pyi','scripts/validate_stereo_group_inversion.py']
    record={
        'schema_version':1,'status':'local-source-candidate','issue':'https://github.com/rdkit/rdkit/issues/9705',
        'scope':'one checked full alternative for an AND/OR group containing tetrahedral atoms/geometric double bonds; no group canonicalization or axial inversion',
        'base':'cac480fb source snapshot plus preserved earlier RDKit/Open Babel candidates; requires checked stereo-group construction candidate',
        'runtime':{'python':sys.version.split()[0],'rdkit':rdBase.rdkitVersion,'settings':'legacy stereo; retain explicit hydrogen atoms; regenerate double-bond neighbour directions after direct Stereo E/Z swaps','chematic_extension_sha256':digest(Path(chematic.__file__).resolve().read_bytes())},
        'source_sha256':{s:digest((ROOT/s).read_bytes()) for s in paths},
        'fixture':{'path':str(FIXTURE.relative_to(ROOT)),'sha256':digest(FIXTURE.read_bytes()),'encoding':'UTF-8 LF',**fixture},
        'results':{'comparisons':len(checked),'membership_mismatches':0,'drawn_configuration_mismatches':0,'indexed_graph_mismatches':0,'double_inversion_mismatches':0,'numberings_per_kind_fixture':24,'kinds':['and','or'],'fixtures':len(fixture['cases']),'cases_sha256':digest(json.dumps(checked,sort_keys=True,separators=(',',':')).encode()),'references':references},
        'rdkit_canonicalization_diagnostic':{'input_atoms':[3],'input_bonds':[1],'output_groups':after,'note':'Group bond loss reproduced; direct atom/bond inversion above is the independent drawn-configuration reference'},
        'limits':'ABS, unsupported/non-tetrahedral/unspecified frames, axial and non-double members, small-ring double bonds, conflicting/unrepresentable direction representations fail atomically. No 3D feasibility or meso reduction; plain SMILES/InChI validate the drawn representative only, not enhanced-group identity. All native atom CIP annotations are cleared after a successful alternative.',
    }
    args.output.parent.mkdir(parents=True,exist_ok=True);args.output.write_text(json.dumps(record,indent=2)+'\n')
    print(json.dumps({k:v for k,v in record['results'].items() if k!='references'}))


if __name__=='__main__':main()
