#!/usr/bin/env python3
"""Checked membership and drawn-geometry diagnostics for RDKit issue 9705."""
import argparse
import hashlib
import json
from pathlib import Path
import random
import sys
import chematic
from rdkit import Chem, RDLogger, rdBase

ROOT=Path(__file__).resolve().parents[2]
FIXTURE=ROOT/'validation/rdkit_issues/stereo/checked-groups.json'


def digest(data):
    return hashlib.sha256(data).hexdigest()


def groups(mol):
    return [(g.GetGroupType().name,[a.GetIdx() for a in g.GetAtoms()],[b.GetIdx() for b in g.GetBonds()]) for g in mol.GetStereoGroups()]


def inchi(smiles):
    mol=Chem.MolFromSmiles(smiles)
    assert mol is not None
    return Chem.MolToInchi(mol)


def orders(n):
    rng=random.Random(9705)
    result=[list(range(n)),list(reversed(range(n)))]
    while len(result)<24:
        order=list(range(n));rng.shuffle(order)
        if order not in result:result.append(order)
    return result


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,default=ROOT/'validation/results/rdkit-checked-stereo-groups-2026-10-15.json')
    args=parser.parse_args()
    if rdBase.rdkitVersion!='2026.09.1':raise RuntimeError(f'expected RDKit 2026.09.1, got {rdBase.rdkitVersion}')
    RDLogger.DisableLog('rdApp.warning')
    fixtures=json.loads(FIXTURE.read_text())
    checked=[]
    for case in fixtures['cases']:
        base=chematic.from_smiles(case['smiles']);coordinates=base.rdkit_2d_coords()
        expected_inchi=inchi(case['smiles'])
        assert inchi(base.smiles)==expected_inchi
        for kind,number,suffix in [('abs',0,'ABS'),('and',1,'RAC1'),('or',1,'REL1')]:
            original=base.with_stereo_groups([(kind,number,case['atoms'],case['bonds'])])
            for order in orders(base.atom_count):
                inverse={old:new for new,old in enumerate(order)}
                out=original.renumber_atoms(order)
                expected=[(kind,number,[inverse[a] for a in case['atoms']],case['bonds'])]
                assert out.stereo_groups==expected
                block=out.to_mol_v3000([coordinates[old] for old in order],strict=True)
                atom_ids=' '.join(str(a+1) for a in expected[0][2])
                bond_ids=' '.join(str(b+1) for b in case['bonds'])
                assert f'MDLV30/STEB{suffix} BONDS=({len(case["bonds"])} {bond_ids})' in block
                if case['atoms']:assert f'MDLV30/STE{suffix} ATOMS=({len(case["atoms"])} {atom_ids})' in block
                reread=chematic.from_mol_v3000(block)
                assert reread.stereo_groups==expected
                assert inchi(out.smiles)==inchi(reread.smiles)==expected_inchi
                checked.append({'case':case['id'],'kind':kind,'order':order,'membership':expected,'drawn_inchi':expected_inchi,'block_sha256':digest(block.encode())})
    upstream=[]
    for case in (fixtures['cases'][0],fixtures['cases'][2]):
        mol=Chem.RWMol(Chem.MolFromSmiles(case['smiles']))
        mol.SetStereoGroups([Chem.CreateStereoGroup(Chem.StereoGroupType.STEREO_AND,mol,case['atoms'],case['bonds'])])
        mol=mol.GetMol()
        before=groups(mol)
        block=Chem.MolToV3KMolBlock(mol)
        reread=Chem.MolFromMolBlock(block)
        assert reread is not None
        after=groups(reread)
        assert sum(len(g[2]) for g in before)==len(case['bonds'])
        assert sum(len(g[2]) for g in after)==0, 'upstream behavior changed; inspect before rebaseline'
        canonical=groups(Chem.CanonicalizeStereoGroups(mol)) if case['atoms'] else None
        if canonical is not None:assert sum(len(g[2]) for g in canonical)==0
        upstream.append({'case':case['id'],'before':before,'after_V3000':after,'after_CanonicalizeStereoGroups':canonical,'collection_lines':[s for s in block.splitlines() if 'MDLV30/' in s],'block_sha256':digest(block.encode())})
    paths=['crates/chematic-core/src/molecule/checked_stereo_groups.rs','crates/chematic-py/src/mol_methods/stereo_groups.rs','crates/chematic-py/python/chematic/__init__.pyi','crates/chematic-core/tests/checked_stereo_groups.rs','crates/chematic-py/tests/test_stereo_groups.py','scripts/validation/validate_checked_stereo_groups.py']
    record={
        'schema_version':1,'status':'local-source-candidate','issue':'https://github.com/rdkit/rdkit/issues/9705',
        'scope':'checked enhanced-group construction and inspection; existing STE/STEB V3000 interchange is reused, not newly implemented',
        'base':'cac480fb snapshot plus preserved earlier RDKit/Open Babel candidates; this patch adds the checked grouping API',
        'runtime':{'python':sys.version.split()[0],'rdkit':rdBase.rdkitVersion,'chematic_extension_sha256':digest(Path(chematic.__file__).resolve().read_bytes())},
        'source_sha256':{s:digest((ROOT/s).read_bytes()) for s in paths},
        'fixture':{'path':str(FIXTURE.relative_to(ROOT)),'sha256':digest(FIXTURE.read_bytes()),**fixtures['provenance']},
        'results':{'cases':len(checked),'membership_mismatches':0,'drawn_configuration_mismatches':0,'cases_sha256':digest(json.dumps(checked,sort_keys=True,separators=(',',':')).encode()),'kinds':['abs','and','or'],'structures':len(fixtures['cases']),'numberings_per_kind_structure':24},
        'rdkit_reference':{'results':upstream,'note':'RDKit is an independent drawn-geometry/InChI reference and a defect diagnostic; its lossy enhanced-group output is not the membership oracle'},
        'limits':'Group structure only; no automatic stereo perception, inversion, isomer enumeration, canonical enhanced-group identity or new CXSMILES syntax. Same-kind same-number STE/STEB correlation reuses the existing chematic reader convention; BIOVIA examples do not explicitly settle mixed namespace correlation. Plain SMILES and other graph projections do not encode these groups.',
    }
    args.output.parent.mkdir(parents=True,exist_ok=True);args.output.write_text(json.dumps(record,indent=2)+'\n')
    print(json.dumps(record['results']))


if __name__=='__main__':main()
