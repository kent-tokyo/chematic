#!/usr/bin/env python3
"""Pinned direct-configuration oracle for correlated enhanced-group alternatives."""
import argparse
import hashlib
import json
from pathlib import Path
import random
import sys
import chematic
from rdkit import Chem, RDLogger, rdBase
from rdkit.Chem.EnumerateStereoisomers import EnumerateStereoisomers, StereoEnumerationOptions

ROOT = Path(__file__).resolve().parents[1]
FIXTURE = ROOT/'validation/rdkit_issues/stereo/group-alternatives.json'


def digest(data):
    return hashlib.sha256(data).hexdigest()


def parse(s):
    p=Chem.SmilesParserParams();p.removeHs=False
    m=Chem.MolFromSmiles(s,p)
    assert m is not None
    return m


def inchi(smiles):
    return Chem.MolToInchi(parse(smiles))


def reference(case, groups, mask, native):
    m=parse(case['smiles'])
    relative=[g for g in groups if g[0]!='abs']
    for bit,(_,_,atoms,bonds) in enumerate(relative):
        if not mask & (1<<bit):continue
        for index in atoms:m.GetAtomWithIdx(index).InvertChirality()
        for index in bonds:
            a,b,*_=native.bond_table[index]
            bond=m.GetBondBetweenAtoms(a,b)
            assert bond is not None
            stereo=bond.GetStereo()
            assert stereo in (Chem.BondStereo.STEREOE,Chem.BondStereo.STEREOZ)
            bond.SetStereo(Chem.BondStereo.STEREOZ if stereo==Chem.BondStereo.STEREOE else Chem.BondStereo.STEREOE)
    for bond in m.GetBonds():bond.SetBondDir(Chem.BondDir.NONE)
    Chem.SetDoubleBondNeighborDirections(m)
    smiles=Chem.MolToSmiles(m)
    expected=Chem.MolToInchi(m)
    assert inchi(smiles)==expected
    return expected


def orders(n):
    rng=random.Random(9705);out=[list(range(n)),list(reversed(range(n)))]
    while len(out)<24:
        order=list(range(n));rng.shuffle(order)
        if order not in out:out.append(order)
    return out


def normalize_bonds(bonds):
    return sorted((min(a,b),max(a,b),str(kind),aromatic) for a,b,kind,aromatic in bonds)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,default=ROOT/'validation/results/rdkit-stereo-group-alternatives-2026-10-15.json')
    args=parser.parse_args()
    if rdBase.rdkitVersion!='2026.09.1':raise RuntimeError(f'expected RDKit 2026.09.1, got {rdBase.rdkitVersion}')
    Chem.SetUseLegacyStereoPerception(True);RDLogger.DisableLog('rdApp.warning')
    fixture=json.loads(FIXTURE.read_text());checked=[];references=[]
    assert digest(json.dumps(fixture['cases'],sort_keys=True,separators=(',',':')).encode())==fixture['sha256']
    for case in fixture['cases']:
        base=chematic.from_smiles(case['smiles']);rd=parse(case['smiles'])
        assert base.atom_count==rd.GetNumAtoms()
        assert [(a[1],a[2],a[3]) for a in base.atom_table]==[(a.GetAtomicNum(),a.GetFormalCharge(),a.GetIsAromatic()) for a in rd.GetAtoms()]
        assert normalize_bonds(base.bond_table)==normalize_bonds([(b.GetBeginAtomIdx(),b.GetEndAtomIdx(),b.GetBondType(),b.GetIsAromatic()) for b in rd.GetBonds()])
        assert inchi(base.smiles)==Chem.MolToInchi(rd)
        group_orders=[list(range(len(case['groups'])))]
        reverse=list(reversed(group_orders[0]))
        if reverse!=group_orders[0]:group_orders.append(reverse)
        for group_order in group_orders:
            specs=[case['groups'][i] for i in group_order]
            relative_indices=[i for i,g in enumerate(specs) if g[0]!='abs']
            expected=[reference(case,specs,mask,base) for mask in range(1<<len(relative_indices))]
            references.append({'id':case['id'],'group_order':group_order,'representatives':len(expected),'distinct_drawn_inchi':len(set(expected)),'mask_inchi':expected})
            for order in orders(base.atom_count):
                source=base.with_stereo_groups([tuple(g) for g in specs]).renumber_atoms(order)
                source.SetIntProp('oracle_trial',9705)
                before=source.smiles;before_groups=source.stereo_groups
                groups,molecules=source.enumerate_stereo_group_alternatives()
                assert groups==relative_indices and len(molecules)==len(expected)
                for mask,(m,wanted) in enumerate(zip(molecules,expected)):
                    assert m.stereo_groups==before_groups
                    assert m.bond_table==source.bond_table and m.atom_table==source.atom_table
                    assert m.GetIntProp('oracle_trial')==9705
                    assert inchi(m.smiles)==wanted
                    # Applying the selected group inversions again restores mask 0.
                    restored=m
                    for bit,index in enumerate(groups):
                        if mask&(1<<bit):restored=restored.invert_stereo_group(index)
                    assert inchi(restored.smiles)==expected[0]
                    checked.append({'id':case['id'],'order':order,'group_order':group_order,'mask':mask,'inchi':wanted})
                assert source.smiles==before and source.stereo_groups==before_groups
    assert all(r["distinct_drawn_inchi"]==1 for r in references if r["id"]=="symmetric-diol-full-alternative")
    # Independently reproduce atom-only flipping by the upstream enumerator.
    case=fixture['cases'][0];native=chematic.from_smiles(case['smiles'])
    rw=Chem.RWMol(parse(case['smiles']));kind={'and':Chem.StereoGroupType.STEREO_AND,'or':Chem.StereoGroupType.STEREO_OR}
    rw.SetStereoGroups([Chem.CreateStereoGroup(kind[k],rw,a,b) for k,_,a,b in case['groups']])
    products=list(EnumerateStereoisomers(rw.GetMol(),options=StereoEnumerationOptions(onlyUnassigned=True,unique=False)))
    observed=[Chem.MolToInchi(p) for p in products]
    intended=[reference(case,case['groups'],mask,native) for mask in range(4)]
    assert len(products)==4 and set(observed)!=set(intended), 'upstream behavior changed; inspect before rebaseline'
    original_rd=parse(case['smiles'])
    fixed_bonds=[]
    for _,_,_,bonds in case['groups']:
        for index in bonds:
            a,b,*_=native.bond_table[index]
            fixed_bonds.append(all(p.GetBondBetweenAtoms(a,b).GetStereo()==original_rd.GetBondBetweenAtoms(a,b).GetStereo() for p in products))
    assert all(fixed_bonds)
    # Atom-only enhanced enumeration is a valid separate RDKit control.
    case=fixture['cases'][2];native=chematic.from_smiles(case['smiles']);rw=Chem.RWMol(parse(case['smiles']))
    kind['abs']=Chem.StereoGroupType.STEREO_ABSOLUTE
    rw.SetStereoGroups([Chem.CreateStereoGroup(kind[k],rw,a,b) for k,_,a,b in case['groups']])
    control={Chem.MolToInchi(p) for p in EnumerateStereoisomers(rw.GetMol(),options=StereoEnumerationOptions(onlyUnassigned=True,unique=False))}
    expected_control={reference(case,case['groups'],mask,native) for mask in range(4)}
    assert control==expected_control
    paths=['crates/chematic-chem/src/stereo_group_alternatives.rs','crates/chematic-py/src/mol_methods/stereo_groups.rs','crates/chematic-chem/tests/stereo_group_alternatives.rs','crates/chematic-py/tests/test_stereo_group_alternatives.py','crates/chematic-py/python/chematic/__init__.pyi','scripts/validate_stereo_group_alternatives.py']
    record={
        'schema_version':1,'status':'local-source-candidate','issue':'https://github.com/rdkit/rdkit/issues/9705',
        'base':'cac480fb source snapshot plus preserved previous candidates through TURBOMOLE; requires checked stereo groups and full group inversion',
        'scope':'bounded full alternatives of AND/OR tetrahedral and geometric groups, preserving membership and every mask; no distinct compound count, meso reduction or canonicalization',
        'runtime':{'python':sys.version.split()[0],'rdkit':rdBase.rdkitVersion,'settings':'legacy stereo; removeHs=False; direct atom/bond inversions with regenerated neighbor directions; native group enumeration is the subject under test','chematic_extension_sha256':digest(Path(chematic.__file__).resolve().read_bytes())},
        'source_sha256':{s:digest((ROOT/s).read_bytes()) for s in paths},
        'fixture':{'path':str(FIXTURE.relative_to(ROOT)),'manifest_sha256':digest(FIXTURE.read_bytes()),**fixture},
        'results':{'comparisons':len(checked),'configuration_mismatches':0,'group_membership_mismatches':0,'indexed_graph_mismatches':0,'input_mutations':0,'restoration_mismatches':0,'fixtures':len(fixture['cases']),'numberings_per_fixture_group_order':24,'group_orders':'stored and reversed (one when equal)','case_sha256':digest(json.dumps(checked,sort_keys=True,separators=(',',':')).encode()),'references':references},
        'upstream_mixed_enumeration_diagnostic':{'products':len(products),'distinct_drawn_inchi':len(set(observed)),'native_expected_inchi':intended,'rdkit_observed_inchi':observed,'all_grouped_double_bonds_remain_fixed':all(fixed_bonds),'sets_differ':True},
        'upstream_atom_only_control':{'products':len(control),'configuration_set_matches':True},
        'limits':{'max_relative_groups':6,'max_representatives':64,'max_atoms':1024,'max_bonds':4096,'no_partial_results':True},
        'limitations':'No undefined stereo assignment, coordinate generation, pure-isomer conversion, group canonicalization, meso reduction, axes, non-tetrahedral members or molecular-identity claim for retained duplicates. Group metadata describes the original record and remains attached to every drawn representative. The same-number mixed grouping rule is the #9705 proposal, not confirmed BIOVIA product interoperability.',
    }
    args.output.parent.mkdir(parents=True,exist_ok=True);args.output.write_text(json.dumps(record,indent=2)+'\n')
    print(json.dumps({k:v for k,v in record['results'].items() if k!='references'}))


if __name__=='__main__':main()
