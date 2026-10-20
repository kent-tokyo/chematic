#!/usr/bin/env python3
"""Ordered radius bonds and minimum atom distances against pinned RDKit."""
import argparse
import hashlib
import itertools
import json
from pathlib import Path
import random

import chematic
from rdkit import Chem, rdBase, RDLogger

ROOT=Path(__file__).resolve().parents[2]
FIXTURE=ROOT/'validation/rdkit_issues/atom-environment/fixtures.json'


def digest(data):return hashlib.sha256(data).hexdigest()


def orders(n):
    out=[list(range(n))]
    for row in (list(reversed(range(n))),random.Random(9708+n).sample(range(n),n)):
        if row not in out:out.append(row)
    return out


def variant(case,added_hs,order):
    native=chematic.from_smiles(case['smiles']) if case['smiles'] else chematic.RWMol().GetMol()
    if added_hs:
        native=native.rdkit_add_hydrogens();reference=Chem.AddHs(Chem.MolFromSmiles(case['smiles']))
    else:
        params=Chem.SmilesParserParams();params.removeHs=False
        reference=Chem.MolFromSmiles(case['smiles'],params)
    assert reference is not None and native.atom_count==reference.GetNumAtoms()
    native=native.renumber_atoms(order)
    reference=Chem.RenumberAtoms(reference,order) if order else Chem.Mol(reference)
    if not added_hs:reference=Chem.RemoveHs(reference)
    return native,reference


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,default=ROOT/'validation/results/rdkit-atom-environment-2026-10-15.json')
    args=parser.parse_args();assert rdBase.rdkitVersion=='2026.09.1'
    RDLogger.DisableLog('rdApp.*')
    manifest=json.loads(FIXTURE.read_text())
    assert manifest['sha256']==digest(json.dumps(manifest['cases'],sort_keys=True,separators=(',',':')).encode())
    checksum=hashlib.sha256();calls=variants=bonds=distances=empty=radius_zero=hs=permuted=0
    records=[]
    for case in manifest['cases']:
        assert case['sha256']==digest(case['smiles'].encode())
        checked=0
        for added_hs in (False,True):
            base=chematic.from_smiles(case['smiles']) if case['smiles'] else chematic.RWMol().GetMol()
            if added_hs:base=base.rdkit_add_hydrogens()
            for order in orders(base.atom_count):
                native,reference=variant(case,added_hs,order);n=reference.GetNumAtoms()
                assert native.rdkit_num_atoms()==n
                assert native.rdkit_adjacency_matrix()==Chem.GetAdjacencyMatrix(reference,force=True).tolist()
                original=native.smiles;native.SetIntProp('environment_sentinel',9708)
                for root,radius,use_hs,enforce in itertools.product(range(n),[*range(min(10,n+3)),2**32-1],[False,True],[False,True]):
                    atom_map={999999:123,root:456}
                    want=list(Chem.FindAtomEnvironmentOfRadiusN(reference,radius,root,useHs=use_hs,enforceSize=enforce,atomMap=atom_map))
                    got=native.rdkit_atom_environment(radius,root,use_hs=use_hs,enforce_size=enforce)
                    assert got['bond_indices']==want,(case['id'],added_hs,order,root,radius,use_hs,enforce,got,want)
                    assert got['atom_distances']==atom_map,(case['id'],root,radius,got,atom_map)
                    assert len(want)==len(set(want))
                    assert all(i<reference.GetNumBonds() for i in want)
                    checksum.update(json.dumps([want,sorted(atom_map.items())],separators=(',',':')).encode())
                    calls+=1;checked+=1;bonds+=len(want);distances+=len(atom_map)
                    empty+=not bool(want);radius_zero+=radius==0;hs+=use_hs
                    permuted+=order!=list(range(len(order)))
                assert native.smiles==original and native.GetIntProp('environment_sentinel')==9708
                variants+=1
        records.append({'id':case['id'],'status':'passed','comparisons':checked})
    # Record the surprising cyclic boundary separately from the chemistry corpus.
    ring=Chem.MolFromSmiles('C1CCC1');ring3={};ring4={}
    assert list(Chem.FindAtomEnvironmentOfRadiusN(ring,3,0,atomMap=ring3))==[0,3,1,2]
    assert list(Chem.FindAtomEnvironmentOfRadiusN(ring,4,0,atomMap=ring4))==[] and ring4=={}
    report={'schema_version':1,'status':'local-source-candidate; uncommitted and unpushed',
        'oracle':{'name':'RDKit','version':rdBase.rdkitVersion,'settings':'Every root; useHs False/True; enforceSize False/True; radii 0..min(9,n+2) and UINT32_MAX; H removal/addition and three distinct atom orders where available; seeded atomMap'},
        'fixture_manifest':str(FIXTURE.relative_to(ROOT)),'fixture_sha256':manifest['sha256'],
        'results':{'fixtures':len(records),'graph_variants':variants,'exact_comparisons':calls,
            'ordered_bond_ids_compared':bonds,'atom_distances_compared':distances,'empty_bond_results':empty,
            'radius_zero_cases':radius_zero,'use_hs_true_cases':hs,'permuted_cases':permuted,
            'unexpected_mismatches':0,'normalized_output_sha256':checksum.hexdigest()},
        'records':records,'cyclic_boundary':{'C1CCC1_root0_radius3_bonds':[0,3,1,2],'radius3_atom_map':ring3,'radius4_atom_map':ring4,
            'explanation':'A final pending-duplicate frontier satisfies enforceSize even if no new bond is added'},
        'scope':['Sanitized RDKit atom and bond indices','Rust/Python only','No caller dictionary mutation; returns a new map',
            'Input <=1024 atoms/16384 bonds','Default <=1000000 visits/100000 queued entries per frontier',
            'Query/zero bonds and other unsupported model inputs fail explicitly']}
    args.output.parent.mkdir(parents=True,exist_ok=True);args.output.write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(report['results'],indent=2))


if __name__=='__main__':main()
