#!/usr/bin/env python3
"""Pinned RDKit #9701 diagnostic for native carbon stereo assignments."""
import argparse
import hashlib
import json
from pathlib import Path
import random
import sys

import chematic
from rdkit import Chem,rdBase
from rdkit.Chem.EnumerateStereoisomers import EnumerateStereoisomers,StereoEnumerationOptions

ROOT=Path(__file__).resolve().parents[2]
CASES=[
    ('aminoadamantane','NC12CC3CC(CC(C3)C1)C2',False,4),
    ('aminoadamantane-explicit-h','NC12CC3CC(CC(C3)C1)C2',True,4),
    ('inositol','OC1C(O)C(O)C(O)C(O)C1O',False,9),
    ('long-alkane','CCCCCCCCCCCCCCCC',False,1),
    ('long-one-center','CCCCCCCCCCCCCCCC(O)F',False,2),
    ('terminal-halogen-duplicate','FC(F)(Br)Cl',False,1),
    ('terminal-methyl-duplicate','CC(C)(Br)Cl',False,1),
    ('explicit-h-duplicate','[H]C([H])(F)Cl',False,1),
    ('isotopic-h','[H]C([2H])(F)Cl',False,2),
    ('mixed-isotopic-h','[2H]C(F)Cl',False,2),
    ('meso-diol','CC(O)C(O)C',False,3),
    ('fixed-and-unassigned','N[C@@H](C)C(O)F',False,2),
]
BASELINE_SOURCE_SHA='674d63b96202f03a819ed3a87aa80e86d81b1519322919b7bb3cb135a3b716f1'


def digest(data):return hashlib.sha256(data).hexdigest()


def inchi_from_smiles(smiles):
    mol=Chem.MolFromSmiles(smiles)
    assert mol is not None,smiles
    result=Chem.MolToInchi(mol)
    assert result,result
    return result


def orders(n):
    result=[list(range(n)),list(reversed(range(n))),list(range(1,n))+[0]]
    rng=random.Random(9701)
    for _ in range(2):
        order=list(range(n));rng.shuffle(order);result.append(order)
    return list(dict.fromkeys(map(tuple,result)))


def evaluate():
    if rdBase.rdkitVersion!='2026.09.1':raise RuntimeError(f'expected RDKit 2026.09.1, got {rdBase.rdkitVersion}')
    opts=StereoEnumerationOptions(onlyUnassigned=True,unique=True,maxIsomers=1024,tryEmbedding=False)
    fixtures=[];runs=[]
    for identifier,smiles,explicit_h,count in CASES:
        Chem.SetUseLegacyStereoPerception(True)
        reference=list(EnumerateStereoisomers(Chem.MolFromSmiles(smiles),options=opts))
        reference_keys={Chem.MolToInchi(r) for r in reference}
        assert len(reference)==len(reference_keys)==count,(identifier,len(reference),len(reference_keys))
        native=chematic.from_smiles(smiles)
        if explicit_h:native=native.add_hydrogens()
        for order in orders(native.atom_count):
            source=native.renumber_atoms(order)
            forms=source.enumerate_stereoisomers()
            strings=[f.smiles for f in forms]
            assert len(forms)==len(set(strings))==count,(identifier,len(forms))
            keys={inchi_from_smiles(s) for s in strings}
            assert keys==reference_keys,identifier
            assert all(chematic.from_smiles(s).smiles==s for s in strings),identifier
            for form in forms:
                assert form.atom_count==source.atom_count
                assert form.bond_table==source.bond_table
            runs.append({'fixture':identifier,'new_to_old':order,'native_keys':sorted(strings),'inchi_set':sorted(keys)})
        Chem.SetUseLegacyStereoPerception(False)
        modern=list(EnumerateStereoisomers(Chem.MolFromSmiles(smiles),options=opts))
        modern_count=len(modern)
        if identifier=='aminoadamantane':assert modern_count==5
        if identifier=='inositol':assert modern_count==13
        fixtures.append({
            'id':identifier,'smiles':smiles,'explicit_h':explicit_h,'sha256':digest((smiles+'\n').encode()),
            'provenance':{'source':'minimal abstract molecular graphs and original source-order permutation recipe','retrieved':'2026-10-15','issue':'https://github.com/rdkit/rdkit/issues/9701','transformations':'issue skeletons plus independently authored H/isotope/terminal-ligand/chain controls; no upstream dataset or source code copied'},
            'license':'MIT OR Apache-2.0','redistribution':{'status':'allowed','rationale':'minimal abstract structures and original test recipe'},
            'minimization':'one candidate-budget or stereo-preservation condition per small graph',
            'oracle':{'engine':'RDKit','version':'2026.09.1','settings':{'onlyUnassigned':True,'unique':True,'maxIsomers':1024,'tryEmbedding':False,'reference_legacy_stereo':True,'comparison':'exact InChI set after reading native SMILES under legacy perception'}},
            'expected':{'assignment_count':count,'reference_inchi_set':sorted(reference_keys),'count_tolerance':0},
            'native_assignment_count':count,'legacy_reference_count':len(reference),'modern_reference_count':modern_count,'source_numberings':len(orders(native.atom_count)),
        })
    source_files=['crates/chematic-chem/src/stereo.rs','crates/chematic-py/src/mol_methods/transformations.rs']
    return {
        'schema_version':1,'status':'local-source-candidate','issue':'https://github.com/rdkit/rdkit/issues/9701',
        'base':'cac480fb source snapshot plus preserved query-H, PQR, atom-permutation and mmCIF residue candidates; patch isolates native enumeration changes',
        'baseline':{'aminoadamantane_observed_assignment_count':0,'stereo_source_sha256':BASELINE_SOURCE_SHA,'cause':'degree-two CH2 sites counted toward six-candidate cap'},
        'scope':'native saturated-carbon assignments, fixed diagnostic graphs, no 3D feasibility/heteroatom/EZ enumeration or general canonicalization claim; RDKit-compatible enumeration remains separate',
        'runtime':{'python':sys.version.split()[0],'rdkit':rdBase.rdkitVersion,'chematic_extension_sha256':digest(Path(chematic.__file__).resolve().read_bytes())},
        'source_sha256':{s:digest((ROOT/s).read_bytes()) for s in source_files},
        'permutation_recipe':{'identity':True,'reverse':True,'cyclic_shift':1,'random_orders':2,'seed':9701,'deduplicated':True},
        'results':{'fixtures':len(fixtures),'source_numberings':len(runs),'compared_assignments':sum(len(r['native_keys']) for r in runs),'inchi_set_mismatches':0,'native_roundtrip_failures':0,'graph_mismatches':0,'cases_sha256':digest(json.dumps(runs,sort_keys=True,separators=(',',':')).encode())},
        'fixtures':fixtures,
    }


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,default=ROOT/'validation/results/rdkit-stereo-candidates-2026-10-15.json')
    args=parser.parse_args();previous=Chem.GetUseLegacyStereoPerception()
    try:record=evaluate()
    finally:Chem.SetUseLegacyStereoPerception(previous)
    args.output.parent.mkdir(parents=True,exist_ok=True);args.output.write_text(json.dumps(record,indent=2)+'\n')
    print(json.dumps(record['results']))


if __name__=='__main__':main()
