#!/usr/bin/env python3
"""Check every returned anchored witness independently of the search engine.

Rust snapshots provide indexed graphs and full 0D frames. A separate Python
verifier checks bijection/colors/edges/odd target and all other frame parities.
RDKit independently checks atom/bond properties and supplies E/Z reference
facts. Matching labels/canonical strings is never used to certify redundancy.
"""
from __future__ import annotations
import argparse
import gzip
import hashlib
import json
import random
import subprocess
import sys
from pathlib import Path
import chematic
from rdkit import Chem,rdBase

ROOT=Path(__file__).resolve().parents[1]
SENTINEL=2**32-1
COLORS=['element','isotope','charge','hydrogens','aromatic','map','tag']
ORDERS={'Single':1.0,'Up':1.0,'Down':1.0,'Double':2.0,'Triple':3.0,'Aromatic':1.5}


def odd(values):
    return sum(a>b for i,a in enumerate(values) for b in values[i+1:])%2==1


def reference(case):
    p=Chem.SmilesParserParams();p.removeHs=False
    m=Chem.MolFromSmiles(case['smiles'],p);assert m is not None
    return Chem.RenumberAtoms(m,case['order'])


def rdkit_graph(snapshot,m):
    assert len(snapshot['atoms'])==m.GetNumAtoms()
    assert len(snapshot['bonds'])==m.GetNumBonds()
    for data,a in zip(snapshot['atoms'],m.GetAtoms()):
        assert [data[k] for k in COLORS[:-1]]==[a.GetAtomicNum(),a.GetIsotope(),a.GetFormalCharge(),a.GetTotalNumHs(),a.GetIsAromatic(),a.GetAtomMapNum()]
    expected={frozenset([b.GetBeginAtomIdx(),b.GetEndAtomIdx()]):b.GetBondTypeAsDouble() for b in m.GetBonds()}
    actual={frozenset([b['begin'],b['end']]):ORDERS[b['order']] for b in snapshot['bonds']}
    assert expected==actual


def geometry(m):
    facts={}
    for b in m.GetBonds():
        stereo=str(b.GetStereo())
        if stereo not in ('STEREOE','STEREOZ','STEREOTRANS','STEREOCIS'):
            continue
        ends=(b.GetBeginAtomIdx(),b.GetEndAtomIdx());refs=tuple(b.GetStereoAtoms())
        assert len(refs)==2
        facts[frozenset(ends)]=(dict(zip(ends,refs)),stereo in ('STEREOE','STEREOTRANS'))
    return facts


def verify_witness(row,m):
    a,b=row['source'],row['inverted'];w=row['report']['witness'];target=row['report']['atom_index']
    assert sorted(w)==list(range(len(a['atoms']))) and w[target]==target
    target_edges={frozenset([e['begin'],e['end']]):ORDERS[e['order']] for e in b['bonds']}
    for source,dest in enumerate(w):
        x,y=a['atoms'][source],b['atoms'][dest]
        assert [x[k] for k in COLORS]==[y[k] for k in COLORS]
        if x['frame'] is None:
            assert y['frame'] is None
            continue
        assert len(x['frame'])==len(y['frame'])==4
        mapped=[SENTINEL if n==SENTINEL else w[n] for n in x['frame']]
        order=[y['frame'].index(n) for n in mapped]
        assert len(set(order))==4
        assert (x['chirality']==y['chirality']) != odd(order)
        if source==target:assert odd(order)
    for e in a['bonds']:
        assert target_edges[frozenset([w[e['begin']],w[e['end']]])]==ORDERS[e['order']]
    # E/Z facts originate from RDKit, not the native helper used by the search.
    facts=geometry(m)
    for ends,(refs,trans) in facts.items():
        mapped=frozenset(w[n] for n in ends);target_refs,target_trans=facts[mapped]
        changed=sum(w[outer]!=target_refs[w[center]] for center,outer in refs.items())%2
        assert trans==(target_trans ^ bool(changed))
    return len(facts)


def verify_output(row):
    a,b=row['source'],row['output'];target=row['report']['atom_index'];removed=row['report']['removed']
    assert a['bonds']==b['bonds']
    assert len(a['atoms'])==len(b['atoms'])
    for i,(x,y) in enumerate(zip(a['atoms'],b['atoms'])):
        expected=x.copy()
        if removed:expected['cip']='None'
        if removed and i==target:expected.update(chirality='None',frame=None)
        assert y==expected


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--native-dump',required=True,type=Path)
    p.add_argument('--wasm-module',required=True,type=Path)
    p.add_argument('--output',type=Path,default=ROOT/'validation/results/tetrahedral-cleanup-2026-10-16')
    args=p.parse_args();assert rdBase.rdkitVersion=='2026.09.1'
    fixture=ROOT/'validation/rdkit_issues/tetrahedral-cleanup/cases.json'
    cases=[];identity=[]
    for c in json.loads(fixture.read_text())['cases']:
        for mirrored in [False,True]:
            s=c['smiles'].replace('@@','TEMP').replace('@','@@').replace('TEMP','@') if mirrored else c['smiles']
            n=chematic.from_smiles(s).atom_count
            for seed in [None,0,11,53]:
                order=list(range(n))
                if seed is not None:random.Random(seed).shuffle(order)
                row=dict(c,smiles=s,id=f"{c['id']}-{mirrored}-{seed}",order=order)
                cases.append(row)
                if seed is None:identity.append(row)
    dumped=subprocess.run([str(args.native_dump)],input=json.dumps(cases),capture_output=True,text=True,check=True)
    rows=json.loads(dumped.stdout);assert len(rows)==len(cases)
    witnesses=0;ez_checks=0;python=[]
    for row,c in zip(rows,cases):
        assert row['input']==c
        assert row['report']['removed']==c['removed']
        m=reference(c);rdkit_graph(row['source'],m);verify_output(row)
        if c['removed']:
            ez_checks+=verify_witness(row,m);witnesses+=1
        else:assert row['report']['witness'] is None
        source=chematic.from_smiles(c['smiles']).renumber_atoms(c['order'])
        a=c['order'].index(c['atom_index']);report=source.remove_redundant_tetrahedral_stereo(a)
        result=report.pop('molecule');assert report==row['report']
        assert source.smiles==row['source_smiles'] and result.smiles==row['output_smiles']
        assert source.bond_table==result.bond_table
        python.append(dict(id=c['id'],report=report,output_smiles=result.smiles))
    node=subprocess.run(['node',str(ROOT/'scripts/validate_tetrahedral_cleanup_node.mjs'),str(args.wasm_module)],input=json.dumps(identity),text=True,capture_output=True,check=True)
    wasm=json.loads(node.stdout);by_id={row['input']['id']:row for row in rows}
    for row in wasm:
        native=by_id[row['id']]
        assert row['report']==native['report']
        assert row['source_smiles']==native['source_smiles'] and row['output_smiles']==native['output_smiles']
    assert len(wasm)==len(identity)
    packet=dict(schema='chematic.tetrahedral-cleanup-validation.v1',native=rows,python=python,wasm=wasm)
    args.output.mkdir(parents=True,exist_ok=True)
    data=gzip.compress(json.dumps(packet,sort_keys=True).encode(),mtime=0)
    (args.output/'comparisons.json.gz').write_bytes(data)
    paths=['crates/chematic-core/src/molecule/metadata.rs','crates/chematic-chem/src/tetrahedral_cleanup.rs','crates/chematic-chem/src/meso_stereo_group.rs','crates/chematic-chem/src/meso_stereo_group/graph.rs','crates/chematic-chem/src/meso_stereo_group/search.rs','crates/chematic-py/src/mol_methods/stereo_cleanup.rs','crates/chematic-wasm/src/stereo_cleanup.rs','crates/chematic-chem/examples/tetrahedral_cleanup_dump.rs','scripts/validate_tetrahedral_cleanup.py','scripts/validate_tetrahedral_cleanup_node.mjs','validation/rdkit_issues/tetrahedral-cleanup/cases.json']
    sha=lambda f:hashlib.sha256(Path(f).read_bytes()).hexdigest()
    summary=dict(schema='chematic.tetrahedral-cleanup-summary.v1',source_base=subprocess.check_output(['git','merge-base','HEAD','origin/main'],cwd=ROOT,text=True).strip(),
        source_sha256={s:sha(ROOT/s) for s in paths},
        runtime=dict(rdkit=rdBase.rdkitVersion,python=sys.version.split()[0],python_extension_sha256=sha(chematic.chematic.__file__),wasm_sha256=sha(args.wasm_module.with_name('chematic_wasm_bg.wasm')),native_dump_sha256=sha(args.native_dump)),
        results=dict(native_decisions=len(rows),python_decision_parity=len(python),wasm_decision_parity=len(wasm),independently_verified_witnesses=witnesses,rdkit_atom_bond_graphs_checked=len(rows),rdkit_geometric_facts_checked=ez_checks,failures=0),
        packet_sha256=hashlib.sha256(data).hexdigest(),scope='Bounded opt-in anchored 0D graph proof and declaration removal. Not general stereo cleanup, CIP-label correctness, 3D feasibility, performance or release proof.')
    (args.output/'summary.json').write_text(json.dumps(summary,indent=2)+'\n')
    print(json.dumps(summary['results'],indent=2))


if __name__=='__main__':main()
