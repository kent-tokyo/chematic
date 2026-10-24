#!/usr/bin/env python3
"""RDKit #9735: strict valence refusals and unchanged valid molecule controls."""
from __future__ import annotations

import argparse
from collections import Counter
from datetime import datetime, timezone
import gzip
import json
from pathlib import Path
import platform
import subprocess
import sys

ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT))
from scripts.checks.check_interchange_preservation import bounded_run, identity, require, strict_json

WORKER='''
import json,sys
sys.path.insert(0,sys.argv[1])
import chematic
rows=[]
for c in json.load(sys.stdin):
 r={'id':c['id']}
 try:
  mol=chematic.from_smiles(c['smiles'])
  try:r['canonical']={'ok':True,'value':mol.rdkit_smiles}
  except ValueError as e:r['canonical']={'ok':False,'error':str(e)}
  if c.get('lane')!='regression':
   r['problems']=chematic.rdkit_detect_chemistry_problems(c['smiles'])
   try:r['checked']={'ok':True,'value':chematic.rdkit_smiles_checked(mol)}
   except ValueError as e:r['checked']={'ok':False,'error':str(e)}
 except ValueError as e:r['canonical']={'ok':False,'error':str(e)}
 rows.append(r)
import chematic.chematic as native
print(json.dumps({'rows':rows,'artifacts':[native.__file__,chematic.__file__]},allow_nan=False))
'''
RDKIT_WORKER='''
import json,sys,rdkit
from rdkit import Chem,rdBase
rdBase.DisableLog('rdApp.*')
rows=[]
for c in json.load(sys.stdin):
 m=Chem.MolFromSmiles(c['smiles']);r={'id':c['id']}
 r['canonical']={'ok':True,'value':Chem.MolToSmiles(m)} if m is not None else {'ok':False,'error':'MolFromSmiles returned None'}
 raw=Chem.MolFromSmiles(c['smiles'],sanitize=False)
 r['problems']=[[p.GetType(),[p.GetAtomIdx()] if hasattr(p,'GetAtomIdx') else list(p.GetAtomIndices())] for p in Chem.DetectChemistryProblems(raw)]
 rows.append(r)
print(json.dumps({'version':rdkit.__version__,'rows':rows,'artifacts':[Chem.rdchem.__file__,rdBase.__file__,rdkit.__file__]}))
'''


def status(result):
    require(type(result) is dict and type(result.get('ok')) is bool,'invalid result envelope')
    if result['ok']:
        require(isinstance(result.get('value'),str) and result['value'],'missing canonical value')
        return 'accepted'
    require('value' not in result and isinstance(result.get('error'),str),'partial refusal')
    message=result['error']
    return 'valence_error' if 'Explicit valence for atom #' in message and 'greater than permitted' in message else 'other_refusal'


def verify_candidate(case,row,diagnostics):
    outcome=status(row['canonical'])
    require(outcome==('accepted' if case['expected']=='accept' else 'valence_error'),f'{case["id"]}: {outcome}')
    if diagnostics:
        if case['expected']=='valence_error':
            require(['AtomValenceException',[case['atom']]] in row['problems'],'missing atom-indexed valence diagnosis')
            require(status(row['checked'])=='valence_error','checked writer returned/wrongly refused')
        else:
            require(row['problems']==[] and row['checked']==row['canonical'],'valid molecule changed through checked writer')
    return outcome


def run(command,data=b''):
    result=bounded_run([str(x) for x in command],data,timeout=120,max_output=8*1024*1024)
    require(result['resource_status']=='completed','process resource budget exceeded')
    return {**result,**{k:result[k].decode('utf8',errors='replace') for k in ('stdout','stderr')}}


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    for name in ['python-root','baseline-python-root','rdkit-python','wasm-module']:
        parser.add_argument('--'+name,type=Path,required=True)
    parser.add_argument('--output',type=Path,default=ROOT/'validation/results/rdkit-9735-2026-10-24.json')
    args=parser.parse_args()
    fixture=ROOT/'validation/rdkit-9735-2026-10-24/cases.json';cases=strict_json(fixture.read_bytes())['cases']
    corpus=ROOT/'scripts/fixtures/chembl_accuracy_corpus_4999.smi'
    sample=[dict(id=f'corpus-{i:04}',smiles=s,lane='regression') for i,s in enumerate(corpus.read_text().splitlines()[:1000])]
    raw={};results={};counts=Counter()
    for surface,command,inputs in [
        ('candidate_python',[sys.executable,'-c',WORKER,args.python_root],cases+sample),
        ('baseline_python',[sys.executable,'-c',WORKER,args.baseline_python_root],cases+sample),
        ('rdkit',[args.rdkit_python,'-c',RDKIT_WORKER],cases),
        ('candidate_wasm',['node',ROOT/'scripts/validation/validate_rdkit_issue9735_node.mjs',args.wasm_module],cases),
    ]:
        result=run(command,json.dumps(inputs).encode());require(result['exit_code']==0,f'{surface} worker failed: {result}')
        raw[surface]=result;payload=strict_json(result['stdout'].encode());rows=payload['rows'];results[surface]=payload
        require(len(rows)==len(inputs),'worker dropped rows')
        for case,row in zip(inputs,rows):
            require(row['id']==case['id'],'worker reordered rows')
            lane=case.get('lane',case.get('expected'))
            outcome=verify_candidate(case,row,surface=='candidate_python') if surface.startswith('candidate') and lane!='regression' else status(row['canonical'])
            counts[f'{surface}.{lane}.{outcome}']+=1
    a,b=results['candidate_python']['rows'],results['baseline_python']['rows']
    for case,new,old in zip(cases+sample,a,b):
        if case.get('lane')=='regression' or case['expected']=='accept':
            require(new['canonical']==old['canonical'],f'control changed: {case["id"]}')
    old_commit='a4a2e732f50a892c8accf29273d8481f3a005e0a';base='38b4c92c430d9b3bfe70e47e97293b9017ced990'
    def git_file(ref,path):return subprocess.check_output(['git','show',f'{ref}:{path}'],cwd=ROOT)
    trees=[subprocess.check_output(['git','rev-parse',f'{ref}:crates/chematic-smiles'],cwd=ROOT).decode().strip() for ref in (base,old_commit)]
    require(trees[0]==trees[1],'baseline SMILES tree differs from pinned main')
    attested=strict_json(git_file(old_commit,'validation/results/orca-checked-2026-10-24.json'))
    for key,path in zip(['python_extension','python_wrapper'],results['baseline_python']['artifacts']):
        require(identity(Path(path))['sha256']==attested['artifacts'][key]['sha256'],'baseline artifact attestation failed')
    packet=args.output.with_suffix('.packet.json.gz');packet.write_bytes(gzip.compress(json.dumps(raw,indent=2,allow_nan=False).encode(),mtime=0))
    sources=['Cargo.lock','crates/chematic-smiles/src/rdkit/mol.rs','crates/chematic-smiles/src/rdkit/mod.rs','scripts/validation/validate_rdkit_issue9735.py','scripts/validation/validate_rdkit_issue9735_node.mjs']
    report=dict(schema='chematic.rdkit-9735-validation.v1',utc=datetime.now(timezone.utc).isoformat(),platform=platform.platform(),candidate_base=base,candidate_state='working tree, not a published release',upstream=identity(ROOT/'validation/rdkit-9735-2026-10-24/upstream.json'),fixture=identity(fixture),packet=identity(packet),counts=dict(sorted(counts.items())),sources=[identity(ROOT/p) for p in sources],baseline=dict(commit=old_commit,smiles_tree=trees[0],scope='Retained previous PR candidate, SMILES source tree equal to base main; artifact attested against its ORCA report. Not a separately rebuilt main/release executable.'),regression=dict(source=identity(corpus),selection='first 1000 physical lines, fixed existing development corpus; no holdout claim',count=len(sample)),artifacts={k:[identity(Path(p)) for p in v.get('artifacts',[])] for k,v in results.items()},wasm_artifacts=[identity(args.wasm_module),identity(args.wasm_module.with_name('chematic_wasm_bg.wasm'))],versions=dict(rdkit=results['rdkit']['version'],python=sys.version,rust=run(['rustc','--version']),node=run(['node','--version'])),scope='Strict compatibility sanitization/diagnostics and SMILES writing. Ordinary graph parsing and non-strict caches remain unchanged; this intentionally differs from RDKit on #9735 invalid inputs.',process_limits=dict(seconds=120,combined_pipe_bytes=8*1024*1024,memory='not imposed'))
    args.output.write_text(json.dumps(report,indent=2,allow_nan=False)+'\n');print(json.dumps(report['counts'],indent=2))


if __name__=='__main__':main()
