#!/usr/bin/env python3
"""Open Babel #1464: isotope-separated totals with authored stoichiometry."""
from __future__ import annotations

import argparse
from collections import Counter
from datetime import datetime,timezone
import gzip
import json
from pathlib import Path
import platform
import re
import subprocess
import sys

ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT))
from scripts.checks.check_interchange_preservation import bounded_run,identity,require,strict_json
from scripts.validation.isotope_formula_observer import outcome

WORKER='''
import json,sys
sys.path.insert(0,sys.argv[1])
import chematic,chematic.chematic as native
rows=[]
for c in json.load(sys.stdin):
 m=chematic.from_smiles(c['smiles']);before=m.to_moljson()
 rows.append({'id':c['id'],'labelled':{'ok':True,'formula':m.formula_with_isotopes},'legacy':{'ok':True,'formula':m.formula},'graph_unchanged':m.to_moljson()==before})
print(json.dumps({'rows':rows,'artifacts':[native.__file__,chematic.__file__]}))
'''
REFERENCE='''
import json,sys,rdkit
from rdkit import Chem,rdBase
from rdkit.Chem.rdMolDescriptors import CalcMolFormula
rdBase.DisableLog('rdApp.*')
rows=[]
for c in json.load(sys.stdin):
 m=Chem.MolFromSmiles(c['smiles'])
 r={'ok':True,'formula':CalcMolFormula(m,separateIsotopes=True,abbreviateHIsotopes=False)} if m is not None else {'ok':False,'error':'MolFromSmiles returned None'}
 rows.append({'id':c['id'],'labelled':r})
print(json.dumps({'version':rdkit.__version__,'rows':rows,'artifacts':[Chem.rdchem.__file__,rdBase.__file__]}))
'''


def run(command,data=b''):
    result=bounded_run([str(x) for x in command],data,timeout=30,max_output=1_048_576)
    require(result['resource_status']=='completed','process resource budget exceeded')
    return {**result,**{k:result[k].decode('utf8',errors='replace') for k in ('stdout','stderr')}}


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    for name in ['python-root','rdkit-python','wasm-module']:parser.add_argument('--'+name,type=Path,required=True)
    parser.add_argument('--obabel',type=Path,default=Path('/opt/homebrew/bin/obabel'))
    parser.add_argument('--output',type=Path,default=ROOT/'validation/results/openbabel-1464-2026-10-24.json')
    args=parser.parse_args();fixture=ROOT/'validation/openbabel-1464-2026-10-24/cases.json'
    cases=strict_json(fixture.read_bytes())['cases'];counts=Counter();raw={};artifacts={};records=[];versions={}
    for surface,command in [
        ('python',[sys.executable,'-c',WORKER,args.python_root]),
        ('wasm',['node',ROOT/'scripts/validation/validate_isotope_formula_node.mjs',args.wasm_module]),
        ('rdkit',[args.rdkit_python,'-c',REFERENCE]),
    ]:
        result=run(command,json.dumps(cases).encode());require(result['exit_code']==0,f'{surface} failed: {result}')
        raw[surface]=result;payload=strict_json(result['stdout'].encode());rows=payload['rows']
        artifacts[surface]=[identity(Path(p)) for p in payload.get('artifacts',[])]
        if 'version' in payload:versions[surface]=payload['version']
        require(len(rows)==len(cases),'worker dropped cases')
        for case,row in zip(cases,rows):
            require(case['id']==row['id'],'worker reordered cases')
            for lane in ['labelled','legacy']:
                if lane not in row:continue
                observed=outcome(case,row[lane]);counts[f'{surface}.{lane}.{observed["status"]}']+=1
                if lane=='labelled':
                    if surface in ('python','wasm') or row[lane]['ok']:
                        require(observed['status']=='composition_match',f'{surface} {case["id"]}: {observed}')
                        spelling_match=row[lane]['formula']==case['expected_formula']
                        if surface in ('python','wasm'):require(spelling_match,'canonical isotope formula spelling changed')
                        elif not spelling_match:counts[f'{surface}.labelled.spelling_difference']+=1
                records.append(dict(surface=surface,lane=lane,id=case['id'],observation=observed))
            if surface in ('python','wasm'):require(row['graph_unchanged'] is True,'formula mutated graph')
    for case in [c for c in cases if c['comparison_lane']=='common']:
        result=run([args.obabel,'-ismi','-omolreport'],(case['smiles']+'\n').encode())
        formulas=re.findall(r'^FORMULA:\s*(\S+)\s*$',result['stdout'],flags=re.M)
        envelope=dict(ok=True,formula=formulas[0]) if result['exit_code']==0 and len(formulas)==1 else dict(ok=False,error=result['stderr'] or 'missing/extra formula records')
        observed=outcome(case,envelope);counts[f'openbabel.common.{observed["status"]}']+=1
        if observed.get('element_totals_match') and observed.get('charge_match') and not observed.get('isotope_counts_match'):counts['openbabel.common.isotope_difference_only']+=1
        records.append(dict(surface='openbabel',id=case['id'],process=result,formula=envelope,observation=observed))
    base='38b4c92c430d9b3bfe70e47e97293b9017ced990';legacy='crates/chematic-chem/src/descriptors/composition.rs'
    require(subprocess.check_output(['git','show',f'{base}:{legacy}'],cwd=ROOT)==(ROOT/legacy).read_bytes(),'legacy composition implementation differs from base main')
    packet=args.output.with_suffix('.packet.json.gz');packet.write_bytes(gzip.compress(json.dumps(dict(worker_outputs=raw,records=records),indent=2,allow_nan=False).encode(),mtime=0))
    sources=['Cargo.lock','crates/chematic-chem/src/isotope_formula.rs','crates/chematic-chem/src/lib.rs',legacy,'crates/chematic-py/src/mol_methods/identity.rs','crates/chematic-py/python/chematic/__init__.pyi','crates/chematic-wasm/src/lib.rs','scripts/validation/validate_isotope_formula.py','scripts/validation/validate_isotope_formula_node.mjs','scripts/validation/isotope_formula_observer.py']
    versions.update(python=sys.version,rust=run(['rustc','--version']),node=run(['node','--version']),openbabel=run([args.obabel,'-V']))
    report=dict(schema='chematic.isotope-formula-validation.v1',utc=datetime.now(timezone.utc).isoformat(),platform=platform.platform(),candidate_base=base,candidate_state='working tree, not a published release',scope='Isotope-resolved totals plus hydrogen and net charge. No isotope abundance, validity, positions, radical-state or structural-identity claim.',fixture=identity(fixture),upstream=identity(ROOT/'validation/openbabel-1464-2026-10-24/upstream.json'),packet=identity(packet),counts=dict(sorted(counts.items())),sources=[identity(ROOT/p) for p in sources],artifacts=artifacts,wasm_artifacts=[identity(args.wasm_module),identity(args.wasm_module.with_name('chematic_wasm_bg.wasm'))],openbabel_artifact=identity(args.obabel.resolve()),versions=versions,legacy_control=dict(scope='Existing isotope-blind formula in the same candidate runtime; complete composition implementation equals pinned main. Not a separate main/release executable.',source_commit=base),comparison_notes=['Expected counts/charge/formula are authored before execution; RDKit is a cross-check, not the oracle.','Open Babel common lane contains 42 same-SMILES molreport pipelines; D/T and repeated-sign charge aliases normalize to equivalent composition. Existing charge behavior is not scored as a new implementation gap.','Five wildcard/mass-label boundary controls are candidate-only; they are not extra Open Babel failures.'],process_limits=dict(seconds=30,combined_pipe_bytes=1_048_576,memory='not imposed'))
    args.output.write_text(json.dumps(report,indent=2,allow_nan=False)+'\n');print(json.dumps(report['counts'],indent=2))


if __name__=='__main__':main()
