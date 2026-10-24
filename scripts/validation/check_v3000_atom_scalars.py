#!/usr/bin/env python3
"""Atom-field interpretation and explicit refusals across Python/WASM/OB.

Authored single-atom fields are the oracle. RDKit observes emitted SMILES,
not the source dialect's unknown attributes. Refusals are separate from retention.
"""
from __future__ import annotations
import argparse
import gzip
import importlib
import json
from pathlib import Path
import subprocess
import sys

ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT))
from scripts.checks.check_interchange_preservation import bounded_run,identity,require,strict_json

NODE="""
const fs=require('node:fs'),w=require(process.argv[1]);
const p=JSON.parse(fs.readFileSync(0,'utf8'));
const capture=fn=>{try{return {accepted:true,output:fn()};}catch(e){return {accepted:false,error:String(e)};}};
const rows=p.cases.map(c=>{const m=w.mol_from_v3000_block(c.input);try{return {id:c.id,smiles:m.canonical_smiles(),moljson:w.to_moljson(m),common:w.convert_common_format(c.input,'v3000','smiles'),source:w.convert_common_format(c.input,'v3000','v3000')};}finally{m.free();}});
const invalid=p.invalid.map(c=>({id:c.id,direct:capture(()=>{const m=w.mol_from_v3000_block(c.input);try{return m.canonical_smiles();}finally{m.free();}}),common:capture(()=>w.convert_common_format(c.input,'v3000','smiles')),source:capture(()=>w.convert_common_format(c.input,'v3000','v3000'))}));
process.stdout.write(JSON.stringify({rows,invalid}));
"""


def observe(expected,smiles,reader=None):
    if reader is None:
        reader=smiles_fields
    try:
        actual=reader(smiles)
        return dict(readable=True,fields_exact=actual==expected and all(type(actual[key]) is int for key in expected),actual=actual)
    except (ValueError,RuntimeError,KeyError,IndexError,OverflowError) as error:
        return dict(readable=False,fields_exact=False,error=str(error))


def smiles_fields(text):
    from rdkit import Chem
    params=Chem.SmilesParserParams();params.sanitize=False;params.removeHs=False
    mol=Chem.MolFromSmiles(text,params)
    require(mol is not None and mol.GetNumAtoms()==1 and mol.GetNumBonds()==0,'expected one unbonded atom')
    a=mol.GetAtomWithIdx(0)
    return dict(number=a.GetAtomicNum(),charge=a.GetFormalCharge(),isotope=a.GetIsotope(),atom_map=a.GetAtomMapNum())


def moljson_fields_exact(expected,text):
    try:
        value=strict_json(text)
        atoms=value['atoms']
        if len(atoms)!=1 or value['bonds']:
            return False
        a=atoms[0]
        return (a['element']=={1:'H',6:'C',7:'N',8:'O'}[expected['number']]
                and type(a['charge']) is int and a['charge']==expected['charge']
                and (a['isotope'] is None if expected['isotope']==0 else type(a['isotope']) is int and a['isotope']==expected['isotope']))
    except (ValueError,KeyError,TypeError,IndexError):
        return False


def terminal(run):
    return {k:v.decode('utf8',errors='replace') if isinstance(v,bytes) else v for k,v in run.items()}


def capture(fn):
    try:return dict(accepted=True,output=fn())
    except ValueError as error:return dict(accepted=False,error=str(error))


def obabel(text,flags):
    command=['/opt/homebrew/bin/obabel','-imol','-osmi','-xn',*flags]
    run=terminal(bounded_run(command,text.encode()));run['command']=command
    run['outcome']=process_outcome(run)
    return run


def process_outcome(run):
    if run['resource_status']!='completed': return 'resource_limit'
    if run['exit_code'] is not None and run['exit_code']<0: return 'signal'
    if run['exit_code']!=0: return 'error'
    if not run['stdout'].strip(): return 'no_output'
    return 'converted'


def production_sources(revision):
    paths=['crates/chematic-mol/src/mol3000.rs','crates/chematic-mol/src/mol3000/atom_scalars.rs']
    present=set(subprocess.check_output(['git','ls-tree','-r','--name-only',revision,'--',*paths],cwd=ROOT,text=True).splitlines()) if revision else set(paths)
    return {path:(identity(ROOT/path,subprocess.check_output(['git','show',f'{revision}:{path}'],cwd=ROOT) if revision else None) if path in present else dict(absent=True)) for path in paths}


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--wasm-module',type=Path,required=True)
    p.add_argument('--output',type=Path,required=True)
    p.add_argument('--baseline',action='store_true')
    p.add_argument('--runtime-revision')
    args=p.parse_args()
    import chematic
    from rdkit import rdBase
    from rdkit.Chem import rdmolfiles
    extension=importlib.import_module('chematic.chematic')
    path=ROOT/'validation/v3000_atom_scalars/cases.json';packet=strict_json(path.read_text())
    command=['node','-e',NODE,str(args.wasm_module.resolve())]
    run=terminal(bounded_run(command,json.dumps(packet).encode()));run['command']=command
    require(run['resource_status']=='completed' and run['exit_code']==0,'Node probe failed')
    wasm=strict_json(run['stdout'])
    require(len(wasm['rows'])==len(packet['cases']) and len(wasm['invalid'])==len(packet['invalid']),'Node row count mismatch')
    rows=[]
    for c,js in zip(packet['cases'],wasm['rows']):
        require(c['id']==js['id'],'Node order mismatch')
        mol=chematic.from_mol_v3000(c['input'])
        py=dict(id=c['id'],smiles=mol.smiles,moljson=mol.to_moljson(),common=chematic.convert_format(c['input'],'v3000','smiles'),source=chematic.convert_format(c['input'],'v3000','v3000'))
        row=dict(id=c['id'],scope=c['scope'],input=c['input'],expected=c['expected'],python=py,wasm=js)
        for name,value in [('python',py),('wasm',js)]:
            row[name+'_observation']=dict(direct=observe(c['expected'],value['smiles']),common=observe(c['expected'],value['common']),moljson_fields_exact=moljson_fields_exact(c['expected'],value['moljson']),source_exact=value['source']==c['input'])
        row['openbabel']={}
        for name,flags in [('default',[]),('atom_classes',['-xa'])]:
            converted=obabel(c['input'],flags)
            row['openbabel'][name]=dict(run=converted,observation=observe(c['expected'],converted['stdout'].strip()) if converted['outcome']=='converted' else dict(readable=False,fields_exact=False))
        rows.append(row)
    invalid=[]
    for c,js in zip(packet['invalid'],wasm['invalid']):
        require(c['id']==js['id'],'Node refusal order mismatch')
        py=dict(direct=capture(lambda:chematic.from_mol_v3000(c['input']).smiles),common=capture(lambda:chematic.convert_format(c['input'],'v3000','smiles')),source=capture(lambda:chematic.convert_format(c['input'],'v3000','v3000')))
        invalid.append(dict(id=c['id'],input=c['input'],field=c['field'],python=py,wasm=js,openbabel=obabel(c['input'],['-xa'])))
    summary={}
    for scope in ['ordinary','quoted_compatibility']:
        group=[r for r in rows if r['scope']==scope]
        summary[scope]=dict(cases=len(group),python=sum(r['python_observation']['direct']['fields_exact'] and r['python_observation']['common']['fields_exact'] for r in group),wasm=sum(r['wasm_observation']['direct']['fields_exact'] and r['wasm_observation']['common']['fields_exact'] for r in group),openbabel_default=sum(r['openbabel']['default']['observation']['fields_exact'] for r in group),openbabel_atom_classes=sum(r['openbabel']['atom_classes']['observation']['fields_exact'] for r in group))
    summary['invalid']=dict(cases=len(invalid),python_field_errors=sum(all(not r['python'][key]['accepted'] and r['field'] in r['python'][key]['error'] for key in ['direct','common','source']) for r in invalid),wasm_field_errors=sum(all(not r['wasm'][key]['accepted'] and r['field'] in r['wasm'][key]['error'] for key in ['direct','common','source']) for r in invalid),openbabel_outcomes={name:sum(r['openbabel']['outcome']==name for r in invalid) for name in ['converted','error','no_output','signal','resource_limit']})
    result=dict(schema='chematic.v3000-atom-scalar-observation.v1',baseline=args.baseline,head=subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip(),runtime_revision=args.runtime_revision or 'working_tree',production_sources=production_sources(args.runtime_revision),packet=identity(path),harness=identity(__file__),python_package=identity(chematic.__file__),python_extension=identity(extension.__file__),wasm_js=identity(args.wasm_module),wasm_binary=identity(args.wasm_module.with_name('chematic_wasm_bg.wasm')),rdkit_version=rdBase.rdkitVersion,rdkit_parser=identity(rdmolfiles.__file__),rdkit_options=dict(sanitize=False,removeHs=False),openbabel=identity('/opt/homebrew/bin/obabel'),openbabel_version=terminal(bounded_run(['/opt/homebrew/bin/obabel','-V'])),node_run=run,summary=summary,cases=rows,invalid=invalid)
    args.output.parent.mkdir(parents=True,exist_ok=True);args.output.write_bytes(gzip.compress((json.dumps(result,indent=2)+'\n').encode(),mtime=0))
    print(json.dumps(summary))
    if not args.baseline:
        require(all(r[e+'_observation']['direct']['fields_exact'] and r[e+'_observation']['common']['fields_exact'] and r[e+'_observation']['moljson_fields_exact'] and r[e+'_observation']['source_exact'] for r in rows for e in ['python','wasm']),'positive field/source contract failed')
        require(summary['invalid']['python_field_errors']==len(invalid) and summary['invalid']['wasm_field_errors']==len(invalid),'unsafe scalar acceptance')


if __name__=='__main__':main()
