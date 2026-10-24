#!/usr/bin/env python3
"""Independent ASE consumer observations of explicitly supplied ExtXYZ frame data."""
from __future__ import annotations
import argparse
import gzip
import importlib
from io import StringIO
import json
import hashlib
from pathlib import Path
import struct
import subprocess
import sys
ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT))
from scripts.checks.check_interchange_preservation import bounded_run, identity, require, strict_json
NODE="""const fs=require('node:fs'),w=require(process.argv[1]),m=w.parse_smiles('CO');
const encode=x=>JSON.stringify(x,(_,v)=>Object.is(v,-0)?'__negative_zero__':v).replaceAll('"__negative_zero__"','-0.0');
process.stdout.write(JSON.stringify(JSON.parse(fs.readFileSync(0,'utf8')).map(c=>{try{
return {output:w.to_extxyz_json(m,encode(c.coords),encode({lattice:c.lattice,properties:c.properties,info:c.info}))}
}catch(e){return {error:String(e)}}})));m.free();"""
def bits(values):return [struct.pack('>d',float(x)).hex() for x in values]
def quote(value):return '"'+value.replace('\\','\\\\').replace('"','\\"')+'"'
def source_text(c):
    header=[]
    if c['lattice'] is not None:header.append('Lattice='+quote(' '.join(repr(x) for x in c['lattice'])))
    schema='species:S:1:pos:R:3'+''.join(f':{name}:R:{len(values[0])}' for name,values in c['properties'].items())
    header.append('Properties='+schema)
    header.extend(f'{k}={quote(v)}' for k,v in c['info'].items())
    lines=[str(len(c['elements'])),' '.join(header)]
    for i,(e,p) in enumerate(zip(c['elements'],c['coords'],strict=True)):
        data=[*p,*(v for values in c['properties'].values() for v in values[i])]
        lines.append(e+' '+' '.join(repr(x) for x in data))
    return '\n'.join(lines)+'\n'
def jsonable(value):
    if isinstance(value,dict):return {k:jsonable(v) for k,v in value.items()}
    if hasattr(value,'tolist'):return jsonable(value.tolist())
    if isinstance(value,(list,tuple)):return [jsonable(v) for v in value]
    return value

def ase_snapshot(text):
    from ase.io import read
    a=read(StringIO(text),format='extxyz')
    arrays={k:dict(shape=list(v.shape),kind=v.dtype.kind,values=bits(v.flat) if v.dtype.kind=='f' else jsonable(v)) for k,v in a.arrays.items() if k not in {'numbers','positions'}}
    return dict(numbers=a.numbers.tolist(),positions=bits(a.positions.flat),cell=bits(a.cell.array.flat),pbc=a.pbc.tolist(),arrays=arrays,info=jsonable(a.info))
def raw_geometry(text):
    lines=text.splitlines();n=int(lines[0]);require(n>=0 and len(lines)==n+2,'frame count/extent')
    rows=[line.split() for line in lines[2:]]
    return dict(elements=[r[0] for r in rows],positions=bits(float(v) for r in rows for v in r[1:4]))
def compare(expected,snapshot):return snapshot==expected

def observe(case,output,expected):
    result=dict(emitted=bool(output.strip()),geometry=False,ase_readable=False,ase_full=False)
    try:
        g=raw_geometry(output);result['geometry']=g==dict(elements=case['elements'],positions=bits(v for row in case['coords'] for v in row))
    except (ValueError,IndexError,OverflowError):pass
    try:
        result['snapshot']=ase_snapshot(output);result['ase_readable']=True;result['ase_full']=compare(expected,result['snapshot'])
    except Exception as e:result['ase_error']=f'{type(e).__name__}: {e}'
    return result

def terminal(raw):return {k:v.decode('utf8',errors='replace') if isinstance(v,bytes) else v for k,v in raw.items()}
def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--wasm-module',type=Path,required=True);p.add_argument('--output',type=Path,required=True)
    p.add_argument('--obabel',type=Path,default=Path('/opt/homebrew/bin/obabel'));p.add_argument('--baseline',action='store_true');args=p.parse_args()
    import ase,ase.io.extxyz,numpy,chematic
    require(ase.__version__=='3.26.0','ASE version must be 3.26.0')
    extension=importlib.import_module('chematic.chematic');fixture=ROOT/'validation/extxyz_output_integrity/cases.json';packet=strict_json(fixture.read_text())
    cases=[dict(c,group=group) for group in ['valid','invalid'] for c in packet[group]]
    node_input=json.dumps(cases,ensure_ascii=False).encode();raw=bounded_run(['node','-e',NODE,str(args.wasm_module.resolve())],node_input)
    require(raw['resource_status']=='completed' and raw['exit_code']==0,'Node failed');node=strict_json(raw['stdout']);require(len(node)==len(cases),'Node count')
    labels=['python','wasm','openbabel_xyz','openbabel_xyz_no_bonds']
    for c,w in zip(cases,node,strict=True):
        c['wasm']=w
        try:c['python']=dict(output=chematic.to_extxyz(chematic.from_smiles('CO'),c['coords'],lattice=c['lattice'],properties=c['properties'],info=c['info']))
        except ValueError as e:c['python']=dict(error=str(e))
        if c['group']=='valid':
            c['source']=source_text(c);expected=ase_snapshot(c['source']);c['expected_ase']=expected
            require(expected['positions']==bits(v for r in c['coords'] for v in r),'authored coordinate oracle mismatch')
            require(expected['cell']==bits(c['lattice'] or [0.]*9),'authored cell convention mismatch')
            for name,rows in c['properties'].items():
                require(expected['arrays'][name]['values']==bits(v for row in rows for v in row),'authored property oracle mismatch')
                width=len(rows[0]);require(expected['arrays'][name]['shape']==([len(rows)] if width==1 else [len(rows),width]),'authored property shape mismatch')
            if c['id']=='typed_external_metadata':
                require(expected['info']=={'step':12,'active':True} and expected['pbc']==[True,False,True],'authored typed metadata mismatch')
            elif c['id']=='explicit_empty_string':require(expected['info']=={'label':''},'authored typed empty string mismatch')
            else:require(expected['info']==c['info'],'authored string metadata mismatch')
            for label,opts in [('openbabel_xyz',[]),('openbabel_xyz_no_bonds',['-ab'])]:
                command=[str(args.obabel),'-ixyz','-oxyz']+opts;r=terminal(bounded_run(command,c['source'].encode()));r['command']=command
                if r['resource_status']=='completed' and r['exit_code']==0:r['output']=r['stdout']
                c[label]=r
            for label in labels:c[label]['observation']=observe(c,c[label].get('output',''),expected)
    summary={}
    for label in labels:
        pos=[c[label]['observation'] for c in cases if c['group']=='valid'];summary[label]=dict(valid_count=len(pos),**{k:sum(r[k] for r in pos) for k in ['emitted','geometry','ase_readable','ase_full']})
        if label in ['python','wasm']:
            neg=[c[label] for c in cases if c['group']=='invalid'];summary[label].update(outside_profile_count=len(neg),explicit_refusal=sum('error' in r and not r.get('output','').strip() for r in neg))
    paths=['crates/chematic-mol/src/xyz.rs','crates/chematic-mol/src/xyz/extxyz.rs','crates/chematic-mol/src/xyz/extxyz_validation.rs','crates/chematic-py/src/formats/xyz.rs','crates/chematic-wasm/src/mol_io.rs','crates/chematic-wasm/src/mol_io/extxyz_output.rs']
    record=dict(schema='chematic.extxyz-output-integrity.v1',baseline=args.baseline,
        scope='Original synthetic writer/ASE consumer geometry-cell-atom-data-info workflow; Open Babel is measured through XYZ, not ASE ExtXYZ or its distinct EXYZ format; no chemical graph/stereo or speed claim',
        head=subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip(),source_files=[] if args.baseline else [identity(ROOT/v) for v in paths],
        fixture=identity(fixture),runner=identity(__file__),python_package=identity(chematic.__file__),python_extension=identity(extension.__file__),
        wasm_js=identity(args.wasm_module),wasm_binary=identity(args.wasm_module.with_name('chematic_wasm_bg.wasm')),
        ase_version=ase.__version__,ase_reader=identity(ase.io.extxyz.__file__),numpy_version=numpy.__version__,
        baseline_runtime_revision='fda4dfac4fd17e30955cf97d4b00b3b8593e9730' if args.baseline else None,
        node_command=['node','-e','<NODE source in hashed runner>',str(args.wasm_module.resolve())],node_input_sha256=hashlib.sha256(node_input).hexdigest(),
        openbabel=identity(args.obabel),openbabel_version=terminal(bounded_run([str(args.obabel),'-V'])),node_run=terminal(raw),summary=summary,cases=cases)
    args.output.parent.mkdir(parents=True,exist_ok=True);args.output.write_bytes(gzip.compress((json.dumps(record,ensure_ascii=False,indent=2)+'\n').encode(),mtime=0));print(json.dumps(summary))
    if not args.baseline:
        for label in ['python','wasm']:
            require(summary[label]['ase_full']==len(packet['valid']) and summary[label]['geometry']==len(packet['valid']),'candidate field loss')
            require(summary[label]['explicit_refusal']==len(packet['invalid']),'candidate failed to refuse')
if __name__=='__main__':main()
