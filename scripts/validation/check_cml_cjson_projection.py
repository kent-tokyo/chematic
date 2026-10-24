#!/usr/bin/env python3
"""Bounded paired CML->CJSON graph/XYZ/name observations; no general toolkit claim."""
from __future__ import annotations
import argparse
import gzip
import importlib
import json
import hashlib
import platform
from pathlib import Path
import struct
import subprocess
import sys
import xml.etree.ElementTree as ET
ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT))
from scripts.checks.check_interchange_preservation import bounded_run, identity, require, strict_json

NODE="""const fs=require('node:fs'),w=require(process.argv[1]);
process.stdout.write(JSON.stringify(JSON.parse(fs.readFileSync(0,'utf8')).map(c=>{
try{return {output:w.convert_common_format(c.input,'cml','cjson')}}catch(e){return {error:String(e)}}})));"""
NUMBERS={'H':1,'C':6,'N':7,'O':8,'Cl':17}
ORDERS={'1':1,'S':1,'single':1,'2':2,'D':2,'double':2,'3':3,'T':3,'triple':3}
def source_fields(text):
    root=ET.fromstring(text)
    mol=root if root.tag.split('}')[-1]=='molecule' else list(root)[0]
    ns=mol.tag[:-8]
    atoms=[];bonds=[]
    for child in mol:
        if child.tag==ns+'atom':atoms.append(child)
        if child.tag==ns+'bond':bonds.append(child)
        if child.tag==ns+'atomArray':atoms.extend(child)
        if child.tag==ns+'bondArray':bonds.extend(child)
    ids={a.attrib['id']:i for i,a in enumerate(atoms)}
    names=[''.join(n.itertext()) for n in mol if n.tag==ns+'name']
    if 'title' in mol.attrib:names.append(mol.attrib['title'])
    require(len(set(names))<=1,'positive fixture has conflicting names')
    xyz=None if not any('x3' in a.attrib for a in atoms) else [[float(a.attrib[k]) for k in ('x3','y3','z3')] for a in atoms]
    return dict(numbers=[NUMBERS[a.attrib['elementType']] for a in atoms],
        charges=[int(a.attrib.get('formalCharge',0)) for a in atoms],
        isotopes=[int(a.attrib.get('isotopeNumber',a.attrib.get('isotope',0))) for a in atoms],
        bonds=sorted([*sorted(ids[r] for r in b.attrib['atomRefs2'].split()),ORDERS[b.attrib['order']]] for b in bonds),
        coordinates=xyz,name=names[0] if names else None)

def bits(values):return [struct.pack('>d',float(v)).hex() for v in values]
def observe(expected,output):
    fields=dict(graph=False,isotopes=False,coordinates=False,name=False)
    try:
        d=strict_json(output);a=d['atoms'];n=len(expected['numbers']);b=d.get('bonds',{})
        numbers=a['elements']['number'];charges=a.get('formalCharges',[0]*n);isotopes=a.get('isotopes',[0]*n)
        require(all(isinstance(v,list) and all(type(x) is int for x in v) for v in [numbers,charges,isotopes]),'invalid atom arrays')
        indices=b.get('connections',{}).get('index',[]);orders=b.get('order',[])
        require(type(indices) is list and all(type(i) is int for i in indices),'invalid indices')
        require(type(orders) is list and all(type(o) in (int,float) for o in orders),'invalid orders')
        require(len(indices)==2*len(orders),'mismatched candidate bonds')
        bonds=sorted([*sorted(indices[i*2:i*2+2]),o] for i,o in enumerate(orders))
        fields['graph']=(a['elements']['number']==expected['numbers'] and a.get('formalCharges',[0]*n)==expected['charges'] and bonds==expected['bonds'])
        fields['isotopes']=a.get('isotopes',[0]*n)==expected['isotopes']
        xyz=a.get('coords',{}).get('3d');want=expected['coordinates']
        require(xyz is None or type(xyz) is list and all(type(v) in (int,float) for v in xyz),'invalid coordinates')
        fields['coordinates']=(xyz is None if want is None else xyz is not None and bits(xyz)==bits(v for p in want for v in p))
        names=[d[k] for k in ['name'] if k in d]
        if 'name' in d.get('properties',{}):names.append(d['properties']['name'])
        fields['name']=(not names if expected['name'] is None else bool(names) and all(v==expected['name'] for v in names))
        return dict(json_output=True,fields=fields,full=all(fields.values()))
    except (ValueError,KeyError,TypeError,AttributeError,OverflowError,struct.error):
        return dict(json_output=False,fields=fields,full=False)

def terminal(result):
    return {**result,**{k:result[k].decode('utf8',errors='replace') for k in ('stdout','stderr')}}
def inputs(packet):
    out=[]
    for group in ('valid','invalid'):
        for c in packet[group]:
            s=c['input'];lf=s.replace('><','>\n<')
            for variant,text in [('compact',s),('lf',lf),('crlf',lf.replace('\n','\r\n'))]:
                out.append(dict(id=c['id']+'/'+variant,group=group,input=text,expected=source_fields(text) if group=='valid' else None))
    return out

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--wasm-module',type=Path,required=True);p.add_argument('--cli',type=Path,required=True)
    p.add_argument('--obabel',type=Path,default=Path('/opt/homebrew/bin/obabel'))
    p.add_argument('--output',type=Path,required=True);p.add_argument('--baseline',action='store_true');p.add_argument('--runtime-revision');args=p.parse_args()
    import chematic
    ext=importlib.import_module('chematic.chematic')
    fixture=ROOT/'validation/cml_cjson_projection/cases.json';cases=inputs(strict_json(fixture.read_text()))
    node_input=json.dumps(cases).encode()
    raw=bounded_run(['node','-e',NODE,str(args.wasm_module.resolve())],node_input)
    require(raw['resource_status']=='completed' and raw['exit_code']==0,'Node batch failed')
    nodes=strict_json(raw['stdout']);require(len(nodes)==len(cases),'Node count mismatch')
    labels=['python','wasm','cli','openbabel','openbabel_single']
    for c,w in zip(cases,nodes,strict=True):
        c['wasm']=w
        try:c['python']={'output':chematic.convert_format(c['input'],'cml','cjson')}
        except ValueError as e:c['python']={'error':str(e)}
        for label,cmd in [('cli',[str(args.cli.resolve()),'convert','--input-format','cml','--output-format','cjson']),('openbabel',[str(args.obabel.resolve()),'-icml','-ocjson']),('openbabel_single',[str(args.obabel.resolve()),'-icml','-ocjson','-xs'])]:
            r=terminal(bounded_run(cmd,c['input'].encode()));r['command']=cmd
            if r['resource_status']=='completed' and r['exit_code']==0:r['output']=r['stdout']
            c[label]=r
        for label in labels:
            r=c[label];r['emitted']=bool(r.get('output','').strip())
            if c['group']=='valid':r['observation']=observe(c['expected'],r.get('output',''))
    summary={}
    for label in labels:
        pos=[c[label] for c in cases if c['group']=='valid'];neg=[c[label] for c in cases if c['group']=='invalid']
        summary[label]=dict(valid_count=len(pos),valid_full=sum(r['observation']['full'] for r in pos),valid_no_output=sum(not r['emitted'] for r in pos),
            fields={k:sum(r['observation']['fields'][k] for r in pos) for k in ['graph','isotopes','coordinates','name']},
            declared_name_count=sum(c['expected']['name'] is not None for c in cases if c['group']=='valid'),
            declared_name_retained=sum(c['expected']['name'] is not None and c[label]['observation']['fields']['name'] for c in cases if c['group']=='valid'),
            outside_profile_count=len(neg),outside_profile_no_output=sum(not r['emitted'] for r in neg),
            valid_process_signal_count=sum(r.get('exit_code',0)<0 for r in pos),
            outside_process_signal_count=sum(r.get('exit_code',0)<0 for r in neg))
    paths=['crates/chematic-mol/src/conversion.rs','crates/chematic-mol/src/lib.rs','crates/chematic-mol/src/cml.rs','crates/chematic-core/src/valence.rs','crates/chematic-py/src/formats.rs','crates/chematic-mol/src/conversion/cml_cjson.rs','crates/chematic-mol/src/conversion/cml_source.rs','crates/chematic-mol/src/cml/geometry.rs','crates/chematic-mol/src/cml/geometry_xml.rs','crates/chematic-mol/src/cjson.rs','crates/chematic-mol/src/cjson_record.rs','crates/chematic-py/src/formats/interchange.rs','crates/chematic-py/python/chematic/__init__.py','crates/chematic-wasm/src/mol_io/common_formats.rs','crates/chematic-cli/src/conversion_format.rs','crates/chematic-cli/src/main.rs']
    record=dict(schema='chematic.cml-cjson-projection.v1',baseline=args.baseline,
        scope='Original synthetic selected graph/XYZ/name fields only; refusals are not preservation wins; no independent holdout, chemistry/stereo/schema or speed claim',
        environment=dict(platform=platform.platform(),python=sys.version,node=terminal(bounded_run(['node','--version']))),
        node_run={**terminal(raw),'command':['node','-e','<NODE source in hashed runner>',str(args.wasm_module.resolve())], 'input_sha256':hashlib.sha256(node_input).hexdigest()},
        runtime_revision=args.runtime_revision,head=subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip(),
        source_files=[] if args.baseline else [identity(ROOT/path) for path in paths],
        runner=identity(__file__),fixture=identity(fixture),python_package=identity(chematic.__file__),python_extension=identity(ext.__file__),
        wasm_js=identity(args.wasm_module),wasm_binary=identity(args.wasm_module.with_name('chematic_wasm_bg.wasm')),cli=identity(args.cli),
        openbabel=identity(args.obabel),openbabel_version=terminal(bounded_run([str(args.obabel),'-V'])),summary=summary,cases=cases)
    data=(json.dumps(record,ensure_ascii=False,indent=2)+'\n').encode()
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_bytes(gzip.compress(data,mtime=0) if args.output.suffix=='.gz' else data)
    print(json.dumps(summary))
    if not args.baseline:
        for label in ['python','wasm','cli']:
            require(all(('error' in c[label] if label!='cli' else c[label]['resource_status']=='completed' and c[label]['exit_code']>0)
                        for c in cases if c['group']=='invalid'),'candidate refusal lacks explicit error')
            require(summary[label]['valid_full']==33 and summary[label]['outside_profile_no_output']==75,'candidate failed')
if __name__=='__main__':main()
