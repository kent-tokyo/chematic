#!/usr/bin/env python3
"""Bounded MOL2 -> CJSON selected-field preservation and refusal development packet."""
from __future__ import annotations
import argparse
import hashlib
import importlib
import json
import math
from pathlib import Path
import struct
import sys
ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT))
from scripts.check_interchange_preservation import bounded_run,identity,require,strict_json
NODE="""
const fs=require('node:fs'),w=require(process.argv[1]);
const cases=JSON.parse(fs.readFileSync(0,'utf8'));
process.stdout.write(JSON.stringify(cases.map(c=>{try{return {id:c.id,output:w.convert_common_format(c.input,'mol2','cjson')}}catch(e){return {id:c.id,error:String(e)}}})));
"""
def bits(values):return [struct.pack('>d',float(v)).hex() for v in values]
def numeric(values):return isinstance(values,list) and all(type(v) in (int,float) and math.isfinite(v) for v in values)
def observe(expected,output):
    result=dict(json_object=False,indexed_graph=False,coordinates=False,name=False,partial_values=False,charge_method=False)
    if output is None:return result
    try:
        value=strict_json(output);require(isinstance(value,dict),'not an object');result['json_object']=True
        atoms=value['atoms'];numbers=atoms['elements']['number'];n=len(expected['numbers'])
        correspondence=isinstance(numbers,list) and all(type(v) is int for v in numbers) and numbers==expected['numbers']
        coords=atoms.get('coords',{}).get('3d')
        result['coordinates']=correspondence and numeric(coords) and len(coords)==3*n and bits(coords)==bits(v for row in expected['coordinates'] for v in row)
        result['name']=value.get('properties',{}).get('name')==expected['name']
        sets=value.get('partialCharges',{})
        if expected['method'] is None:
            result['partial_values']=result['charge_method']=sets=={}
        elif isinstance(sets,dict):
            result['charge_method']=list(sets)==[expected['method']]
            if len(sets)==1:
                values=next(iter(sets.values()))
                result['partial_values']=correspondence and numeric(values) and bits(values)==bits(expected['partial_charges'])
        bonds=value.get('bonds',{});indices=bonds.get('connections',{}).get('index',[]);orders=bonds.get('order',[])
        require(isinstance(indices,list) and all(type(v) is int and 0<=v<n for v in indices),'bad bond indices')
        require(numeric(orders) and len(indices)==2*len(orders),'bad bond orders')
        observed=sorted([min(indices[2*i:2*i+2]),max(indices[2*i:2*i+2]),order] for i,order in enumerate(orders))
        wanted=sorted([min(a,b),max(a,b),order] for a,b,order in expected['bonds'])
        charges=atoms.get('formalCharges',[0]*n)
        result['indexed_graph']=correspondence and isinstance(charges,list) and all(type(v) is int for v in charges) and charges==expected['formal_charges'] and observed==wanted
    except (ValueError,KeyError,TypeError,IndexError,OverflowError):pass
    return result

def terminal(raw):return {k:v.decode('utf-8',errors='replace') if isinstance(v,bytes) else v for k,v in raw.items()}
def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--wasm-module',type=Path,required=True);parser.add_argument('--obabel',type=Path,default=Path('/opt/homebrew/bin/obabel'));parser.add_argument('--output',type=Path,required=True);parser.add_argument('--baseline',action='store_true');args=parser.parse_args()
    import chematic
    extension=importlib.import_module('chematic.chematic');packet_path=ROOT/'validation/mol2_cjson_common/cases.json';packet=strict_json(packet_path.read_text(encoding='utf-8'));cases=[]
    for c in packet['cases']:
        require(hashlib.sha256(c['input'].encode()).hexdigest()==c['sha256'],'fixture hash mismatch')
        for label,text in [('lf',c['input']),('crlf',c['input'].replace('\n','\r\n')),('no-final-newline',c['input'].rstrip('\n'))]:cases.append(dict(id=c['id']+'/'+label,input=text,expected=c['expected'],outcome=c['outcome']))
    raw=bounded_run(['node','-e',NODE,str(args.wasm_module.resolve())],json.dumps(cases).encode());require(raw['resource_status']=='completed' and raw['exit_code']==0,'Node failed');node=strict_json(raw['stdout']);require(len(node)==len(cases),'Node count mismatch');results=[]
    for c,w in zip(cases,node):
        require(c['id']==w['id'],'Node order mismatch');row=dict(id=c['id'],expected_outcome=c['outcome'],input_sha256=hashlib.sha256(c['input'].encode()).hexdigest(),wasm=w)
        try:row['python']=dict(output=chematic.convert_format(c['input'],'mol2','cjson'))
        except ValueError as error:row['python']=dict(error=str(error))
        for label,options in [('openbabel',[]),('openbabel_single',['-xs'])]:
            command=[str(args.obabel.resolve()),'-imol2','-ocjson']+options;r=terminal(bounded_run(command,c['input'].encode()));r['command']=command
            if r['resource_status']=='completed' and r['exit_code']==0:r['output']=r['stdout']
            row[label]=r
        for label in ['python','wasm','openbabel','openbabel_single']:
            output=row[label].get('output');row[label]['emitted_nonempty_text']=bool(output)
            if c['outcome']=='preserved':row[label]['observation']=observe(c['expected'],output)
        results.append(row)
    labels=['python','wasm','openbabel','openbabel_single'];summary={}
    for label in labels:
        good=[r for r in results if r['expected_outcome']=='preserved'];metrics={metric:sum(r[label]['observation'][metric] for r in good) for metric in good[0][label]['observation']}
        summary[label]=dict(preserved_input_count=len(good),selected_fields=metrics,refusals={kind:dict(input_count=sum(r['expected_outcome']==kind for r in results),no_output=sum(r['expected_outcome']==kind and not r[label]['emitted_nonempty_text'] for r in results)) for kind in ['refused','invalid']})
    record=dict(schema='chematic.mol2-cjson-common-results.v1',baseline=args.baseline,packet=identity(packet_path),python_package=identity(chematic.__file__),python_extension=identity(extension.__file__),wasm_module=identity(args.wasm_module),wasm_binary=identity(args.wasm_module.with_name('chematic_wasm_bg.wasm')),openbabel=identity(args.obabel),openbabel_version=terminal(bounded_run([str(args.obabel.resolve()),'-V'])),case_count=len(results),summary=summary,cases=results)
    args.output.write_text(json.dumps(record,indent=2)+'\n',encoding='utf-8');print(json.dumps(summary))
    if not args.baseline:
        for r in results:
            for label in labels[:2]:require(all(r[label]['observation'].values()) if r['expected_outcome']=='preserved' else not r[label]['emitted_nonempty_text'],'candidate preservation/refusal failure')
if __name__=='__main__':main()
