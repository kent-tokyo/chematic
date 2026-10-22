#!/usr/bin/env python3
"""Bounded development comparison of V3000 -> CJSON coordinates and indexed graph."""
from __future__ import annotations
import argparse
import importlib
import json
import math
from pathlib import Path
import struct
import sys
ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
from scripts.check_interchange_preservation import bounded_run, identity, require, strict_json

NODE = """
const fs=require('node:fs'), w=require(process.argv[1]);
const cases=JSON.parse(fs.readFileSync(0,'utf8'));
process.stdout.write(JSON.stringify(cases.map(c=>{try{return {id:c.id,output:w.convert_common_format(c.input,'v3000','cjson')}}catch(e){return {id:c.id,error:String(e)}}})));
"""

def bits(values):
    return [struct.pack('>d', float(v)).hex() for v in values]

def observe(expected, output):
    result = dict(json_object=False, indexed_graph=False, coordinates_present=False, coordinate_bits=False, isotopes=False)
    if output is None: return result
    try:
        value = strict_json(output)
        require(isinstance(value, dict), 'not an object')
        result['json_object'] = True
        atoms = value['atoms']; numbers = atoms['elements']['number']; n = len(numbers)
        require(isinstance(numbers, list) and all(type(v) is int for v in numbers), 'invalid atomic numbers')
        correspondence = numbers == expected['numbers']
        coords = atoms.get('coords', {}).get('3d')
        result['coordinates_present'] = isinstance(coords, list) and len(coords) == 3*n
        finite = result['coordinates_present'] and all(type(v) in (int, float) and math.isfinite(v) for v in coords)
        result['coordinate_bits'] = correspondence and finite and bits(coords) == bits(v for row in expected['coordinates'] for v in row)
        isotopes = atoms.get('isotopes', [0]*n)
        result['isotopes'] = correspondence and isinstance(isotopes, list) and all(type(v) is int for v in isotopes) and isotopes == expected['isotopes']
        bonds = value.get('bonds', {}); indices = bonds.get('connections', {}).get('index', []); orders = bonds.get('order', [])
        require(isinstance(indices, list) and all(type(v) is int and 0 <= v < n for v in indices), 'invalid bond indices')
        require(isinstance(orders, list) and all(type(v) in (int, float) and math.isfinite(v) for v in orders), 'invalid bond orders')
        require(len(indices) == 2 * len(orders), 'bond count mismatch')
        observed_bonds = sorted([min(indices[2*i:2*i+2]), max(indices[2*i:2*i+2]), order] for i, order in enumerate(orders))
        expected_bonds = sorted([min(a,b), max(a,b), order] for a,b,order in expected['bonds'])
        charges = atoms.get('formalCharges', [0]*n)
        result['indexed_graph'] = correspondence and isinstance(charges, list) and all(type(v) is int for v in charges) and charges == expected['formal_charges'] and observed_bonds == expected_bonds
    except (ValueError, KeyError, TypeError, IndexError, OverflowError): pass
    return result

def terminal(raw):
    return {k:v.decode('utf-8', errors='replace') if isinstance(v,bytes) else v for k,v in raw.items()}

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--wasm-module', type=Path, required=True)
    parser.add_argument('--obabel', type=Path, default=Path('/opt/homebrew/bin/obabel'))
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--baseline', action='store_true')
    args=parser.parse_args()
    import chematic
    extension=importlib.import_module('chematic.chematic')
    packet_path=ROOT/'validation/v3000_coordinate_projection/cases.json'
    packet=strict_json(packet_path.read_text());cases=[]
    import hashlib
    for case in packet['cases']:
        require(hashlib.sha256(case['input'].encode()).hexdigest()==case['sha256'],'fixture hash mismatch')
        for label,text in [('lf',case['input']),('crlf',case['input'].replace('\n','\r\n')),('no-final-newline',case['input'].rstrip('\n'))]:
            cases.append(dict(id=case['id']+'/'+label,input=text,expected=case['expected']))
    raw=bounded_run(['node','-e',NODE,str(args.wasm_module.resolve())],json.dumps(cases).encode())
    require(raw['resource_status']=='completed' and raw['exit_code']==0,'Node runner failed')
    node=strict_json(raw['stdout']); require(len(node)==len(cases),'Node count mismatch')
    results=[]
    for case,w in zip(cases,node):
        require(case['id']==w['id'],'Node order mismatch')
        row=dict(id=case['id'],input_sha256=hashlib.sha256(case['input'].encode()).hexdigest(),expected=case['expected'],wasm=w)
        try:row['python']=dict(output=chematic.convert_format(case['input'],'v3000','cjson'))
        except ValueError as error:row['python']=dict(error=str(error))
        for label,options in [('openbabel',[]),('openbabel_single',['-xs'])]:
            command=[str(args.obabel.resolve()),'-imol','-ocjson']+options
            result=terminal(bounded_run(command,case['input'].encode()));result['command']=command
            if result['resource_status']=='completed' and result['exit_code']==0:result['output']=result['stdout']
            row[label]=result
        for label in ['python','wasm','openbabel','openbabel_single']:row[label]['observation']=observe(case['expected'],row[label].get('output'))
        results.append(row)
    labels=['python','wasm','openbabel','openbabel_single']
    summary={label:{metric:sum(r[label]['observation'][metric] for r in results) for metric in results[0][label]['observation']} for label in labels}
    record=dict(schema='chematic.v3000-coordinate-projection-results.v1',baseline=args.baseline,packet=identity(packet_path),python_package=identity(chematic.__file__),python_extension=identity(extension.__file__),wasm_module=identity(args.wasm_module),wasm_binary=identity(args.wasm_module.with_name('chematic_wasm_bg.wasm')),openbabel=identity(args.obabel),openbabel_version=terminal(bounded_run([str(args.obabel.resolve()),'-V'])),case_count=len(cases),summary=summary,cases=results)
    args.output.write_text(json.dumps(record,indent=2)+'\n');print(json.dumps(summary))
    if not args.baseline:require(all(all(r[label]['observation'].values()) for r in results for label in labels[:2]),'candidate projection failure')

if __name__=='__main__':main()
