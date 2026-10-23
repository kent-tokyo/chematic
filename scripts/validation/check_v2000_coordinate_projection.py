#!/usr/bin/env python3
"""Bounded source-coordinate comparison for V2000 -> CJSON, not full-format parity."""
from __future__ import annotations

import argparse
import hashlib
import importlib
import json
from pathlib import Path
import platform
import struct
import subprocess
import sys

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
from scripts.checks.check_interchange_preservation import bounded_run, identity, require, strict_json
from scripts.validation.check_v3000_coordinate_projection import observe, terminal

NODE = """
const fs=require('node:fs'),w=require(process.argv[1]);
process.stdout.write(JSON.stringify(JSON.parse(fs.readFileSync(0,'utf8')).map(c=>{
try{return {id:c.id,output:w.convert_common_format(c.input,c.format,'cjson')}}
catch(e){return {id:c.id,error:String(e)}}})));
"""


def cases(packet):
    result=[]
    for c in packet['cases']:
        # Decode fixed columns independently of chematic's readers.
        rows=c['input'].splitlines()[4:4+len(c['expected']['numbers'])]
        decoded=[[float(row[a:a+10]) for a in (0,10,20)] for row in rows]
        bits=lambda ps:[struct.pack('>d',v) for row in ps for v in row]
        require(bits(decoded)==bits(c['expected']['coordinates']),'source-column fixture mismatch')
        for variant,text in [('lf',c['input']),('crlf',c['input'].replace('\n','\r\n')),('no-final-newline',c['input'].rstrip('\n'))]:
            for fmt in ['mol','sdf']:
                source=text+'\n$$$$\n' if fmt=='sdf' else text
                result.append(dict(id=f"{c['id']}/{variant}/{fmt}",input=source,format=fmt,expected=c['expected']))
    return result


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--wasm-module',type=Path,required=True)
    parser.add_argument('--obabel',type=Path,default=Path('/opt/homebrew/bin/obabel'))
    parser.add_argument('--output',type=Path,required=True)
    parser.add_argument('--baseline',action='store_true')
    args=parser.parse_args()
    import chematic
    native=importlib.import_module('chematic.chematic')
    fixture=ROOT/'validation/v2000_coordinate_projection/cases.json'
    inputs=cases(strict_json(fixture.read_text()))
    raw=bounded_run(['node','-e',NODE,str(args.wasm_module.resolve())],json.dumps(inputs).encode())
    require(raw['resource_status']=='completed' and raw['exit_code']==0,'Node runner failed')
    node=strict_json(raw['stdout']);require(len(node)==len(inputs),'Node count mismatch')
    results=[]
    for c,w in zip(inputs,node,strict=True):
        require(c['id']==w['id'],'Node order mismatch')
        row=dict(id=c['id'],input_sha256=hashlib.sha256(c['input'].encode()).hexdigest(),expected=c['expected'],wasm=w)
        try: row['python']=dict(output=chematic.convert_format(c['input'],c['format'],'cjson'))
        except ValueError as e: row['python']=dict(error=str(e))
        for label,options in [('openbabel',[]),('openbabel_single',['-xs'])]:
            command=[str(args.obabel.resolve()),'-i'+c['format'],'-ocjson']+options
            value=terminal(bounded_run(command,c['input'].encode()))
            value['command']=command
            if value['resource_status']=='completed' and value['exit_code']==0:value['output']=value['stdout']
            row[label]=value
        for label in ('python','wasm','openbabel','openbabel_single'):
            row[label]['observation']=observe(c['expected'],row[label].get('output'))
        results.append(row)
    labels=('python','wasm','openbabel','openbabel_single')
    summary={label:{metric:sum(r[label]['observation'][metric] for r in results) for metric in results[0][label]['observation']} for label in labels}
    sources=['crates/chematic-mol/src/mol2000/coordinates.rs','crates/chematic-mol/src/mol2000/reader.rs',
             'crates/chematic-py/src/formats/v2000_coordinates.rs','crates/chematic-py/python/chematic/__init__.py',
             'crates/chematic-wasm/src/mol_io/common_formats.rs',
             'scripts/validation/check_v3000_coordinate_projection.py',
             'scripts/validation/check_v2000_coordinate_projection.py']
    record=dict(schema='chematic.v2000-coordinate-projection-results.v1',baseline=args.baseline,
        scope='source xyz bits and indexed elements/formal charges/isotopes/bond orders; no stereo/metadata/performance claim',
        head=subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip(),
        source_state='dirty source candidate' if not args.baseline else 'prior local package; see artifact hashes',
        host=platform.platform(),python=platform.python_version(),packet=identity(fixture),
        source_files=[] if args.baseline else [identity(ROOT/p) for p in sources],
        python_package=identity(chematic.__file__),python_extension=identity(native.__file__),
        wasm_module=identity(args.wasm_module),wasm_binary=identity(args.wasm_module.with_name('chematic_wasm_bg.wasm')),
        openbabel=identity(args.obabel),openbabel_version=terminal(bounded_run([str(args.obabel.resolve()),'-V'])),
        case_count=len(results),summary=summary,cases=results)
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(record,indent=2)+'\n')
    print(json.dumps(summary))
    if not args.baseline:
        require(all(all(r[label]['observation'].values()) for r in results for label in labels[:2]),'candidate projection failure')


if __name__=='__main__':main()
