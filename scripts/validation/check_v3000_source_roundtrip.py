#!/usr/bin/env python3
"""Bounded same-format V3000 source and coordinate retention development packet."""
from __future__ import annotations
import argparse
import hashlib
import importlib
import json
import math
from pathlib import Path
import struct
import sys
import tempfile
ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT))
from scripts.check_interchange_preservation import bounded_run,identity,require,strict_json

NODE="""
const fs=require('node:fs'),w=require(process.argv[1]);
const cases=JSON.parse(fs.readFileSync(0,'utf8'));
const call=f=>{try{return {output:f()}}catch(e){return {error:String(e)}}};
process.stdout.write(JSON.stringify(cases.map(c=>({id:c.id,
 common:call(()=>w.convert_common_format(c.input,'v3000','v3000')),
 roundtrip:call(()=>w.roundtrip_mol_v3000_block(c.input))}))));
"""


def atom_coordinates(text):
    # Independent limited observer: token-boundary continuations and source atom
    # order. It does not interpret CFG, groups, query semantics or chemistry.
    logical=[];pending=''
    for line in text.splitlines():
        if not line.startswith('M  V30 '):continue
        payload=line[7:]
        pending+=payload
        if pending.endswith('-'):
            pending=pending[:-1];continue
        logical.append(pending.split());pending=''
    require(not pending,'unterminated observed continuation')
    atoms=[];inside=False
    for parts in logical:
        if parts[:2]==['BEGIN','ATOM']:inside=True;continue
        if parts[:2]==['END','ATOM']:inside=False;continue
        if inside:
            require(len(parts)>=6,'short observed atom row')
            p=[float(v) for v in parts[2:5]]
            require(all(math.isfinite(v) for v in p),'nonfinite observed coordinates')
            atoms.append(dict(id=int(parts[0]),element=parts[1],coords=p))
    require(bool(atoms),'no observed atoms')
    return atoms


def observe(source,output):
    if output is None:return dict(atom_table_observed=False,exact_source=False,coordinates_exact=False,source_ids=False)
    before=atom_coordinates(source)
    try: after=atom_coordinates(output)
    except ValueError:return dict(atom_table_observed=False,exact_source=False,coordinates_exact=False,source_ids=False)
    same_order=[a['element'] for a in before]==[a['element'] for a in after]
    bits=lambda rows:[struct.pack('>d',v) for a in rows for v in a['coords']]
    return dict(atom_table_observed=True,exact_source=source==output,coordinates_exact=same_order and bits(before)==bits(after),
                source_ids=same_order and [a['id'] for a in before]==[a['id'] for a in after])


def terminal(raw):
    return {k:v.decode('utf-8',errors='replace') if isinstance(v,bytes) else v for k,v in raw.items()}


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--wasm-module',required=True,type=Path)
    parser.add_argument('--output',required=True,type=Path)
    parser.add_argument('--obabel',default=Path('/opt/homebrew/bin/obabel'),type=Path)
    parser.add_argument('--baseline',action='store_true')
    args=parser.parse_args()
    import chematic
    extension=importlib.import_module('chematic.chematic')
    packet_path=ROOT/'validation/v3000_source_roundtrip/cases.json'
    packet=strict_json(packet_path.read_text());cases=[]
    for c in packet['cases']:
        require(hashlib.sha256(c['input'].encode()).hexdigest()==c['sha256'],'fixture hash mismatch')
        for label,text in [('lf',c['input']),('crlf',c['input'].replace('\n','\r\n')),('no-final-newline',c['input'].rstrip('\n'))]:
            cases.append(dict(id=c['id']+'/'+label,input=text))
    run=bounded_run(['node','-e',NODE,str(args.wasm_module.resolve())],json.dumps(cases).encode())
    require(run['resource_status']=='completed' and run['exit_code']==0,'Node failed')
    wasm=strict_json(run['stdout']);require(len(wasm)==len(cases),'Node count mismatch')
    results=[]
    for c,w in zip(cases,wasm):
        require(w['id']==c['id'],'Node order mismatch')
        row=dict(id=c['id'],input_sha256=hashlib.sha256(c['input'].encode()).hexdigest())
        try: row['python']=dict(output=chematic.convert_format(c['input'],'v3000','v3000'))
        except ValueError as error:row['python']=dict(error=str(error))
        row['wasm_common']=w['common'];row['wasm_roundtrip']=w['roundtrip']
        for label,options in [('openbabel',['-x3']),('openbabel_atomclass',['-x3','-xa'])]:
            command=[str(args.obabel.resolve()),'-imol','-omol']+options
            raw=terminal(bounded_run(command,c['input'].encode()));raw['command']=command
            if raw['resource_status']=='completed' and raw['exit_code']==0:raw['output']=raw['stdout']
            row[label]=raw
        # Raw-copy reference uses a seekable, explicitly SDF-framed carrier.
        # Title-only reading avoids claiming chemistry validation on this route.
        newline='\r\n' if '\r\n' in c['input'] else '\n'
        framed=c['input']+('' if c['input'].endswith('\n') else newline)+'$$$$'+newline
        with tempfile.TemporaryDirectory(prefix='v3000-copy-') as directory:
            infile=Path(directory)/'source.sdf';infile.write_bytes(framed.encode())
            command=[str(args.obabel.resolve()),'-isdf',str(infile),'-ocopy','-aT']
            raw=terminal(bounded_run(command));raw['command']=command
            raw['input_sha256']=hashlib.sha256(framed.encode()).hexdigest()
            raw['workflow']='SDF-framed seekable raw-copy control, title-only read'
            if raw['resource_status']=='completed' and raw['exit_code']==0:raw['output']=raw['stdout']
            raw['observation']=observe(framed,raw.get('output'))
            row['openbabel_copy_sdf']=raw
        for label in ['python','wasm_common','wasm_roundtrip','openbabel','openbabel_atomclass']:
            row[label]['observation']=observe(c['input'],row[label].get('output'))
        results.append(row)
    labels=['python','wasm_common','wasm_roundtrip','openbabel','openbabel_atomclass','openbabel_copy_sdf']
    summary={e:{metric:sum(r[e]['observation'][metric] for r in results) for metric in ['atom_table_observed','exact_source','coordinates_exact','source_ids']} for e in labels}
    record=dict(schema='chematic.v3000-source-roundtrip-validation.v1',baseline=args.baseline,
        packet=identity(packet_path),python_package=identity(chematic.__file__),python_extension=identity(extension.__file__),
        wasm_module=identity(args.wasm_module),wasm_binary=identity(args.wasm_module.with_name('chematic_wasm_bg.wasm')),
        openbabel=identity(args.obabel),openbabel_version=terminal(bounded_run([str(args.obabel.resolve()),'-V'])),
        case_count=len(results),summary=summary,cases=results)
    args.output.parent.mkdir(parents=True,exist_ok=True);args.output.write_text(json.dumps(record,indent=2)+'\n')
    print(json.dumps(summary))
    if not args.baseline:
        require(all(all(r[e]['observation'].values()) for r in results for e in labels[:3]),'source retention failed')


if __name__=='__main__':main()
