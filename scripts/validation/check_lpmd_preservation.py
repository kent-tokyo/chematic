#!/usr/bin/env python3
"""Observe all-frame LPMD preservation through real Python, Node/WASM and CLI.

The independent observer checks record framing, all source fields and a separate
Cartesian projection. Open Babel default self-roundtrip and source-to-XYZ reads
are separate observations; successful exit alone never establishes preservation.
"""
from __future__ import annotations
import argparse
import hashlib
import importlib
import json
import math
from pathlib import Path
import subprocess
import sys
import tempfile

ROOT = Path(__file__).resolve().parents[2]
PACKET = ROOT / 'validation/lpmd_preservation/cases.json'

def digest(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as f:
        for block in iter(lambda: f.read(1 << 20), b''): h.update(block)
    return h.hexdigest()

def command(args, text=None):
    return subprocess.run(list(map(str,args)), input=text, text=True, encoding='utf-8', capture_output=True, timeout=30, check=False)

def physical(record):
    """Independent row-vector convention; no call into the candidate parser."""
    idx = [record['columns'].index(c) for c in ('X','Y','Z')]
    output = []
    for frame in record['frames']:
        cell = [[float(x) for x in v] for v in frame['cell']]
        atoms = []
        for row in frame['atoms']:
            fraction = [float(row[i]) for i in idx]
            position = []
            for j in range(3):
                terms = [fraction[i]*cell[i][j] for i in range(3)]
                if any(not math.isfinite(v) or (v == 0 and fraction[i] != 0 and cell[i][j] != 0) for i,v in enumerate(terms)):
                    raise ValueError('unsafe Cartesian product')
                result = terms[0]+terms[1]+terms[2]
                if not math.isfinite(result): raise ValueError('unsafe Cartesian sum')
                position.append(result)
            atoms.append(position)
        output.append(atoms)
    return output

def bits(value):
    if isinstance(value,list): return [bits(v) for v in value]
    return float(value).hex()

def read_layout(text):
    """Framing inspector, independent of chematic; refuses a missing SYM column."""
    lines = text.split('\n')
    while lines and not lines[-1].strip(): lines.pop()
    if len(lines)<2 or lines[0].split()!=['LPMD','2.0','L']: raise ValueError('not LPMD 2.0 L')
    header = lines[1].split()
    if not header or header[0]!='HDR': raise ValueError('missing HDR')
    columns = header[1:]
    if len(columns)!=len(set(columns)) or any(c not in columns for c in ('SYM','X','Y','Z')): raise ValueError('missing or duplicate required columns')
    frames=[]; i=2
    while i<len(lines):
        if not lines[i].strip().isdigit(): raise ValueError('invalid frame count')
        count=int(lines[i]);i+=1
        if not count or i>=len(lines): raise ValueError('empty or truncated frame')
        values=lines[i].split();i+=1
        if len(values)!=9: raise ValueError('invalid cell width')
        rows=[line.split() for line in lines[i:i+count]];i+=count
        if len(rows)!=count or any(len(row)!=len(columns) for row in rows): raise ValueError('atom row width differs from HDR')
        frames.append({'cell':[values[n:n+3] for n in (0,3,6)],'atoms':rows})
    if not frames: raise ValueError('no observable frames')
    return {'columns':columns,'frames':frames}

def xyz(text):
    lines=text.split('\n');i=0;frames=[]
    while i<len(lines) and lines[i].strip():
        count=int(lines[i]);i+=2
        rows=[line.split() for line in lines[i:i+count]];i+=count
        if not count or len(rows)!=count or any(len(row)!=4 for row in rows): raise ValueError('truncated XYZ')
        positions=[[float(v) for v in r[1:]] for r in rows]
        if any(not math.isfinite(v) for row in positions for v in row): raise ValueError('nonfinite XYZ')
        frames.append({'elements':[r[0] for r in rows],'positions':positions})
    if not frames: raise ValueError('empty XYZ')
    if any(line.strip() for line in lines[i:]): raise ValueError('trailing XYZ data')
    return frames

def geometry(record, output):
    try: frames=xyz(output)
    except (ValueError,IndexError) as e: return {'observable':False,'error':str(e)}
    idx=record['columns'].index('SYM')
    elements=[[r[idx] for r in f['atoms']] for f in record['frames']]
    expected=physical(record)
    shape=[len(f['positions']) for f in frames]==[len(f) for f in expected]
    actual=[f['positions'] for f in frames]
    return {'observable':True,'frame_and_atom_counts':shape,'elements':elements==[f['elements'] for f in frames],
            'cartesian_bits':shape and bits(actual)==bits(expected),
            'max_abs_coordinate_error_angstrom':max((abs(a-b) for frame_a,frame_b in zip(actual,expected) for row_a,row_b in zip(frame_a,frame_b) for a,b in zip(row_a,row_b)),default=0) if shape else None}

def checkpoint(module, case):
    record=module.parse_lpmd(case['text']); states=[]; positions=[]
    for _ in range(4):
        states.append(record);positions.append(module.lpmd_positions_angstrom(record))
        record=module.parse_lpmd(module.write_lpmd(record))
    return {'states':states,'positions':positions}

NODE = r'''
const fs=require('node:fs');const w=require(process.argv[1]);
const p=JSON.parse(fs.readFileSync(0,'utf8'));const result={valid:[],refused:[],unsafe:[]};
for(const c of p.valid){let r=w.parse_lpmd_json(c.text);let states=[],positions=[];for(let i=0;i<4;i++){states.push(JSON.parse(r));positions.push(w.lpmd_positions_angstrom_json(r));r=w.parse_lpmd_json(w.write_lpmd_from_json(r));}result.valid.push({states,positions});}
for(const c of p.refused){try{w.parse_lpmd_json(c.text);result.refused.push(false);}catch(e){result.refused.push(String(e).startsWith('LPMD line '));}}
for(const c of p.projection_refused){let r=w.parse_lpmd_json(c.text);let raw=JSON.parse(w.parse_lpmd_json(w.write_lpmd_from_json(r)));let refused=false;try{w.lpmd_positions_angstrom_json(r);}catch(e){refused=String(e).startsWith('LPMD line ');}result.unsafe.push({raw,refused});}
process.stdout.write(JSON.stringify(result));
'''

def observe(args):
    packet=json.loads(PACKET.read_text(encoding='utf-8'))
    sys.path.insert(0,str(args.python_path.resolve()));module=importlib.import_module('chematic')
    py={'valid':[checkpoint(module,c) for c in packet['valid']], 'refused':[], 'unsafe':[]}
    for c in packet['refused']:
        try: module.parse_lpmd(c['text']);py['refused'].append(False)
        except ValueError as e: py['refused'].append(str(e).startswith('LPMD line '))
    for c in packet['projection_refused']:
        r=module.parse_lpmd(c['text']);entry={'raw':module.parse_lpmd(module.write_lpmd(r)),'refused':False}
        try: module.lpmd_positions_angstrom(r)
        except ValueError as e: entry['refused']=str(e).startswith('LPMD line ')
        py['unsafe'].append(entry)
    node=command(['node','-e',NODE,args.wasm_module.resolve()],json.dumps(packet))
    if node.returncode: raise RuntimeError(node.stderr)
    wasm=json.loads(node.stdout)
    for entry in wasm['valid']: entry['positions']=[json.loads(s) for s in entry['positions']]
    cli={'valid':[],'refused':[],'unsafe':[]}
    for c in packet['valid']:
        text=c['text'];entry={'states':[],'positions':[]}
        for _ in range(4):
            r=command([args.cli,'lpmd','--inspect'],text);p=command([args.cli,'lpmd','--angstrom'],text)
            if r.returncode or p.returncode: raise RuntimeError(r.stderr+p.stderr)
            entry['states'].append(json.loads(r.stdout));entry['positions'].append(json.loads(p.stdout))
            out=command([args.cli,'lpmd','--from-json'],r.stdout)
            if out.returncode: raise RuntimeError(out.stderr)
            text=out.stdout
        cli['valid'].append(entry)
    for c in packet['refused']:
        out=command([args.cli,'lpmd'],c['text']);cli['refused'].append(out.returncode==2 and not out.stdout and 'LPMD line ' in out.stderr)
    for c in packet['projection_refused']:
        out=command([args.cli,'lpmd'],c['text']);projection=command([args.cli,'lpmd','--angstrom'],c['text'])
        cli['unsafe'].append({'raw':read_layout(out.stdout),'refused':projection.returncode==2 and not projection.stdout})
    candidates={}
    for name,result in [('python',py),('node_wasm',wasm),('cli',cli)]:
        if len(result['valid'])!=len(packet['valid']) or len(result['refused'])!=len(packet['refused']) or len(result['unsafe'])!=len(packet['projection_refused']): raise AssertionError('incomplete runtime response')
        entries=[]
        for c,entry in zip(packet['valid'],result['valid']):
            entries.append({'id':c['id'],'all_fields_all_four_states':len(entry['states'])==4 and all(r==c['record'] for r in entry['states']),
                            'projection_bits_all_four_states':len(entry['positions'])==4 and all(bits(p)==bits(physical(c['record'])) for p in entry['positions'])})
        safe=[entry['raw']==c['record'] and entry['refused'] for c,entry in zip(packet['projection_refused'],result['unsafe'])]
        candidates[name]={'valid':entries,'refused':result['refused'],'raw_retained_unsafe_projection_refused':safe}
    external=[]
    with tempfile.TemporaryDirectory(prefix='lpmd-observer-') as directory:
        d=Path(directory)
        for c in packet['valid']:
            source=d/'in.lpmd';selfout=d/'out.lpmd';source_xyz=d/'source.xyz';back_xyz=d/'back.xyz'
            for f in (selfout,source_xyz,back_xyz): f.unlink(missing_ok=True)
            source.write_text(c['text'],encoding='utf-8')
            receipts=[]
            for sourcefile,fmt,target in [(source,'lpmd',selfout),(source,'xyz',source_xyz),(selfout,'xyz',back_xyz)]:
                out=command([args.openbabel,'-ilpmd',sourcefile,'-o'+fmt,'-O',target,'-ab'])
                receipts.append({'exit':out.returncode,'stderr':out.stderr,'output':target.read_text(encoding='utf-8') if target.exists() else ''})
            try: layout=read_layout(receipts[0]['output']);valid=True;failure=None
            except (ValueError,IndexError) as e: layout=None;valid=False;failure=str(e)
            external.append({'id':c['id'],'default_self_output_valid_declared_schema':valid,'schema_error':failure,
                'raw_record_equal':layout==c['record'],'source_read_xyz':geometry(c['record'],receipts[1]['output']),
                'self_roundtrip_xyz':geometry(c['record'],receipts[2]['output']),'receipts':receipts})
    source_paths=['crates/chematic-mol/src/lpmd.rs','crates/chematic-mol/src/lpmd/numeric.rs','crates/chematic-mol/src/lpmd/reader.rs','crates/chematic-py/src/formats/lpmd.rs','crates/chematic-wasm/src/lpmd.rs','crates/chematic-cli/src/lpmd.rs','scripts/validation/check_lpmd_preservation.py','validation/lpmd_preservation/cases.json']
    result={'schema_version':1,'scope':'local source candidate, bounded LPMD 2.0 L; no release or toolkit-wide superiority claim',
            'candidate_base':command(['git','-C',ROOT,'merge-base','HEAD','origin/main']).stdout.strip(),
            'source_sha256':{p:digest(ROOT/p) for p in source_paths},
            'runtime_sha256':{'cli':digest(args.cli),'wasm_js':digest(args.wasm_module),'wasm':digest(args.wasm_module.with_name('chematic_wasm_bg.wasm')),
                'python_extension':{str(p.name):digest(p) for p in args.python_path.rglob('*.so')}},
            'versions':{'python':sys.version,'node':command(['node','--version']).stdout.strip(),'openbabel':command([args.openbabel,'-V']).stdout.strip()},
            'counts':{k:len(packet[k]) for k in ('valid','refused','projection_refused')},'candidates':candidates,'openbabel':external}
    args.output.parent.mkdir(parents=True,exist_ok=True);args.output.write_text(json.dumps(result,indent=2,ensure_ascii=False)+'\n',encoding='utf-8')
    for name,result in candidates.items():
        assert all(e['all_fields_all_four_states'] and e['projection_bits_all_four_states'] for e in result['valid']),name
        assert all(result['refused']) and all(result['raw_retained_unsafe_projection_refused']),name
    print(json.dumps({'counts':{k:len(packet[k]) for k in ('valid','refused','projection_refused')},'candidate_surfaces':list(candidates),'openbabel_schema_valid':sum(e['default_self_output_valid_declared_schema'] for e in external),'report':str(args.output)}))

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--python-path',type=Path,required=True);p.add_argument('--wasm-module',type=Path,required=True);p.add_argument('--cli',type=Path,required=True)
    p.add_argument('--openbabel',default='/opt/homebrew/bin/obabel');p.add_argument('--output',type=Path,default=ROOT/'validation/lpmd_preservation/observed.json')
    observe(p.parse_args())
if __name__=='__main__': main()
