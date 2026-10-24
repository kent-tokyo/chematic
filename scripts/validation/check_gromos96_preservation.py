#!/usr/bin/env python3
"""Observe GROMOS96 records through Python, actual Node/WASM and CLI.

External ASE geometry and Open Babel's one-way writer are separate lanes. File
unit context is supplied explicitly; atom names never establish element truth.
"""
from __future__ import annotations
import argparse
from decimal import Decimal
import hashlib
import importlib
import io
import json
import math
from pathlib import Path
import subprocess
import sys
import tempfile

ROOT=Path(__file__).resolve().parents[2]
PACKET=ROOT/'validation/gromos96_preservation/cases.json'

def digest(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as f:
        for chunk in iter(lambda:f.read(1<<20),b''):h.update(chunk)
    return h.hexdigest()

def run(args,text=None):
    return subprocess.run(list(map(str,args)),input=text,text=True,encoding='utf-8',capture_output=True,timeout=30,check=False)

def numbers(line,full=False,decimal=False):
    payload=line[24:] if full else line
    width={45:15,48:16}.get(len(payload))
    if width is None:raise ValueError('unsupported numeric row width')
    constructor=Decimal if decimal else float
    return [constructor(payload[width*i:width*(i+1)].strip()) for i in range(3)]

def project(record):
    unit=record['length_unit']
    if unit not in ('nm','angstrom'):raise ValueError('unit is unspecified')
    factor=10. if unit=='nm' else 1.
    frames=[]
    for b in record['blocks']:
        data=[s for s in b['lines'] if not s.startswith('#')]
        if b['kind'] in ('POSITION','POSITIONRED'):
            full=b['kind']=='POSITION'
            identities=[dict(zip(('residue_number','residue_name','atom_name','atom_number'),[s[:5].strip(),s[6:11].strip(),s[12:17].strip(),s[17:24].strip()])) for s in data] if full else None
            frames.append({'identities':identities,'positions_angstrom':[[v*factor for v in numbers(s,full)] for s in data],'box_vectors_angstrom':None})
        elif b['kind']=='BOX':
            if len(data)!=1 or len(data[0]) not in (45,135):raise ValueError('unsupported BOX')
            v=[float(data[0][i:i+15])*factor for i in range(0,len(data[0]),15)]
            cell=[[v[0],0.,0.],[0.,v[1],0.],[0.,0.,v[2]]]
            if len(v)==9:cell[0][1],cell[0][2],cell[1][0],cell[1][2],cell[2][0],cell[2][1]=v[3:]
            frames[-1]['box_vectors_angstrom']=cell
    for frame in frames:
        values=frame['positions_angstrom']+(frame['box_vectors_angstrom'] or [])
        if any(not math.isfinite(v) for row in values for v in row):raise ValueError('projection is nonfinite')
    return frames

def bits(value):
    if isinstance(value,list):return [bits(v) for v in value]
    if isinstance(value,dict):return {k:bits(v) for k,v in value.items()}
    return float(value).hex() if isinstance(value,(float,int)) else value

def inspect_text(text,unit):
    record={'length_unit':unit,'blocks':[],'trailing':[]};active=None;leading=[]
    lines=text.split('\n')
    if lines[-1]=='':lines.pop()
    for line in lines:
        line=line.removesuffix('\r')
        if active is not None:
            if line.strip()=='END':record['blocks'].append(active);active=None
            else:active['lines'].append(line)
        elif not line.strip() or line.startswith('#'):leading.append(line)
        elif line.strip() in ('TITLE','TIMESTEP','POSITION','POSITIONRED','VELOCITY','VELOCITYRED','BOX'):
            active={'kind':line.strip(),'leading':leading,'lines':[]};leading=[]
        else:raise ValueError('unknown block')
    if active is not None or not record['blocks'] or record['blocks'][0]['kind']!='TITLE':raise ValueError('incomplete document')
    if not any(b['kind'] in ('POSITION','POSITIONRED') for b in record['blocks']):raise ValueError('no coordinate frames')
    record['trailing']=leading
    return record

def ase_geometry(read,text,record):
    expected=project(record)
    if len(expected)!=1:return {'observable':False,'reason':'single-Atoms reader does not expose the frame sequence'}
    try:
        atoms=read(io.StringIO(text),format='gromos')
        # ASE assumes nm. Adjust only with the declared caller context; never guess.
        correction=1. if record['length_unit']=='nm' else .1
        positions=(atoms.positions*correction).tolist();cell=(atoms.cell.array*correction).tolist()
        target=expected[0];shape=len(positions)==len(target['positions_angstrom'])
        errors=[abs(a-b) for row_a,row_b in zip(positions,target['positions_angstrom']) for a,b in zip(row_a,row_b)]
        box_errors=[abs(a-b) for row_a,row_b in zip(cell,target['box_vectors_angstrom']) for a,b in zip(row_a,row_b)] if target['box_vectors_angstrom'] is not None else []
        finite=all(math.isfinite(v) for row in positions+cell for v in row)
        scale=max([1.]+[abs(v) for row in target['positions_angstrom']+(target['box_vectors_angstrom'] or []) for v in row])
        maximum=max(errors+box_errors,default=0.)
        return {'observable':True,'count_match':shape,'finite':finite,'max_abs_geometry_error_angstrom':maximum,'within_2e_14_relative':shape and finite and maximum<=2e-14*scale,
                'scope':'positions and declared BOX only; no element/label, title, time, velocity or frame-sequence conclusion'}
    except (ValueError,RuntimeError,IndexError,KeyError) as e:return {'observable':False,'reason':str(e)}

def decimal_positions(record):
    if record['length_unit'] not in ('nm','angstrom'):raise ValueError('unit is unspecified')
    factor=Decimal(10) if record['length_unit']=='nm' else Decimal(1)
    return [[[v*factor for v in numbers(s,b['kind']=='POSITION',True)] for s in b['lines'] if not s.startswith('#')] for b in record['blocks'] if b['kind'] in ('POSITION','POSITIONRED')]

def writer_record(probe,unit):
    def field(value):
        value=Decimal(value)/(Decimal(10) if unit=='nm' else Decimal(1))
        choices=[str(value),format(value,'f'),format(value.normalize(),'E')]
        for candidate in list(choices):
            if candidate.startswith('0.'):choices.append(candidate[1:])
            elif candidate.startswith(('-0.','+0.')):choices.append(candidate[0]+candidate[2:])
        token=next((s for s in choices if len(s)<=15),None)
        if token is None:raise ValueError('probe cannot fit without rounding')
        return token.rjust(15)
    rows=[]
    for i,(element,position) in enumerate(zip(probe['elements'],probe['positions_angstrom']),1):
        rows.append(f'{1:5d} {"MOL":<5} {element:<5}{i:7d}'+''.join(field(v) for v in position))
    return {'length_unit':unit,'blocks':[{'kind':'TITLE','leading':[],'lines':['Original synthetic geometry']},{'kind':'POSITION','leading':[],'lines':rows}],'trailing':[]}

NODE=r'''
const fs=require('node:fs'),w=require(process.argv[1]),p=JSON.parse(fs.readFileSync(0,'utf8'));
const out={valid:[],refused:[],unsafe:[]};
for(const c of p.valid){let r=w.parse_gromos96_json(c.text,c.record.length_unit),states=[],geometry=[];for(let i=0;i<4;i++){states.push(JSON.parse(r));geometry.push(w.gromos96_geometry_angstrom_json(r));r=w.parse_gromos96_json(w.write_gromos96_from_json(r),c.record.length_unit);}out.valid.push({states,geometry});}
for(const c of p.refused){try{w.parse_gromos96_json(c.text,'unspecified');out.refused.push(false);}catch(e){out.refused.push(String(e).startsWith('GROMOS96 line '));}}
for(const c of p.projection_refused){let r=w.parse_gromos96_json(c.text,c.record.length_unit),refused=false;try{w.gromos96_geometry_angstrom_json(r);}catch(e){refused=String(e).startsWith('GROMOS96 line ');}out.unsafe.push({raw:JSON.parse(w.parse_gromos96_json(w.write_gromos96_from_json(r),c.record.length_unit)),refused});}
process.stdout.write(JSON.stringify(out));
'''

def observe(args):
    packet=json.loads(PACKET.read_text(encoding='utf-8'))
    sys.path[:0]=[str(args.python_path.resolve()),str(args.ase_path.resolve())]
    module=importlib.import_module('chematic');ase=importlib.import_module('ase');read=importlib.import_module('ase.io').read
    py={'valid':[],'refused':[],'unsafe':[]};cli={'valid':[],'refused':[],'unsafe':[]}
    for c in packet['valid']:
        unit=c['record']['length_unit'];r=module.parse_gromos96(c['text'],length_unit=unit);p={'states':[],'geometry':[]};q={'states':[],'geometry':[]};text=c['text']
        for _ in range(4):
            p['states'].append(r);p['geometry'].append(module.gromos96_geometry_angstrom(r));r=module.parse_gromos96(module.write_gromos96(r),length_unit=unit)
            state=run([args.cli,'gromos96','--inspect','--length-unit',unit],text);geometry=run([args.cli,'gromos96','--angstrom','--length-unit',unit],text)
            if state.returncode or geometry.returncode:raise RuntimeError(state.stderr+geometry.stderr)
            q['states'].append(json.loads(state.stdout));q['geometry'].append(json.loads(geometry.stdout))
            written=run([args.cli,'gromos96','--from-json'],state.stdout)
            if written.returncode:raise RuntimeError(written.stderr)
            text=written.stdout
        py['valid'].append(p);cli['valid'].append(q)
    for c in packet['refused']:
        try:module.parse_gromos96(c['text']);py['refused'].append(False)
        except ValueError as e:py['refused'].append(str(e).startswith('GROMOS96 line '))
        out=run([args.cli,'gromos96'],c['text']);cli['refused'].append(out.returncode==2 and not out.stdout and 'GROMOS96 line ' in out.stderr)
    for c in packet['projection_refused']:
        unit=c['record']['length_unit'];r=module.parse_gromos96(c['text'],length_unit=unit);entry={'raw':module.parse_gromos96(module.write_gromos96(r),length_unit=unit),'refused':False}
        try:module.gromos96_geometry_angstrom(r)
        except ValueError as e:entry['refused']=str(e).startswith('GROMOS96 line ')
        py['unsafe'].append(entry);out=run([args.cli,'gromos96','--inspect','--length-unit',unit],c['text']);projection=run([args.cli,'gromos96','--angstrom','--length-unit',unit],c['text'])
        cli['unsafe'].append({'raw':json.loads(out.stdout),'refused':projection.returncode==2 and not projection.stdout})
    node=run(['node','-e',NODE,args.wasm_module.resolve()],json.dumps(packet))
    if node.returncode:raise RuntimeError(node.stderr)
    wasm=json.loads(node.stdout)
    for entry in wasm['valid']:entry['geometry']=[json.loads(s) for s in entry['geometry']]
    candidates={}
    for name,result in [('python',py),('node_wasm',wasm),('cli',cli)]:
        if any(len(result[k])!=len(packet[p]) for k,p in [('valid','valid'),('refused','refused'),('unsafe','projection_refused')]):raise AssertionError('incomplete runtime receipt')
        entries=[]
        for c,e in zip(packet['valid'],result['valid']):entries.append({'id':c['id'],'all_fields_all_four_states':len(e['states'])==4 and all(s==c['record'] for s in e['states']),'geometry_bits_all_four_states':len(e['geometry'])==4 and all(bits(g)==bits(project(c['record'])) for g in e['geometry'])})
        candidates[name]={'valid':entries,'refused':result['refused'],'raw_retained_projection_refused':[e['raw']==c['record'] and e['refused'] for e,c in zip(result['unsafe'],packet['projection_refused'])]}
    ase_records=[{'id':c['id'],'source':ase_geometry(read,c['text'],c['record']),'candidate':ase_geometry(read,module.write_gromos96(c['record']),c['record'])} for c in packet['valid']]
    writers=[]
    with tempfile.TemporaryDirectory(prefix='gromos96-observer-') as directory:
        d=Path(directory)
        for probe in packet['writer_probes']:
            source=d/'in.xyz';output=d/'out.g96';xyz=f'{len(probe["elements"])}\nOriginal synthetic geometry\n'+''.join(e+' '+' '.join(p)+'\n' for e,p in zip(probe['elements'],probe['positions_angstrom']))
            source.write_text(xyz,encoding='utf-8')
            for unit in ('angstrom','nm'):
                output.unlink(missing_ok=True);cmd=[args.openbabel,'-ixyz',source,'-ogr96','-O',output,'-ab']+(['-xn'] if unit=='nm' else [])
                out=run(cmd);text=output.read_text(encoding='utf-8') if output.exists() else ''
                candidate=writer_record(probe,unit);candidate_text=module.write_gromos96(candidate)
                expected=[[[Decimal(v) for v in row] for row in probe['positions_angstrom']]]
                try:external=inspect_text(text,unit);equal=decimal_positions(external)==expected;observed=True;error=None
                except (ValueError,IndexError) as e:external=None;equal=False;observed=False;error=str(e)
                writers.append({'id':probe['id'],'declared_unit':unit,'candidate_decimal_positions_equal':decimal_positions(candidate)==expected,
                    'candidate_ase_geometry':ase_geometry(read,candidate_text,candidate),'openbabel':{'exit':out.returncode,'stderr':out.stderr,'output':text,'observable':observed,'error':error,'decimal_positions_equal':equal,
                    'ase_geometry':ase_geometry(read,text,external) if external else None},'scope':'common author-supplied elements and Cartesian geometry, no source metadata or round-trip claim'})
    paths=['crates/chematic-mol/src/gromos96.rs','crates/chematic-mol/src/gromos96/fields.rs','crates/chematic-mol/src/gromos96/frames.rs','crates/chematic-mol/src/gromos96/text.rs','crates/chematic-py/src/formats/gromos96.rs','crates/chematic-wasm/src/gromos96.rs','crates/chematic-cli/src/gromos96.rs','scripts/validation/check_gromos96_preservation.py','validation/gromos96_preservation/cases.json']
    report={'schema_version':1,'scope':'bounded local GROMOS96 source candidate; caller-owned units; no release/toolkit-wide/engine/speed claim','candidate_base':run(['git','-C',ROOT,'merge-base','HEAD','origin/main']).stdout.strip(),
        'source_sha256':{p:digest(ROOT/p) for p in paths},'runtime_sha256':{'cli':digest(args.cli),'wasm_js':digest(args.wasm_module),'wasm':digest(args.wasm_module.with_name('chematic_wasm_bg.wasm')),'python_extension':{p.name:digest(p) for p in args.python_path.rglob('*.so')},'ase_gromos_reader':digest(Path(ase.__file__).parent/'io/gromos.py')},
        'versions':{'python':sys.version,'node':run(['node','--version']).stdout.strip(),'ase':ase.__version__,'openbabel':run([args.openbabel,'-V']).stdout.strip()},'counts':{k:len(packet[k]) for k in ('valid','refused','projection_refused','writer_probes')},'candidates':candidates,'ase_record_crossread':ase_records,'common_geometry_writers':writers}
    args.output.parent.mkdir(parents=True,exist_ok=True);args.output.write_text(json.dumps(report,indent=2,ensure_ascii=False)+'\n',encoding='utf-8')
    for name,r in candidates.items():
        assert all(e['all_fields_all_four_states'] and e['geometry_bits_all_four_states'] for e in r['valid']),name
        assert all(r['refused']) and all(r['raw_retained_projection_refused']),name
    assert all(w['candidate_decimal_positions_equal'] for w in writers)
    print(json.dumps({'counts':report['counts'],'surfaces':list(candidates),'common_writer_cases':len(writers),'openbabel_decimal_exact':sum(w['openbabel']['decimal_positions_equal'] for w in writers),'report':str(args.output)}))

def main():
    p=argparse.ArgumentParser(description=__doc__)
    for flag in ('python-path','wasm-module','cli','ase-path'):p.add_argument('--'+flag,type=Path,required=True)
    p.add_argument('--openbabel',default='/opt/homebrew/bin/obabel');p.add_argument('--output',type=Path,default=ROOT/'validation/gromos96_preservation/observed.json')
    observe(p.parse_args())
if __name__=='__main__':main()
