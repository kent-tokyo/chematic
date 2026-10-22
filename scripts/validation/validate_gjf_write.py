#!/usr/bin/env python3
"""Paired Cartesian GJF writing; independently observe decimal values and flags."""
from __future__ import annotations
import argparse
from collections import Counter
from datetime import datetime, timezone
import gzip
import hashlib
import importlib
import json
import math
from pathlib import Path
import re
import struct
import subprocess
import sys

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / 'scripts'))
from check_interchange_preservation import bounded_run, identity, require


def observe(text):
    """Read this writer profile independently of chematic; no bond perception."""
    sections, current = [], []
    for line in text.splitlines():
        if not line.strip():
            if current: sections.append(current); current=[]
        else: current.append(line)
    if current: sections.append(current)
    require(len(sections)==3 and len(sections[0])==1 and len(sections[1])==1, 'unexpected GJF sections')
    route, title = sections[0][0], sections[1][0]
    require(route.startswith('#') and '#' not in route[1:] and '!' not in route, 'invalid route')
    units = re.findall(r'\bUnits\s*=\s*(Ang|AU)\b', route, flags=re.I)
    require(len(units)<=1 and len(re.findall(r'\bUnits\b',route,flags=re.I))==len(units), 'conflicting or unsupported unit keywords')
    state = sections[2][0].split()
    require(len(state)==2 and all(re.fullmatch(r'[+-]?\d+',v) for v in state), 'invalid charge/spin')
    charge, multiplicity = map(int,state)
    require(-2**31<=charge<2**31 and 0<multiplicity<2**32, 'out-of-profile state')
    elements, coords, codes = [], [], []
    for row in sections[2][1:]:
        fields = row.split()
        require(len(fields) in (4,5) and re.fullmatch(r'[A-Z][a-z]?',fields[0]), 'invalid atom row')
        code=None
        if len(fields)==5:
            require(re.fullmatch(r'[+-]?\d+',fields[1]), 'invalid freeze code')
            code=int(fields[1]); require(-2**31<=code<2**31, 'freeze code outside i32')
        point = list(map(float,fields[-3:]))
        require(all(math.isfinite(v) for v in point), 'nonfinite coordinate')
        elements.append(fields[0]); coords.append(point); codes.append(code)
    require(elements, 'no atoms')
    return dict(elements=elements,coords=coords,freeze_codes=codes,charge=charge,multiplicity=multiplicity,
                units='bohr' if units and units[0].lower()=='au' else 'angstrom',title=title,route=route)


def same_numbers(actual, expected):
    return (len(actual)==len(expected) and all(len(a)==len(b)==3 and
            all(struct.pack('>d',x)==struct.pack('>d',y) for x,y in zip(a,b)) for a,b in zip(actual,expected)))


def run(command,text):
    r=bounded_run([str(x) for x in command],text.encode(),timeout=10,max_output=1_048_576)
    require(r['resource_status']=='completed','process budget exceeded')
    return {**r, **{key:r[key].decode(errors='replace') for key in ('stdout','stderr')}}


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--obabel',type=Path,default=Path('/opt/homebrew/bin/obabel'))
    parser.add_argument('--output',type=Path,default=ROOT/'validation/results/gjf-write-2026-10-22.json')
    args=parser.parse_args()
    import chematic
    fixture=ROOT/'validation/gjf-cartesian-2026-10-22/write-cases.json'
    cases=json.loads(fixture.read_text())
    records=[];counts=Counter();maximum=Counter()
    for case in cases['common']:
        xyz=str(len(case['elements']))+'\noriginal geometry\n'+''.join(f'{e} {x!r} {y!r} {z!r}\n' for e,(x,y,z) in zip(case['elements'],case['coords']))
        mol,_=chematic.from_xyz(xyz)
        for engine in ('checked','legacy_api','openbabel'):
            raw=None
            if engine=='checked':text=chematic.write_gjf_cartesian(case['elements'],case['coords'],charge=0,multiplicity=1,units='angstrom',route='# HF/STO-3G')
            elif engine=='legacy_api':text=chematic.write_gjf(mol,case['coords'],0,1,'HF/STO-3G','geometry')
            else:
                raw=run([args.obabel,'-ixyz','-ogjf','-ab','-xk','# HF/STO-3G Units=Ang'],xyz)
                require(raw['exit_code']==0,'Open Babel writer failed');text=raw['stdout']
            observed=observe(text)
            require(observed['elements']==case['elements'],'ordered species mismatch')
            exact=same_numbers(observed['coords'],case['coords'])
            counts[f'{engine}.common.bitwise_match' if exact else f'{engine}.common.rounded']+=1
            error=max(abs(a-b) for got,want in zip(observed['coords'],case['coords']) for a,b in zip(got,want))
            maximum[engine]=max(maximum[engine],error)
            if engine=='checked': require(exact and observed['units']=='angstrom','checked writer changed numeric values/units')
            records.append(dict(id=case['id'],lane='common',engine=engine,input_xyz=xyz,output=text,bitwise_match=exact,max_coordinate_error=error,process=raw))
    for hexadecimal in cases['numeric_hex']:
        value=float.fromhex(hexadecimal);point=[[value,-value,0.]]
        text=chematic.write_gjf_cartesian(['C'],point,charge=0,multiplicity=1,units='angstrom')
        require(same_numbers(observe(text)['coords'],point),'numeric control changed f64 bits')
        counts['checked.numeric_profile.bitwise_match']+=1
        records.append(dict(id=hexadecimal,lane='numeric_profile',output=text))
    for units in ('angstrom','bohr'):
        for codes in (None,[None,None,None],[0,-1,None],[-2**31,2**31-1,0]):
            points=[[1.2345678901234567,-0.,1e-200],[2.,3.,4.],[-4.,5.,6.]]
            text=chematic.write_gjf_cartesian(['C','O','H'],points,charge=-1,multiplicity=2,units=units,freeze_codes=codes,title='明示 geometry',route='#p HF/STO-3G')
            got=observe(text);expected=codes if codes is not None else [None]*3
            require(same_numbers(got['coords'],points) and got['freeze_codes']==expected and got['units']==units and got['charge']==-1 and got['multiplicity']==2 and got['title']=='明示 geometry','metadata control mismatch')
            # A second read checks composition of the two public APIs.
            back=chematic.parse_gjf_geometry(text)
            require(same_numbers(back['coords'],points) and back['freeze_codes']==expected and back['elements']==['C','O','H'],'public readback mismatch')
            counts['checked.metadata_profile.match']+=1
            records.append(dict(id=f'{units}:{codes}',lane='metadata_profile',output=text))
    invalid=[('missing-coordinates',dict(coords=[])),('extra-coordinates',dict(coords=[[1,2,3],[4,5,6]])),('missing-code',dict(freeze_codes=[])),('extra-code',dict(freeze_codes=[None,0])),('nan',dict(coords=[[math.nan,0,0]])),('infinity',dict(coords=[[0,math.inf,0]])),('zero-spin',dict(multiplicity=0)),('newline-title',dict(title='title\n\n0 1')),('units-route',dict(route='# Units=AU')),('comment-route',dict(route='# ! comment'))]
    for name,overrides in invalid:
        kwargs=dict(charge=0,multiplicity=1,units='angstrom');kwargs.update(overrides)
        points=kwargs.pop('coords',[[1,2,3]])
        try: output=chematic.write_gjf_cartesian(['C'],points,**kwargs)
        except ValueError as error:
            counts['checked.invalid.refused']+=1;records.append(dict(id=name,lane='invalid',error=str(error),output=None))
        else: raise ValueError(f'unsafe checked output for {name}: {output}')
    # Attest the legacy control to main by hashing the unchanged function bodies,
    # rather than describing this current-candidate runtime as a main rebuild.
    base='654993b832f1431ce4137675102368340b049981'
    body_hashes={}
    for path,start,end in [('crates/chematic-mol/src/gaussian.rs','pub fn write_gjf(','// Gaussian LOG parser'),('crates/chematic-py/src/formats.rs','fn write_gjf(','/// Return True if the SMILES')]:
        baseline=subprocess.check_output(['git','show',f'{base}:{path}'],cwd=ROOT,text=True)
        current=(ROOT/path).read_text()
        old=baseline.split(start,1)[1].split(end,1)[0];new=current.split(start,1)[1].split(end,1)[0]
        require(old==new,'legacy writer code changed; review baseline description')
        body_hashes[path]=hashlib.sha256((start+new).encode()).hexdigest()
    packet_path=args.output.with_suffix('.packet.json.gz')
    packet_path.write_bytes(gzip.compress(json.dumps(records,indent=2,allow_nan=False).encode(),mtime=0))
    report=dict(schema='chematic.gjf-write-validation.v1',utc=datetime.now(timezone.utc).isoformat(),scope='Local candidate geometry writing. Not graph/job/source preservation, Gaussian execution, or toolkit-wide superiority.',counts=dict(counts),max_coordinate_error=dict(maximum),fixture=identity(fixture),packet=identity(packet_path),legacy_control=dict(scope='Unchanged legacy writer API in candidate runtime; function bodies equal source main, not a separate main binary.',source_main=base,function_body_sha256=body_hashes),artifacts=dict(extension=identity(importlib.import_module('chematic.chematic').__file__),openbabel=identity(args.obabel.resolve())),sources=[identity(ROOT/path) for path in ['crates/chematic-mol/src/gaussian/writer.rs','crates/chematic-py/src/formats/gaussian/writer.rs','scripts/validation/validate_gjf_write.py']],openbabel_version=run([args.obabel,'-V'],''),limits=dict(external_seconds=10,combined_pipe_bytes=1_048_576,memory='not imposed'),comparison_notes=['Common lane uses identical ordered elements/coordinate values and angstrom units; legacy/Open Babel precision differences are observed independently.','Freeze/unit/signed-zero/subnormal controls are separate candidate-only profile lanes; not scored as Open Babel source-retention failures.','Malformed input is a candidate-only refusal lane. Open Babel is not used as the semantic oracle.'])
    args.output.write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(dict(counts=report['counts'],max_coordinate_error=report['max_coordinate_error']),indent=2))


if __name__=='__main__': main()
