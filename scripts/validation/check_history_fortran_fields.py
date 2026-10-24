#!/usr/bin/env python3
"""Read the writer's real fields with independent Fortran F/E/G descriptors.

Requires a Fortran compiler. This checks real-field numeric compatibility, not
execution of a DL_POLY simulator or its postprocessing programs.
"""
import argparse
import copy
import hashlib
import json
from pathlib import Path
import struct
import subprocess
import tempfile

ROOT=Path(__file__).resolve().parents[2]


def fields(text, record):
    rows=iter(text.splitlines())
    next(rows);next(rows)
    result=[]
    def take(row, offset, count):
        result.extend(row[offset+12*i:offset+12*(i+1)] for i in range(count))
    for frame in record['frames']:
        take(next(rows),48,1)
        if frame['cell_vectors'] is not None:
            for _ in range(3):take(next(rows),0,3)
        for _ in frame['atoms']:
            take(next(rows),18,2)
            for _ in range(record['keytrj']+1):take(next(rows),0,3)
    if list(rows) or any(len(s)!=12 for s in result):
        raise ValueError('unexpected fixed-width output layout')
    return result


def values(record):
    result=[]
    for frame in record['frames']:
        result.append(frame['timestep'])
        for v in frame['cell_vectors'] or []:result.extend(v)
        for a in frame['atoms']:
            result.extend([a['mass'],a['charge'],*a['position']])
            if a['velocity'] is not None:result.extend(a['velocity'])
            if a['force'] is not None:result.extend(a['force'])
    return result


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--cli',type=Path,required=True)
    parser.add_argument('--compiler',default='gfortran')
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args()
    packet=json.loads((ROOT/'validation/dlpoly_history/cases.json').read_text(encoding='utf-8'))
    records=[c['expected'] for c in packet['valid']]
    for value in (-0.0,1e-100,1e100):
        r=copy.deepcopy(records[0]);r['frames'][0]['atoms'][0]['position'][0]=value;records.append(r)
    raw=[];expected=[]
    for record in records:
        p=subprocess.run([str(args.cli.resolve()),'dlpoly-history','--from-json'],input=json.dumps(record,allow_nan=False),text=True,capture_output=True,check=True)
        raw.extend(fields(p.stdout,record));expected.extend(values(record))
    with tempfile.TemporaryDirectory(prefix='chematic-fortran-fields-') as directory:
        binary=Path(directory)/'probe'
        source=ROOT/'validation/dlpoly_history/real_field_probe.f90'
        subprocess.run([args.compiler,'-std=f2008','-O0','-o',str(binary),str(source)],capture_output=True,text=True,check=True)
        p=subprocess.run([str(binary)],input='\n'.join(raw)+'\n',text=True,capture_output=True,check=True)
    observed=[[float(v) for v in line.split()] for line in p.stdout.splitlines()]
    if len(observed)!=len(expected):raise ValueError('Fortran result count mismatch')
    for index,(want,got) in enumerate(zip(expected,observed)):
        if len(got)!=3 or any(struct.pack('>d',v)!=struct.pack('>d',want) for v in got):
            raise ValueError(f'Fortran F/E/G numeric mismatch at {index}: {want!r} -> {got!r}')
    report=dict(schema_version=1,scope='Independent real-field numeric compatibility under Fortran F12.6/E12.4/G12.4; not a simulation or full postprocessor test',
                records_checked=len(records),real_fields_checked=len(expected),descriptors=['F12.6','E12.4','G12.4'],all_f64_bits_equal=True,
                formatted_reads_checked=3*len(expected),
                packet_sha256=hashlib.sha256((ROOT/'validation/dlpoly_history/cases.json').read_bytes()).hexdigest(),
                probe_sha256=hashlib.sha256((ROOT/'validation/dlpoly_history/real_field_probe.f90').read_bytes()).hexdigest(),
                script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                cli_sha256=hashlib.sha256(args.cli.resolve().read_bytes()).hexdigest(),
                real_field_input_sha256=hashlib.sha256(('\n'.join(raw)+'\n').encode()).hexdigest(),
                compiler=subprocess.check_output([args.compiler,'--version'],text=True).splitlines()[0])
    args.output.write_text(json.dumps(report,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(report,indent=2))


if __name__=='__main__':raise SystemExit(main())
