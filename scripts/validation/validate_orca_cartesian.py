#!/usr/bin/env python3
"""Unit-aware ORCA geometry and refusals; opaque job settings are not validated."""
from __future__ import annotations
import argparse
from collections import Counter
from datetime import datetime,timezone
import gzip
import json
from pathlib import Path
import platform
import struct
import subprocess
import sys

ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT))
from scripts.checks.check_interchange_preservation import bounded_run,identity,require,strict_json

WORKER='''
import json,sys
sys.path.insert(0,sys.argv[1])
import chematic
records=[]
for c in json.load(sys.stdin):
 r={"id":c["id"]}
 functions={"core":chematic.parse_orca_input}
 if sys.argv[2]=="candidate":functions["checked"]=chematic.parse_orca_cartesian
 for name,call in functions.items():
  try:r[name]={"ok":True,"value":call(c["source"])}
  except (ValueError,TypeError) as e:r[name]={"ok":False,"error":str(e)}
 records.append(r)
print(json.dumps(records,ensure_ascii=True,allow_nan=False))
'''


def run(command,data=b''):
    r=bounded_run([str(x) for x in command],data,timeout=30,max_output=1048576)
    require(r['resource_status']=='completed','process budget exceeded')
    return r


def serial(r):return {**r,**{k:r[k].decode('utf8',errors='replace') for k in ['stdout','stderr']}}


def returned(result):
    require(type(result) is dict and type(result.get('ok')) is bool,'invalid worker envelope')
    if not result['ok']:
        require('value' not in result and isinstance(result.get('error'),str),'refusal contains partial value')
        return None
    value=result.get('value')
    if type(value) is str:value=strict_json(value.encode())
    require(type(value) is dict,'success lacks a result object')
    return value


def geometry_result(case,result):
    actual=returned(result)
    if actual is None:return 'refused'
    if actual.get('schema')!='chematic.orca-cartesian.v1':return 'untagged'
    expected=case.get('fields')
    if expected is None:return 'unexpected_geometry'
    if type(actual.get('units_explicit')) is not bool or any(type(actual.get(k)) is not int for k in ['charge','multiplicity']):return 'mismatched'
    for key in ['coordinate_units','units_explicit','elements','charge','multiplicity','frozen','atom_extras']:
        if actual.get(key)!=expected[key]:return 'mismatched'
    for row in actual['frozen']:
        if type(row) is not list or len(row)!=3 or any(type(v) is not bool for v in row):return 'mismatched'
    coordinates=actual.get('coords')
    if type(coordinates) is not list or len(coordinates)!=len(expected['coordinate_bits']):return 'mismatched'
    for row,bits in zip(coordinates,expected['coordinate_bits']):
        if type(row) is not list or len(row)!=3:return 'mismatched'
        for value,wanted in zip(row,bits):
            if type(value) not in (float,int) or struct.pack('>d',value).hex()!=wanted:return 'mismatched'
    return 'verified'


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    for name in ['python-root','wasm-module','baseline-python-root','baseline-wasm-module']:
        parser.add_argument('--'+name,type=Path,required=True)
    parser.add_argument('--obabel',type=Path,default=Path('/opt/homebrew/bin/obabel'))
    parser.add_argument('--output',type=Path,default=ROOT/'validation/results/orca-cartesian-2026-10-22.json')
    args=parser.parse_args();fixture=ROOT/'validation/orca-cartesian-2026-10-22/cases.json';cases=strict_json(fixture.read_bytes())['cases']
    counts=Counter();records=[];raw={};core_results={}
    for engine,command in [
        ('candidate_python',[sys.executable,'-c',WORKER,args.python_root,'candidate']),
        ('baseline_python',[sys.executable,'-c',WORKER,args.baseline_python_root,'baseline']),
        ('candidate_wasm',['node',ROOT/'scripts/validation/validate_orca_cartesian_node.mjs',args.wasm_module,'candidate']),
        ('baseline_wasm',['node',ROOT/'scripts/validation/validate_orca_cartesian_node.mjs',args.baseline_wasm_module,'baseline']),
    ]:
        r=run(command,json.dumps(cases).encode());require(r['exit_code']==0,f'worker failed: {serial(r)}');raw[engine]=serial(r)
        result=strict_json(r['stdout']);require(type(result) is list and len(result)==len(cases),'worker dropped cases');core_results[engine]={}
        for case,row in zip(cases,result):
            require(row['id']==case['id'],'worker changed order')
            core=returned(row['core']);core_status='accepted' if core is not None else 'refused';core_results[engine][case['id']]=core_status
            counts[f'{engine}.core.{core_status}']+=1
            if engine.startswith('candidate'):
                status=geometry_result(case,row['checked']);counts[f'{engine}.checked.{status}']+=1
                require(status==('verified' if case['expected']=='accept' else 'refused'),f'{engine}: {case["id"]}: {status}')
                if case.get('core_refuse'):require(core_status=='refused','core selected/accepted invalid geometry')
            if 'legacy_coords' in row:
                legacy=returned(row['legacy_coords']);status='returned_without_unit_tag' if legacy is not None else 'refused'
                counts[f'{engine}.legacy_coords.{status}']+=1
                if engine=='candidate_wasm' and case['expected']=='accept':
                    require((legacy is None)==(case['fields']['coordinate_units']=='bohr'),'legacy Bohr guard failed')
            records.append(dict(engine=engine,**row))
    fixed={surface:[c['id'] for c in cases if c.get('core_refuse') and core_results['baseline_'+surface][c['id']]=='accepted' and core_results['candidate_'+surface][c['id']]=='refused'] for surface in ['python','wasm']}
    old_commit='891bb83fe2966b374ba6d8537d9a429d20ca07d7';base='978669e4264683bb73f5a7a3389a945b5fdf2b87'
    def old_file(commit,path):return subprocess.check_output(['git','show',f'{commit}:{path}'],cwd=ROOT)
    require(old_file(old_commit,'crates/chematic-mol/src/orca.rs')==old_file(base,'crates/chematic-mol/src/orca.rs'),'baseline ORCA source differs from pinned main')
    attested=strict_json(old_file(old_commit,'validation/results/pqr-checked-2026-10-22.json'))
    baseline_ext=next((args.baseline_python_root/'chematic').glob('chematic*.so'))
    baseline_artifacts=dict(python_extension=baseline_ext,wasm_js=args.baseline_wasm_module,wasm_binary=args.baseline_wasm_module.with_name('chematic_wasm_bg.wasm'))
    for name,path in baseline_artifacts.items():require(identity(path)['sha256']==attested['artifacts'][name]['sha256'],'baseline runtime attestation failed')
    require((args.baseline_python_root/'chematic/__init__.py').read_bytes()==old_file(old_commit,'crates/chematic-py/python/chematic/__init__.py'),'baseline wrapper drift')
    catalog=run([args.obabel,'-L','orcainp']);attempt=run([args.obabel,'-iorcainp','-oxyz'],cases[0]['source'].encode())
    require(b'Write-only' in catalog['stdout']+catalog['stderr'] and attempt['exit_code']!=0 and b'Invalid input format' in attempt['stderr'],'Open Babel ORCA input capability changed')
    packet=args.output.with_suffix('.packet.json.gz');packet.write_bytes(gzip.compress(json.dumps(dict(worker_outputs=raw,records=records,openbabel_catalog=serial(catalog),openbabel_attempt=serial(attempt)),indent=2,allow_nan=False).encode(),mtime=0))
    sources=['Cargo.lock','crates/chematic-mol/src/orca.rs','crates/chematic-mol/src/orca/cartesian.rs','crates/chematic-mol/src/lib.rs','crates/chematic-py/src/formats/orca_cartesian.rs','crates/chematic-py/src/formats/materials.rs','crates/chematic-py/src/formats.rs','crates/chematic-wasm/src/format_io.rs','scripts/validation/validate_orca_cartesian.py','scripts/validation/validate_orca_cartesian_node.mjs']
    report=dict(schema='chematic.orca-cartesian-validation.v1',utc=datetime.now(timezone.utc).isoformat(),platform=platform.platform(),python_version=sys.version,node_version=serial(run(['node','--version'])),rust_version=serial(run(['rustc','--version'])),scope='Single embedded Cartesian geometry with input numeric units. Full structured baseline retains raw unit keywords; only the legacy WASM coordinate-only view lacked a unit tag. No conversion, full job validation or toolkit-wide ranking.',counts=dict(sorted(counts.items())),new_core_refusals=fixed,fixture=identity(fixture),packet=identity(packet),sources=[identity(ROOT/p) for p in sources],candidate=dict(base_commit=base,state='working tree at measurement time, not a published release'),artifacts=dict(python_extension=identity(next((args.python_root/'chematic').glob('chematic*.so'))),python_wrapper=identity(args.python_root/'chematic/__init__.py'),wasm_js=identity(args.wasm_module),wasm_binary=identity(args.wasm_module.with_name('chematic_wasm_bg.wasm')),openbabel=identity(args.obabel.resolve())),baseline=dict(commit=old_commit,scope='Retained previous PR candidate, not a main/release rebuild; entire ORCA source verified equal to pinned main. Dependency/runtime identity remains separate.',artifacts={name:identity(path) for name,path in baseline_artifacts.items()},attestation=f'{old_commit}:validation/results/pqr-checked-2026-10-22.json'),openbabel=dict(version=serial(run([args.obabel,'-V'])),scope='ORCA input unsupported: coverage only, not per-field loss or malformed-input quality.'),process_limits=dict(seconds=30,combined_pipe_bytes=1048576,memory='not imposed'))
    args.output.write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report['counts'],indent=2))


if __name__=='__main__':main()
