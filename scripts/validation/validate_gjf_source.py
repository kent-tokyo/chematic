#!/usr/bin/env python3
"""Bounded source identity measurements; opaque retention is not job validation."""
from __future__ import annotations
import argparse
from collections import Counter
from datetime import datetime,timezone
import gzip
import json
import platform
from pathlib import Path
import subprocess
import sys
import tempfile

ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT))
from scripts.checks.check_interchange_preservation import bounded_run,identity,require,strict_json

PYTHON_WORKER='''
import sys,json
sys.path.insert(0,sys.argv[1])
import chematic
output=[]
for case in json.load(sys.stdin):
 record={"id":case["id"]}
 functions={"common":lambda:chematic.convert_format(case["source"],".COM","gjf")}
 if hasattr(chematic,"roundtrip_gjf_record"):
  functions["direct"]=lambda:chematic.roundtrip_gjf_record(case["source"])
 for name,function in functions.items():
  try:record[name]={"ok":True,"output":function()}
  except (ValueError,TypeError) as error:record[name]={"ok":False,"error":str(error)}
 output.append(record)
print(json.dumps(output,ensure_ascii=True,allow_nan=False))
'''


def run(command,data=b'',max_output=1_048_576):
    r=bounded_run([str(x) for x in command],data,timeout=30,max_output=max_output)
    require(r['resource_status']=='completed',f'process budget exceeded: {command}')
    return r


def serial(r):
    return {**r,**{key:r[key].decode('utf8',errors='replace') for key in ('stdout','stderr')}}


def source_result(source,result):
    """Errors, partial JSON and altered bytes must never look like retention."""
    require(type(result) is dict and type(result.get('ok')) is bool,'invalid worker result')
    if not result['ok']:
        require('output' not in result and isinstance(result.get('error'),str),'refusal must have no output')
        return 'refused'
    require(isinstance(result.get('output'),str),'successful result must contain source text')
    return 'identical' if result['output'].encode('utf8')==source.encode('utf8') else 'altered'


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--python-root',type=Path,required=True)
    parser.add_argument('--baseline-python-root',type=Path,required=True)
    parser.add_argument('--wasm-module',type=Path,required=True)
    parser.add_argument('--cli',type=Path,required=True)
    parser.add_argument('--obabel',type=Path,default=Path('/opt/homebrew/bin/obabel'))
    parser.add_argument('--output',type=Path,default=ROOT/'validation/results/gjf-source-2026-10-22.json')
    args=parser.parse_args()
    fixture=ROOT/'validation/gjf-cartesian-2026-10-22/source-cases.json'
    cases=strict_json(fixture.read_bytes())['cases'];positive=[c for c in cases if c['expected']=='retain']
    counts=Counter();records=[];raw={}
    for name,command,selected in [
        ('candidate_python',[sys.executable,'-c',PYTHON_WORKER,args.python_root],cases),
        ('baseline_python',[sys.executable,'-c',PYTHON_WORKER,args.baseline_python_root],positive),
        ('candidate_wasm',['node',ROOT/'scripts/validation/validate_gjf_source_node.mjs',args.wasm_module],cases)]:
        r=run(command,json.dumps(selected).encode());require(r['exit_code']==0,f'worker failure: {serial(r)}')
        observed=strict_json(r['stdout']);require(len(observed)==len(selected),'worker dropped cases')
        raw[name]=serial(r)
        for case,result in zip(selected,observed):
            require(result['id']==case['id'],'worker changed order')
            for api in ('common','direct') if name!='baseline_python' else ('common',):
                status=source_result(case['source'],result[api])
                counts[f'{name}.{api}.{status}']+=1
                if name!='baseline_python':require(status==('identical' if case['expected']=='retain' else 'refused'),f'{name}: {case["id"]}: {status}')
            records.append(dict(engine=name,**result))
    cli_records=[]
    for case in cases:
        r=run([args.cli,'convert','--input-format','.COM','--output-format','gjf'],case['source'].encode())
        if case['expected']=='retain':
            require(r['exit_code']==0 and r['stdout']==case['source'].encode(),'CLI source identity failed')
            counts['candidate_cli.common.identical']+=1
        else:
            require(r['exit_code']!=0 and not r['stdout'],'CLI refusal emitted source')
            counts['candidate_cli.common.refused']+=1
        cli_records.append(dict(id=case['id'],**serial(r)))
    # Real CLI files: exact bytes on success, and no replacement of a prior
    # output file when conversion or the requested input budget fails.
    files=[]
    with tempfile.TemporaryDirectory(prefix='chematic-gjf-source-') as temporary:
        directory=Path(temporary);destination=directory/'output.gjf'
        for case in positive:
            source=directory/'input.com';source.write_bytes(case['source'].encode())
            r=run([args.cli,'convert','--input-format','com','--output-format','gjf','--input',source,'--output',destination])
            require(r['exit_code']==0 and not r['stdout'] and destination.read_bytes()==source.read_bytes(),'CLI file identity failed')
        files.append(dict(check='all 21 fixture files retain exact UTF-8 bytes',passed=True))
        destination.write_bytes(b'prior output must survive')
        for data,extra in [(b'not a GJF',[]),(positive[0]['source'].encode(),['--max-input-bytes','8'])]:
            r=run([args.cli,'convert','--input-format','gjf','--output-format','com','--output',destination,*extra],data)
            require(r['exit_code']!=0 and not r['stdout'] and destination.read_bytes()==b'prior output must survive','CLI refusal replaced prior file')
            files.append(dict(check='failure leaves prior output untouched',process=serial(r),passed=True))
    catalog=run([args.obabel,'-L','gjf']);ob=run([args.obabel,'-igjf','-ogjf'],positive[0]['source'].encode())
    require(b'Write-only' in catalog['stdout']+catalog['stderr'] and ob['exit_code']!=0 and b'Invalid input format' in ob['stderr'],'Open Babel GJF input capability changed; review comparison')
    # Check the retained pre-change wrapper and binary against the previous PR
    # candidate, not against a guessed release or a different source runtime.
    baseline_commit='9c391f606a84e8cbab8fa11c93f7336c1a2a46bc'
    old=strict_json(subprocess.check_output(['git','show',f'{baseline_commit}:validation/results/gjf-write-2026-10-22.json'],cwd=ROOT))
    wrapper=args.baseline_python_root/'chematic/__init__.py'
    require(wrapper.read_bytes()==subprocess.check_output(['git','show',f'{baseline_commit}:crates/chematic-py/python/chematic/__init__.py'],cwd=ROOT),'baseline wrapper drift')
    baseline_ext=next((args.baseline_python_root/'chematic').glob('chematic*.so'))
    require(identity(baseline_ext)['sha256']==old['artifacts']['extension']['sha256'],'baseline binary attestation mismatch')
    packet_path=args.output.with_suffix('.packet.json.gz')
    packet_path.write_bytes(gzip.compress(json.dumps(dict(worker_outputs=raw,observations=records,cli=cli_records,files=files,openbabel_catalog=serial(catalog),openbabel_attempt=serial(ob)),indent=2,allow_nan=False).encode(),mtime=0))
    report=dict(schema='chematic.gjf-source-validation.v1',utc=datetime.now(timezone.utc).isoformat(),platform=platform.platform(),python_version=sys.version,node_version=serial(run(['node','--version'])),rust_version=serial(run(['rustc','--version'])),scope='Exact source retention for the bounded single leading-Cartesian job profile. Opaque sections are not validated; no Gaussian execution or toolkit-wide superiority claim.',counts=dict(sorted(counts.items())),fixture=identity(fixture),packet=identity(packet_path),baseline=dict(commit=baseline_commit,scope='Retained previous PR candidate, not main or a release. Baseline positives only measure unavailable same-format coverage.',wrapper=identity(wrapper),extension=identity(baseline_ext),attestation=f'{baseline_commit}:validation/results/gjf-write-2026-10-22.json'),artifacts=dict(python_extension=identity(next((args.python_root/'chematic').glob('chematic*.so'))),python_wrapper=identity(args.python_root/'chematic/__init__.py'),wasm_js=identity(args.wasm_module),wasm_binary=identity(args.wasm_module.with_name('chematic_wasm_bg.wasm')),cli=identity(args.cli),openbabel=identity(args.obabel.resolve())),sources=[identity(ROOT/path) for path in ['crates/chematic-mol/src/gaussian/source.rs','crates/chematic-mol/src/gaussian.rs','crates/chematic-mol/src/lib.rs','crates/chematic-py/src/formats/gaussian.rs','crates/chematic-py/python/chematic/__init__.py','crates/chematic-wasm/src/mol_io/common_formats.rs','crates/chematic-cli/src/main.rs','crates/chematic-cli/src/conversion_format.rs','scripts/validation/validate_gjf_source.py','scripts/validation/validate_gjf_source_node.mjs']],openbabel=dict(version=serial(run([args.obabel,'-V'])),comparison='GJF input unsupported: a coverage difference, not a per-field preservation or malformed-input validation score. Usage text on stdout is not molecule output.'),process_limits=dict(seconds=30,combined_pipe_bytes=1_048_576,memory='not imposed'),files=files)
    args.output.write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report['counts'],indent=2))


if __name__=='__main__':main()
