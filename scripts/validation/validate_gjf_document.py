#!/usr/bin/env python3
"""Exact all-job source and fixture-authored geometry, not job execution validity."""
from __future__ import annotations
import argparse
from collections import Counter
from datetime import datetime, timezone
import gzip
import hashlib
import json
from pathlib import Path
import platform
import struct
import subprocess
import sys
import tempfile

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT/'scripts/validation'))
from validate_gjf_source import run, serial, source_result
from check_interchange_preservation import identity, require, strict_json

WORKER = '''
import json,sys
sys.path.insert(0,sys.argv[1])
import chematic
records=[]
for c in json.load(sys.stdin):
 try:
  if sys.argv[2]=="baseline":
   r={"ok":True,"output":chematic.convert_format(c["source"],"com","gjf")}
  else:
   output=chematic.roundtrip_gjf_document(c["source"])
   r={"ok":True,"output":output,"inspection":chematic.parse_gjf_document_json(c["source"])}
 except (ValueError,TypeError) as e:r={"ok":False,"error":str(e)}
 records.append(dict(id=c["id"],**r))
print(json.dumps(records,ensure_ascii=True,allow_nan=False))
'''


def document_result(case, result):
    status = source_result(case['source'], result)
    if status != 'identical':
        require(status != 'refused' or 'inspection' not in result, 'refusal contains partial jobs')
        return status
    report = result.get('inspection')
    if type(report) is str: report = strict_json(report.encode('utf8'))
    require(type(report) is dict and report.get('schema') == 'chematic.gjf-document.v1', 'invalid document inspection')
    jobs = report.get('jobs')
    require(type(jobs) is list and type(case.get('jobs')) is list, 'invalid job arrays')
    if len(jobs) != len(case['jobs']): return 'mismatched'
    for index, (actual, expected) in enumerate(zip(jobs, case['jobs'])):
        if type(actual) is not dict or type(actual.get('index')) is not int or actual['index'] != index: return 'mismatched'
        if any(type(actual.get(k)) is not int for k in ['source_start','source_end','charge','multiplicity']): return 'mismatched'
        if type(actual.get('freeze_codes')) is not list or any(v is not None and type(v) is not int for v in actual['freeze_codes']): return 'mismatched'
        for key in ['source_start', 'source_end', 'elements', 'freeze_codes', 'charge', 'multiplicity']:
            if actual.get(key) != expected[key]: return 'mismatched'
        coordinates = actual.get('coords')
        if type(coordinates) is not list or len(coordinates) != len(expected['coords']): return 'mismatched'
        for actual_row, expected_row in zip(coordinates, expected['coords']):
            if type(actual_row) is not list or len(actual_row) != 3: return 'mismatched'
            for got, wanted in zip(actual_row, expected_row):
                if type(got) not in (float, int) or struct.pack('>d', got) != struct.pack('>d', wanted): return 'mismatched'
    return 'identical_verified'


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ['python-root', 'wasm-module', 'cli', 'baseline-python-root', 'baseline-cli']:
        parser.add_argument('--'+name, type=Path, required=True)
    parser.add_argument('--obabel', type=Path, default=Path('/opt/homebrew/bin/obabel'))
    parser.add_argument('--output', type=Path, default=ROOT/'validation/results/gjf-document-2026-10-22.json')
    args = parser.parse_args()
    fixture = ROOT/'validation/gjf-document-2026-10-22/cases.json'
    cases = strict_json(fixture.read_bytes())['cases']
    positive = [c for c in cases if c['expected']=='retain']
    counts = Counter(); records = []; raw = {}
    for engine, command, selected in [
        ('python', [sys.executable,'-c',WORKER,args.python_root,'candidate'], cases),
        ('wasm', ['node',ROOT/'scripts/validation/validate_gjf_document_node.mjs',args.wasm_module], cases),
        ('baseline_python', [sys.executable,'-c',WORKER,args.baseline_python_root,'baseline'], positive),
    ]:
        process = run(command, json.dumps(selected).encode())
        require(process['exit_code']==0, f'worker failed: {serial(process)}')
        observed = strict_json(process['stdout']); require(type(observed) is list and len(observed)==len(selected), 'worker dropped cases')
        raw[engine] = serial(process)
        for case, result in zip(selected, observed):
            require(result['id']==case['id'], 'worker changed case order')
            status = source_result(case['source'], result) if engine.startswith('baseline') else document_result(case,result)
            if not engine.startswith('baseline'):
                require(status==('identical_verified' if case['expected']=='retain' else 'refused'), f'{engine}: {case["id"]}: {status}')
            counts[f'{engine}.{status}'] += 1
            records.append(dict(engine=engine, **result))
    for case in cases:
        retained = run([args.cli,'gjf-document'],case['source'].encode())
        inspected = run([args.cli,'gjf-document','--inspect'],case['source'].encode())
        if case['expected']=='retain':
            require(retained['exit_code']==inspected['exit_code']==0,'CLI operation failed')
            result = dict(ok=True,output=retained['stdout'].decode(),inspection=strict_json(inspected['stdout']))
        else:
            require(retained['exit_code']==inspected['exit_code']==2 and not retained['stdout'] and not inspected['stdout'],'CLI emitted partial success or crashed')
            result = dict(ok=False,error=retained['stderr'].decode())
        status = document_result(case,result)
        require(status==('identical_verified' if case['expected']=='retain' else 'refused'),'CLI field/source mismatch')
        counts[f'cli.{status}'] += 1
        records.append(dict(engine='cli',id=case['id'],retained_process=serial(retained),inspection_process=serial(inspected)))
    for case in positive:
        r = run([args.baseline_cli,'convert','--input-format','com','--output-format','gjf'],case['source'].encode())
        status = 'identical' if r['exit_code']==0 and r['stdout']==case['source'].encode() else 'refused' if r['exit_code']==2 and not r['stdout'] else 'unexpected'
        require(status!='unexpected','baseline CLI crashed or changed unexpectedly')
        counts[f'baseline_cli.{status}'] += 1
        records.append(dict(engine='baseline_cli',id=case['id'],process=serial(r)))
    file_checks = []
    with tempfile.TemporaryDirectory(prefix='chematic-gjf-document-') as temp:
        destination = Path(temp)/'out.gjf'; source = Path(temp)/'in.com'
        for case in positive:
            source.write_bytes(case['source'].encode())
            r = run([args.cli,'gjf-document','--input',source,'--output',destination])
            require(r['exit_code']==0 and not r['stdout'] and destination.read_bytes()==source.read_bytes(),'file identity failed')
            file_checks.append(dict(id=case['id'],check='file_identity',process=serial(r),output_sha256=hashlib.sha256(destination.read_bytes()).hexdigest()))
        for case in [c for c in cases if c['expected']=='refuse']:
            destination.write_bytes(b'existing output')
            r = run([args.cli,'gjf-document','--output',destination],case['source'].encode())
            require(r['exit_code']==2 and not r['stdout'] and destination.read_bytes()==b'existing output','failed document replaced output')
            file_checks.append(dict(id=case['id'],check='refusal_preserves_destination',process=serial(r),preserved=True))
        for extra in [['--max-jobs','1'],['--max-input-bytes','8']]:
            destination.write_bytes(b'existing output')
            r = run([args.cli,'gjf-document','--output',destination,*extra],positive[0]['source'].encode())
            require(r['exit_code']==2 and not r['stdout'] and destination.read_bytes()==b'existing output','limit failure replaced output')
            file_checks.append(dict(options=extra,check='limit_preserves_destination',process=serial(r),preserved=True))
    baseline_commit = '296ddd3f13607501456df9fef6efac98c5568d26'
    attested = strict_json(subprocess.check_output(['git','show',f'{baseline_commit}:validation/results/gjf-source-2026-10-22.json'],cwd=ROOT))
    baseline_extension = next((args.baseline_python_root/'chematic').glob('chematic*.so'))
    require(identity(baseline_extension)['sha256']==attested['artifacts']['python_extension']['sha256'],'baseline extension changed')
    require(identity(args.baseline_cli)['sha256']==attested['artifacts']['cli']['sha256'],'baseline CLI changed')
    baseline_wrapper = args.baseline_python_root/'chematic/__init__.py'
    require(baseline_wrapper.read_bytes()==subprocess.check_output(['git','show',f'{baseline_commit}:crates/chematic-py/python/chematic/__init__.py'],cwd=ROOT),'baseline wrapper changed')
    catalog = run([args.obabel,'-L','gjf']); attempt = run([args.obabel,'-igjf','-ogjf'],positive[0]['source'].encode())
    require(b'Write-only' in catalog['stdout']+catalog['stderr'] and attempt['exit_code']!=0 and b'Invalid input format' in attempt['stderr'],'Open Babel input capability changed')
    packet = args.output.with_suffix('.packet.json.gz')
    packet.write_bytes(gzip.compress(json.dumps(dict(records=records,worker_outputs=raw,file_checks=file_checks,openbabel_catalog=serial(catalog),openbabel_attempt=serial(attempt)),indent=2,allow_nan=False).encode(),mtime=0))
    sources = ['crates/chematic-mol/src/gaussian/document.rs','crates/chematic-mol/src/gaussian/source.rs','crates/chematic-mol/src/gaussian.rs','crates/chematic-mol/src/lib.rs','crates/chematic-py/src/formats/gaussian.rs','crates/chematic-py/src/formats.rs','crates/chematic-wasm/src/mol_io/gjf_document.rs','crates/chematic-wasm/src/mol_io.rs','crates/chematic-cli/src/main.rs','crates/chematic-cli/src/conversion_format.rs','crates/chematic-cli/src/dispatch.rs','scripts/validation/validate_gjf_document.py','scripts/validation/validate_gjf_document_node.mjs']
    report = dict(schema='chematic.gjf-document-validation.v1',utc=datetime.now(timezone.utc).isoformat(),platform=platform.platform(),python_version=sys.version,node_version=serial(run(['node','--version'])),rust_version=serial(run(['rustc','--version'])),scope='All explicit Cartesian jobs: exact UTF-8 source and fixture-authored ordered geometry. No checkpoint resolution, Gaussian execution or toolkit-wide superiority claim.',counts=dict(sorted(counts.items())),fixture=identity(fixture),packet=identity(packet),sources=[identity(ROOT/p) for p in sources],artifacts=dict(python_extension=identity(next((args.python_root/'chematic').glob('chematic*.so'))),python_wrapper=identity(args.python_root/'chematic/__init__.py'),wasm_js=identity(args.wasm_module),wasm_binary=identity(args.wasm_module.with_name('chematic_wasm_bg.wasm')),cli=identity(args.cli),openbabel=identity(args.obabel.resolve())),baseline=dict(commit=baseline_commit,scope='Retained previous PR candidate, not main/release. Common same-format route compared with explicit document API; unavailable multi-job input is a coverage difference.',extension=identity(baseline_extension),wrapper=identity(baseline_wrapper),cli=identity(args.baseline_cli),attestation=f'{baseline_commit}:validation/results/gjf-source-2026-10-22.json'),openbabel=dict(version=serial(run([args.obabel,'-V'])),scope='GJF input unsupported: coverage only, not per-field loss or malformed-job quality.'),files=dict(positive_identity=len(positive),refusal_preserves_destination=len(cases)-len(positive),limit_failure_preserves_destination=2),process_limits=dict(seconds=30,combined_pipe_bytes=1048576,memory='not imposed'))
    args.output.write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(report['counts'],indent=2))


if __name__ == '__main__': main()
