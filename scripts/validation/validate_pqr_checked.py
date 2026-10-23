#!/usr/bin/env python3
"""Observe PQR atom fields independently; distinguish writer and pipeline lanes."""
from __future__ import annotations
import argparse
from collections import Counter
from datetime import datetime,timezone
import gzip
import json
import math
from pathlib import Path
import platform
import re
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
out=[]
for c in json.load(sys.stdin):
 if c.get("surfaces")==["wasm"]:continue
 atoms=json.loads(c["records_json"])
 row={"id":c["id"]}
 functions={"checked":chematic.write_pqr_checked}
 if c["expected"]=="write":functions["legacy"]=chematic.write_pqr
 for name,call in functions.items():
  try:row[name]={"ok":True,"output":call(atoms)}
  except (ValueError,TypeError,OverflowError) as e:row[name]={"ok":False,"error":str(e)}
 out.append(row)
print(json.dumps(out,ensure_ascii=True,allow_nan=False))
'''


def run(command,data=b''):
    r=bounded_run([str(x) for x in command],data,timeout=30,max_output=1048576)
    require(r['resource_status']=='completed','process budget exceeded')
    return r


def serial(r):
    return {**r,**{k:r[k].decode('utf8',errors='replace') for k in ['stdout','stderr']}}


def observe(text,element_extension=False):
    """APBS whitespace fields; OB's final element extension is explicit, not guessed."""
    records=[]
    for line in text.splitlines():
        tokens=line.split()
        if not tokens or tokens[0] not in ('ATOM','HETATM'):continue
        element=None
        if element_extension:
            if len(tokens) not in (11,12):raise ValueError('invalid element-extension record')
            element=tokens.pop()
        if len(tokens) not in (10,11):raise ValueError('invalid APBS atom record')
        chain=tokens[4] if len(tokens)==11 else None
        offset=6 if chain is not None else 5
        residue=re.fullmatch(r'([+-]?\d+)([A-Za-z]?)',tokens[offset-1])
        if residue is None:raise ValueError('invalid residue token')
        values=[float(v) for v in tokens[offset:]]
        if not all(math.isfinite(v) for v in values):raise ValueError('nonfinite field')
        record=dict(group_pdb=tokens[0],serial=str(int(tokens[1])),atom_name=tokens[2],res_name=tokens[3],chain_id=chain,res_seq=str(int(residue[1])),icode=residue[2] or None,bits={k:struct.pack('>d',v).hex() for k,v in zip(['x','y','z','charge','radius'],values)})
        if element_extension:record['element']=element
        records.append(record)
    if not records:raise ValueError('no atom output')
    return records


def outcome(case,result,element_extension=False):
    require(type(result) is dict and type(result.get('ok')) is bool,'invalid result envelope')
    if not result['ok']:
        require('output' not in result and isinstance(result.get('error'),str),'refusal contains output')
        return dict(status='refused')
    require(isinstance(result.get('output'),str),'success lacks text')
    try:actual=observe(result['output'],element_extension)
    except (ValueError,OverflowError) as error:return dict(status='unreadable',error=str(error))
    expected=case.get('fields',[])
    if len(actual)!=len(expected):return dict(status='mismatched_count',actual=len(actual),expected=len(expected))
    scalar={k:all(a['bits'][k]==e['bits'][k] for a,e in zip(actual,expected)) for k in ['x','y','z','charge','radius']}
    metadata=all(all(a[k]==e[k] for k in ['group_pdb','serial','atom_name','res_name','chain_id','res_seq','icode']) for a,e in zip(actual,expected))
    elements=all(a['element']==e['element'] for a,e in zip(actual,expected)) if element_extension else None
    errors={k:max(abs(struct.unpack('>d',bytes.fromhex(a['bits'][k]))[0]-struct.unpack('>d',bytes.fromhex(e['bits'][k]))[0]) for a,e in zip(actual,expected)) for k in scalar}
    return dict(status='numeric_match' if all(scalar.values()) else 'numeric_mismatch',scalar_bits=scalar,metadata_match=metadata,element_column_match=elements,max_absolute_error=errors,observed=actual)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--python-root',type=Path,required=True)
    parser.add_argument('--wasm-module',type=Path,required=True)
    parser.add_argument('--obabel',type=Path,default=Path('/opt/homebrew/bin/obabel'))
    parser.add_argument('--output',type=Path,default=ROOT/'validation/results/pqr-checked-2026-10-22.json')
    args=parser.parse_args();fixture=ROOT/'validation/pqr-checked-2026-10-22/cases.json'
    packet=strict_json(fixture.read_bytes());cases=packet['cases'];counts=Counter();observations=[];raw={}
    controls=[dict(id=c['id'],category='element-extension-control',expected='write',records_json=json.dumps(c['atoms'],allow_nan=False),fields=c['fields']) for c in packet['openbabel_element_extension_controls']]
    cases=cases+controls
    for engine,command,selected in [('python',[sys.executable,'-c',WORKER,args.python_root],[c for c in cases if c.get('surfaces')!=['wasm']]),('wasm',['node',ROOT/'scripts/validation/validate_pqr_checked_node.mjs',args.wasm_module],cases)]:
        r=run(command,json.dumps(selected).encode());require(r['exit_code']==0,f'worker failed: {serial(r)}');raw[engine]=serial(r)
        results=strict_json(r['stdout']);require(type(results) is list and len(results)==len(selected),'worker dropped cases')
        for case,result in zip(selected,results):
            require(case['id']==result['id'],'worker changed order')
            for api in ['checked','legacy'] if case['expected']=='write' else ['checked']:
                record=outcome(case,result[api]);counts[f'{engine}.{api}.{case["category"]}.{record["status"]}']+=1
                if api=='checked':
                    require(record['status']==('numeric_match' if case['expected']=='write' else 'refused'),f'{engine}: {case["id"]}: {record}')
                    if case['expected']=='write':require(record['metadata_match'],'checked output changed atom identity')
                observations.append(dict(engine=engine,api=api,id=case['id'],**record))
    for lane,selected in [('apbs_pipeline',[c for c in cases if c['category']=='paired']),('element_extension_pipeline',packet['openbabel_element_extension_controls'])]:
        for case in selected:
            source=case['pqr_input'] if lane=='apbs_pipeline' else case['source']
            # All paired inputs must first match the authored expectations under
            # the independent observer, before either engine's output is scored.
            before=outcome(case,dict(ok=True,output=source),lane!='apbs_pipeline')
            require(before['status']=='numeric_match' and before['metadata_match'],'fixture/source disagreement')
            r=run([args.obabel,'-ipqr','-opqr','-ab'],source.encode())
            if r['exit_code']==0:record=outcome(case,dict(ok=True,output=r['stdout'].decode()),True)
            else:record=dict(status='process_error',exit_code=r['exit_code'])
            counts[f'openbabel.{lane}.{record["status"]}']+=1
            observations.append(dict(engine='openbabel',lane=lane,id=case['id'],input_observation=before,process=serial(r),**record))
    baseline='e9e578e242729636535d8f701328fddbab3737c7'
    def legacy_body(text):return text[text.index('pub fn write_pqr('):text.index('\n// ----',text.index('pub fn write_pqr('))]
    require(legacy_body((ROOT/'crates/chematic-mol/src/pqr.rs').read_text())==legacy_body(subprocess.check_output(['git','show',f'{baseline}:crates/chematic-mol/src/pqr.rs'],cwd=ROOT,text=True)),'legacy writer body changed')
    paired_summary={}
    paired_ids={c['id'] for c in packet['cases'] if c['category']=='paired'}
    for engine in ['python','wasm']:
        for api in ['checked','legacy']:
            selected=[r for r in observations if r['engine']==engine and r.get('api')==api and r['id'] in paired_ids]
            paired_summary[f'{engine}.{api}']=dict(numeric_match=sum(r['status']=='numeric_match' for r in selected),cases=len(selected),metadata_match=sum(r.get('metadata_match',False) for r in selected),max_absolute_error={k:max(r['max_absolute_error'][k] for r in selected) for k in ['x','y','z','charge','radius']})
    packet_path=args.output.with_suffix('.packet.json.gz');packet_path.write_bytes(gzip.compress(json.dumps(dict(worker_outputs=raw,observations=observations),indent=2,allow_nan=False).encode(),mtime=0))
    sources=['crates/chematic-mol/src/pqr/checked.rs','crates/chematic-mol/src/pqr.rs','crates/chematic-mol/src/lib.rs','crates/chematic-py/src/formats/pqr_checked.rs','crates/chematic-py/src/formats/materials.rs','crates/chematic-py/src/formats.rs','crates/chematic-wasm/src/format_io/pqr_checked.rs','crates/chematic-wasm/src/format_io.rs','scripts/validation/validate_pqr_checked.py','scripts/validation/validate_pqr_checked_node.mjs']
    report=dict(schema='chematic.pqr-checked-validation.v1',utc=datetime.now(timezone.utc).isoformat(),platform=platform.platform(),python_version=sys.version,node_version=serial(run(['node','--version'])),rust_version=serial(run(['rustc','--version'])),scope='Ordered explicit PQR atom fields. Writer-only checked/legacy API comparison; separate Open Babel read/write pipeline lanes. No original-document, charge-calculation or toolkit-wide superiority claim.',candidate=dict(base_commit=baseline,source_state='Working tree with checked PQR output at measurement time; not a published release'),counts=dict(sorted(counts.items())),paired_writer_summary=paired_summary,fixture=identity(fixture),packet=identity(packet_path),sources=[identity(ROOT/p) for p in sources],artifacts=dict(python_extension=identity(next((args.python_root/'chematic').glob('chematic*.so'))),python_wrapper=identity(args.python_root/'chematic/__init__.py'),wasm_js=identity(args.wasm_module),wasm_binary=identity(args.wasm_module.with_name('chematic_wasm_bg.wasm')),openbabel=identity(args.obabel.resolve())),legacy_control=dict(main_commit=baseline,scope='Legacy API in the same candidate runtime; Rust function body verified equal to pinned main, not a separately rebuilt main/release.'),openbabel=dict(version=serial(run([args.obabel,'-V'])),scope='Same APBS-whitespace input read/write pipeline; element-column dialect controls are separate. Reader and writer effects are not separated. No invalid-input ranking.'),process_limits=dict(seconds=30,combined_pipe_bytes=1048576,memory='not imposed'))
    args.output.write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report['counts'],indent=2))


if __name__=='__main__':main()
