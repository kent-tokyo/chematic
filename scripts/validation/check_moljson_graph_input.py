#!/usr/bin/env python3
"""Observe the local MolJSON graph profile; no Open Babel MolJSON score."""
from __future__ import annotations
import argparse
import hashlib
import importlib
import json
from pathlib import Path
import subprocess
import sys

NODE = """const fs=require('node:fs'),w=require(process.argv[1]);
const cases=JSON.parse(fs.readFileSync(0,'utf8'));
const failure=e=>({status:String(e).startsWith('MolJSON:')?'refused':'error',error:String(e)});
process.stdout.write(JSON.stringify(cases.map(c=>{
 const r={};
 try{const m=w.mol_from_moljson(c.input);try{r.wasm_direct={status:'accepted',output:w.to_moljson(m)}}finally{m.free()}}
 catch(e){r.wasm_direct=failure(e)}
 try{r.wasm_common={status:'accepted',output:w.convert_common_format(c.input,'moljson','moljson')}}
 catch(e){r.wasm_common=failure(e)}
 return r;
})));"""

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
from scripts.checks.check_interchange_preservation import bounded_run, identity, require, strict_json


def snapshot(text):
    value = strict_json(text)
    atoms = value['atoms']
    ids = {a['id']: i for i, a in enumerate(atoms)}
    require(len(ids) == len(atoms), 'duplicate emitted atom IDs')
    for a in atoms:
        require(type(a['id']) is str and bool(a['id']), 'atom ID type')
        require(type(a['charge']) is int and -128 <= a['charge'] <= 127, 'charge type/range')
        require(a['isotope'] is None or type(a['isotope']) is int and 0 <= a['isotope'] <= 65535, 'isotope type/range')
        require(type(a['hydrogens']) is int and 0 <= a['hydrogens'] <= 255, 'H type/range')
        require(type(a['aromatic']) is bool, 'atom aromatic type')
    for b in value['bonds']:
        require(type(b['order']) in {int, float} and b['order'] in {1, 1.5, 2, 3, 4}, 'bond order type/value')
        require(type(b['aromatic']) is bool and b['aromatic'] == (b['order'] == 1.5), 'bond aromatic type/agreement')
    aa = [{k: a[k] for k in ['element', 'charge', 'isotope', 'hydrogens', 'aromatic']} for a in atoms]
    bb = [dict(source=ids[b['source_id']], target=ids[b['target_id']],
               order=b['order'], aromatic=b['aromatic']) for b in value['bonds']]
    return dict(atoms=aa, bonds=bb)


def observe(output, expected):
    try:
        actual = snapshot(output)
        return dict(preserved=actual == expected, snapshot=actual)
    except (ValueError, KeyError, TypeError) as e:
        return dict(preserved=False, error=f'{type(e).__name__}: {e}')


def cli_status(process):
    if process['resource_status'] != 'completed':
        return 'error'
    if process['exit_code'] == 0:
        return 'accepted'
    if (process['exit_code'] == 2 and not process['stdout']
            and process['stderr'].startswith('chematic: MolJSON:')):
        return 'refused'
    return 'error'


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--output', type=Path, required=True)
    p.add_argument('--baseline-ref')
    p.add_argument('--parent-ref', default='HEAD')
    p.add_argument('--cli', type=Path, required=True)
    p.add_argument('--wasm-module', type=Path, required=True)
    args = p.parse_args()
    import chematic
    fixture = ROOT / 'validation/moljson_graph_input/cases.json'
    packet = strict_json(fixture.read_text(encoding='utf-8'))
    cases = packet['valid'] + packet['invalid']
    node_command = ['node', '-e', NODE, str(args.wasm_module.resolve())]
    node_process = bounded_run(node_command, json.dumps(cases, ensure_ascii=False).encode('utf-8'))
    require(node_process['resource_status'] == 'completed' and node_process['exit_code'] == 0, 'Node worker failed')
    node_results = strict_json(node_process['stdout'])
    require(len(node_results) == len(cases), 'Node record count')
    records = []
    for group in ['valid', 'invalid']:
        for case in packet[group]:
            record = dict(case, group=group, **node_results[len(records)])
            for lane, run in [('python_direct', lambda: chematic.from_moljson(case['input']).to_moljson()),
                              ('python_common', lambda: chematic.convert_format(case['input'], 'moljson', 'moljson'))]:
                try:
                    output = run()
                    record[lane] = dict(status='accepted', output=output)
                except ValueError as e:
                    record[lane] = dict(status='refused' if str(e).startswith('MolJSON:') else 'error', error=str(e))
            command = [str(args.cli.resolve()), 'convert', '--input-format', 'moljson', '--output-format', 'moljson']
            process = bounded_run(command, case['input'].encode('utf-8'))
            process = {k: v.decode('utf-8', errors='replace') if isinstance(v, bytes) else v for k, v in process.items()}
            status = cli_status(process)
            record['cli'] = dict(status=status, output=process['stdout'], command=command, process=process)
            if group == 'valid':
                for lane in ['python_direct', 'python_common', 'wasm_direct', 'wasm_common', 'cli']:
                    record[lane]['observation'] = observe(record[lane].get('output', ''), case['expected'])
            records.append(record)
    summary = {}
    for lane in ['python_direct', 'python_common', 'wasm_direct', 'wasm_common', 'cli']:
        pos = [r[lane] for r in records if r['group'] == 'valid']
        neg = [r[lane] for r in records if r['group'] == 'invalid']
        summary[lane] = dict(preserved=sum(r['observation']['preserved'] for r in pos), valid_count=len(pos),
                             explicit_refusals=sum(r['status'] == 'refused' for r in neg), refusal_count=len(neg))
    source = 'crates/chematic-mol/src/moljson.rs'
    parent = subprocess.check_output(['git', 'rev-parse', args.parent_ref], cwd=ROOT, text=True).strip()
    if args.baseline_ref:
        ref = subprocess.check_output(['git', 'rev-parse', args.baseline_ref], cwd=ROOT, text=True).strip()
        data = subprocess.check_output(['git', 'show', f'{ref}:{source}'], cwd=ROOT)
        require(data == subprocess.check_output(['git', 'show', f'{parent}:{source}'], cwd=ROOT), 'baseline native MolJSON source drift')
        sources = dict(reference=ref, matches_parent_source=True, files=[dict(path=source, sha256=hashlib.sha256(data).hexdigest())])
    else:
        sources = dict(reference='working tree', files=[identity(ROOT/source), identity(ROOT/'crates/chematic-mol/src/moljson/fields.rs'), identity(ROOT/'crates/chematic-mol/src/cjson_record.rs')])
    obabel = Path('/opt/homebrew/bin/obabel')
    inventory = subprocess.check_output([str(obabel), '-L', 'formats'], text=True)
    require('moljson' not in inventory.lower(), 'Open Babel inventory now advertises MolJSON; update comparison scope')
    result = dict(schema='chematic.moljson-graph-input-observation.v1', parent_commit=parent,
                  scope='Original known-field local MolJSON graph packet. Refusals separate from preservation. Not oxpig schema, source metadata, stereo, general writer, speed or toolkit superiority evidence.',
                  baseline=bool(args.baseline_ref), fixture=identity(fixture), native_sources=sources,
                  python_runtime=identity(importlib.import_module('chematic.chematic').__file__),
                  python_version=sys.version, cli=identity(args.cli), harness=identity(Path(__file__)),
                  wasm_module=identity(args.wasm_module),
                  wasm_binary=identity(args.wasm_module.with_name('chematic_wasm_bg.wasm')),
                  node_command=node_command, node_version=subprocess.check_output(['node','--version'],text=True).strip(),
                  node_process={k:v.decode('utf-8',errors='replace') if isinstance(v,bytes) else v for k,v in node_process.items()},
                  process_runner=identity(ROOT/'scripts/checks/check_interchange_preservation.py'),
                  openbabel=dict(binary=identity(obabel), inventory=inventory, same_format_score=None),
                  records=records, summary=summary)
    args.output.write_text(json.dumps(result, ensure_ascii=False, indent=2, allow_nan=False)+'\n', encoding='utf-8')
    print(json.dumps(summary))
    if not args.baseline_ref:
        require(all(r['preserved'] == r['valid_count'] and r['explicit_refusals'] == r['refusal_count'] for r in summary.values()), 'candidate contract failure')


if __name__ == '__main__':
    main()
