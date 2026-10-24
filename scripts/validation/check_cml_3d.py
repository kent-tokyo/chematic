#!/usr/bin/env python3
"""Bounded CML-to-CJSON source-order Cartesian development comparison."""
from __future__ import annotations

import argparse
import hashlib
import importlib
import json
import math
from pathlib import Path
import struct
import sys

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
from scripts.checks.check_interchange_preservation import bounded_run, identity, require, strict_json

NODE = """
const fs=require('node:fs'),w=require(process.argv[1]);
const cases=JSON.parse(fs.readFileSync(0,'utf8'));
process.stdout.write(JSON.stringify(cases.map(c=>{try{return {id:c.id,output:w.convert_common_format(c.input,'cml','cjson')}}catch(e){return {id:c.id,error:String(e)}}})));
"""


def layouts(text):
    return [('lf', text), ('crlf', text.replace('\n', '\r\n')),
            ('compact', ''.join(line.strip() for line in text.split('\n')))]


def coordinate_bits(values):
    require(isinstance(values, list), 'coordinate array missing')
    require(all(isinstance(v, (int, float)) and not isinstance(v, bool) for v in values),
            'invalid coordinate type')
    require(all(math.isfinite(v) for v in values), 'nonfinite coordinates')
    return [struct.pack('>d', value).hex() for value in values]


def graph_matches(data, case):
    atoms = data.get('atoms', {})
    bonds = data.get('bonds', {})
    indices = bonds.get('connections', {}).get('index', [])
    orders = bonds.get('order', [])
    if len(indices) != len(orders) * 2:
        return False
    def edges(ids, values):
        return sorted((min(a,b), max(a,b), value) for a,b,value in zip(ids[::2], ids[1::2], values))
    return (atoms.get('elements', {}).get('number') == case['numbers']
            and edges(indices, orders) == edges(case['bond_indices'], case['bond_orders']))


def observe(case, result):
    text = result.get('output')
    if text is None or not text.strip():
        return {'outcome': 'no_output', 'error_status': 'error' in result or result.get('exit_code',0) != 0,
                'coordinate_bits': False, 'selected_graph': False, 'absence_retained': False}
    try:
        data = strict_json(text)
        require(isinstance(data, dict), 'output is not a JSON object')
        coords = data.get('atoms', {}).get('coords', {}).get('3d')
        if case['group'] == 'invalid':
            return {'outcome': 'accepted_outside_profile', 'coordinate_bits': False,
                    'selected_graph': False, 'absence_retained': False}
        flat = [v for row in case['coords'] or [] for v in row]
        return {'outcome': 'json_output', 'coordinate_bits': coords is not None and coordinate_bits(coords) == coordinate_bits(flat),
                'selected_graph': graph_matches(data, case), 'absence_retained': coords is None}
    except (ValueError, KeyError, TypeError, AttributeError, OverflowError) as error:
        return {'outcome': 'invalid_output', 'reason': str(error), 'coordinate_bits': False,
                'selected_graph': False, 'absence_retained': False}


def terminal(raw):
    return {k:v.decode('utf-8', errors='replace') if isinstance(v, bytes) else v for k,v in raw.items()}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--wasm-module', type=Path, required=True)
    parser.add_argument('--obabel', type=Path, default=Path('/opt/homebrew/bin/obabel'))
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--baseline', action='store_true')
    args = parser.parse_args()
    import chematic
    extension = importlib.import_module('chematic.chematic')
    packet_path = ROOT/'validation/cml_3d/cases.json'
    packet = strict_json(packet_path.read_text())
    cases = []
    for group in ['valid', 'absent', 'invalid']:
        for case in packet[group]:
            require(hashlib.sha256(case['input'].encode()).hexdigest() == case['sha256'], 'fixture SHA mismatch')
            for layout, text in layouts(case['input']):
                cases.append({**case, 'group': group, 'id':case['id']+'/'+layout, 'input':text})
    command = ['node', '-e', NODE, str(args.wasm_module.resolve())]
    node_input = json.dumps(cases).encode()
    raw = bounded_run(command, node_input)
    require(raw['resource_status']=='completed' and raw['exit_code']==0, 'Node failed')
    nodes = strict_json(raw['stdout'])
    require(len(nodes)==len(cases), 'Node count mismatch')
    results = []
    for case, node in zip(cases,nodes):
        require(case['id']==node['id'], 'Node order mismatch')
        row = {'id':case['id'], 'group':case['group'], 'input_sha256':hashlib.sha256(case['input'].encode()).hexdigest(), 'wasm':node}
        try:
            row['python'] = {'output':chematic.convert_format(case['input'],'cml','cjson')}
        except ValueError as error:
            row['python'] = {'error':str(error)}
        command = [str(args.obabel.resolve()), '-icml', '-ocjson']
        result = terminal(bounded_run(command, case['input'].encode()))
        result['command'] = command
        if result['resource_status']=='completed':
            result['output'] = result['stdout']
        row['openbabel'] = result
        for label in ['python','wasm','openbabel']:
            row[label]['observation'] = observe(case,row[label])
        if not args.baseline:
            for label in ['python','wasm']:
                observation = row[label]['observation']
                if case['group']=='invalid' or case['id'].split('/')[0] in {'both_2d_and_3d','only_2d'}:
                    require(observation['outcome']=='no_output' and observation['error_status'], 'outside-profile input not refused')
                else:
                    require(observation['selected_graph'], 'selected graph changed')
                    metric = 'coordinate_bits' if case['group']=='valid' else 'absence_retained'
                    require(observation[metric], 'geometry changed')
        results.append(row)
    summary = {}
    for label in ['python','wasm','openbabel']:
        summary[label] = {}
        for group in ['valid','absent','invalid']:
            rows = [row[label]['observation'] for row in results if row['group']==group]
            summary[label][group] = {'cases':len(rows), **{metric:sum(bool(row.get(metric)) for row in rows) for metric in ['coordinate_bits','selected_graph','absence_retained','error_status']},
                                    'outcomes':{outcome:sum(row['outcome']==outcome for row in rows) for outcome in sorted({row['outcome'] for row in rows})}}
    record = {'schema':'chematic.cml-3d-development-results.v1', 'baseline':args.baseline,
              'packet':identity(packet_path), 'python_package':identity(chematic.__file__),
              'python_extension':identity(extension.__file__), 'wasm_module':identity(args.wasm_module),
              'wasm_binary':identity(args.wasm_module.with_name('chematic_wasm_bg.wasm')),
              'openbabel':identity(args.obabel), 'openbabel_version':terminal(bounded_run([str(args.obabel.resolve()),'-V'])),
              'node_run':{**terminal(raw),'command':['node','-e','<NODE source in hashed runner>',str(args.wasm_module.resolve())], 'input_sha256':hashlib.sha256(node_input).hexdigest()},
              'summary':summary, 'cases':results}
    args.output.write_text(json.dumps(record,indent=2)+'\n')
    print(json.dumps(summary))


if __name__=='__main__':
    main()
