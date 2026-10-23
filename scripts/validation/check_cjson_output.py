#!/usr/bin/env python3
"""Bounded CML -> CJSON graph projection/refusal comparison; no toolkit ranking."""
from __future__ import annotations
import argparse
import gzip
import hashlib
import importlib
import json
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
from scripts.checks.check_interchange_preservation import bounded_run, identity, require, strict_json
from scripts.validation.check_v3000_coordinate_projection import observe, terminal

NODE = """
const fs=require('node:fs'),w=require(process.argv[1]);
process.stdout.write(JSON.stringify(JSON.parse(fs.readFileSync(0,'utf8')).map(c=>{
try{return {id:c.id,output:w.convert_common_format(c.input,'cml','cjson')}}
catch(e){return {id:c.id,error:String(e)}}})));
"""


def classify(case, result):
    output = result.get('output', '')
    try:
        value = strict_json(output)
        valid = isinstance(value, dict) and isinstance(value['atoms']['elements']['number'], list)
    except (ValueError, KeyError, TypeError):
        valid = False
    result['molecular_json'] = valid
    if case['accepted']:
        result['observation'] = observe(case['expected'], output)
        result['outcome'] = 'preserved' if all(result['observation'].values()) else 'changed' if valid else 'no_output'
    else:
        result['outcome'] = 'unsupported_output' if valid else 'refused' if result.get('error') else 'no_output'
        if valid:
            result['emitted_bond_orders'] = value.get('bonds', {}).get('order')
        result['typed_refusal'] = not valid and 'cannot preserve bond order ' + case['error_order'] in result.get('error', '')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--wasm-module', type=Path, required=True)
    parser.add_argument('--obabel', type=Path, default=Path('/opt/homebrew/bin/obabel'))
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    import chematic
    extension = importlib.import_module('chematic.chematic')
    packet_path = ROOT / 'validation/cjson_output/cases.json'
    cases = strict_json(packet_path.read_text())['cases']
    for case in cases:
        require(hashlib.sha256(case['input'].encode()).hexdigest() == case['sha256'], 'fixture hash mismatch')
    node_raw = terminal(bounded_run(['node', '-e', NODE, str(args.wasm_module.resolve())], json.dumps(cases).encode()))
    require(node_raw['resource_status'] == 'completed' and node_raw['exit_code'] == 0, 'Node failed')
    node = strict_json(node_raw['stdout'])
    require(len(node) == len(cases), 'Node count mismatch')
    rows = []
    for case, wasm in zip(cases, node):
        require(wasm['id'] == case['id'], 'Node order mismatch')
        row = dict(case=case, wasm=wasm)
        mol, coords = chematic.from_cml_with_3d(case['input'])
        row['legacy_python'] = dict(output=mol.to_cjson(coords or []))
        try:
            row['python'] = dict(output=chematic.convert_format(case['input'], 'cml', 'cjson'))
        except ValueError as error:
            row['python'] = dict(error=str(error))
        for label, flags in [('openbabel', []), ('openbabel_single', ['-xs'])]:
            command = [str(args.obabel.resolve()), '-icml', '-ocjson'] + flags
            result = terminal(bounded_run(command, case['input'].encode()))
            result['command'] = command
            if result['resource_status'] == 'completed' and result['exit_code'] == 0:
                result['output'] = result['stdout']
            row[label] = result
        for label in ['python', 'wasm', 'legacy_python', 'openbabel', 'openbabel_single']:
            classify(case, row[label])
        rows.append(row)
    labels = ['python', 'wasm', 'legacy_python', 'openbabel', 'openbabel_single']
    summary = {label: {outcome: sum(r[label]['outcome'] == outcome for r in rows)
                      for outcome in ['preserved', 'changed', 'refused', 'unsupported_output', 'no_output']}
               for label in labels}
    sources = [Path(__file__), ROOT/'crates/chematic-mol/src/cjson.rs', ROOT/'crates/chematic-mol/src/cjson_output.rs',
               ROOT/'crates/chematic-py/src/mol_methods/fingerprints.rs', ROOT/'crates/chematic-py/python/chematic/__init__.py',
               ROOT/'crates/chematic-wasm/src/mol_io/common_formats.rs']
    record = dict(schema='chematic.cjson-output-results.v1', packet=identity(packet_path),
                  sources=[identity(p) for p in sources], python_package=identity(chematic.__file__),
                  python_extension=identity(extension.__file__), wasm_module=identity(args.wasm_module),
                  wasm_binary=identity(args.wasm_module.with_name('chematic_wasm_bg.wasm')),
                  openbabel=identity(args.obabel), openbabel_version=terminal(bounded_run([str(args.obabel.resolve()), '-V'])),
                  node_terminal=node_raw, summary=summary, cases=rows)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_bytes(gzip.compress((json.dumps(record, indent=2) + '\n').encode(), mtime=0))
    print(json.dumps(summary))
    for row in rows:
        for label in ['python', 'wasm']:
            require(row[label]['outcome'] == ('preserved' if row['case']['accepted'] else 'refused'), 'candidate failed')
            if not row['case']['accepted']:
                require(row[label]['typed_refusal'], 'wrong candidate refusal')


if __name__ == '__main__':
    main()
