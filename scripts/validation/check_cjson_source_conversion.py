#!/usr/bin/env python3
"""Paired CJSON same-format source/JSON-field retention development diagnostic.

Use a built package on PYTHONPATH and a matching Node/WASM module. This is an
immutable source workflow, not editable/cross-format chemistry verification.
"""
from __future__ import annotations

import argparse
import hashlib
import importlib
import json
import math
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
from scripts.check_interchange_preservation import bounded_run, identity, strict_json, require

NODE = """
const fs=require('node:fs'); const wasm=require(process.argv[1]);
const cases=JSON.parse(fs.readFileSync(0,'utf8'));
process.stdout.write(JSON.stringify(cases.map(c=>{
 try {return {id:c.id,output:wasm.convert_common_format(c.input,'cjson','cjson')};}
 catch(error) {return {id:c.id,error:String(error)};}
})));
"""


def source_fields_equal(source, candidate):
    """Keep source object fields; additional candidate object fields are allowed."""
    if isinstance(source, dict):
        return isinstance(candidate, dict) and all(
            k in candidate and source_fields_equal(v, candidate[k]) for k, v in source.items())
    if isinstance(source, list):
        return (isinstance(candidate, list) and len(source) == len(candidate)
                and all(source_fields_equal(a, b) for a, b in zip(source, candidate)))
    if isinstance(source, bool) or source is None or isinstance(source, str):
        return type(source) is type(candidate) and source == candidate
    if type(candidate) not in (int, float) or source != candidate:
        return False
    return source != 0 or math.copysign(1, source) == math.copysign(1, candidate)


def changed_source_paths(source, candidate, path=''):
    if source_fields_equal(source, candidate):
        return []
    if isinstance(source, dict) and isinstance(candidate, dict):
        paths = []
        for key, value in source.items():
            child = path + '/' + key.replace('~', '~0').replace('/', '~1')
            paths.extend(changed_source_paths(value, candidate[key], child) if key in candidate else [child])
        return paths
    return [path or '/']


def variants(packet):
    for case in packet['cases']:
        require(hashlib.sha256(case['input'].encode()).hexdigest() == case['sha256'], 'fixture hash mismatch')
        for layout in packet['variants']:
            text = case['input']
            if layout == 'crlf':
                text = text.replace('\n', '\r\n')
            elif layout == 'minified':
                text = json.dumps(strict_json(text), ensure_ascii=True, separators=(',', ':'))
            yield dict(id=case['id'] + '/' + layout, input=text)


def observe(text, output):
    if output is None:
        return dict(exact_source=False, source_fields_retained=False)
    expected = strict_json(text)
    try:
        value = strict_json(output)
    except (ValueError, RuntimeError):
        return dict(exact_source=False, source_fields_retained=False, invalid_json_output=True)
    return dict(exact_source=output == text,
                source_fields_retained=source_fields_equal(expected, value),
                changed_source_paths=changed_source_paths(expected, value))


def terminal(raw):
    return {k: v.decode('utf-8', errors='replace') if isinstance(v, bytes) else v for k, v in raw.items()}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--wasm-module', type=Path, required=True)
    parser.add_argument('--obabel', type=Path, default=Path('/opt/homebrew/bin/obabel'))
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--baseline', action='store_true')
    args = parser.parse_args()
    import chematic
    extension = importlib.import_module('chematic.chematic')
    packet_path = ROOT / 'validation/cjson_source_conversion/cases.json'
    cases = list(variants(strict_json(packet_path.read_text(encoding='utf-8'))))
    raw = bounded_run(['node', '-e', NODE, str(args.wasm_module.resolve())], json.dumps(cases).encode())
    require(raw['resource_status'] == 'completed' and raw['exit_code'] == 0, 'WASM process failed')
    wasm = strict_json(raw['stdout'])
    require(len(wasm) == len(cases), 'WASM result count mismatch')
    results = []
    for case, wasm_row in zip(cases, wasm):
        require(case['id'] == wasm_row['id'], 'WASM order mismatch')
        row = dict(id=case['id'], input_sha256=hashlib.sha256(case['input'].encode()).hexdigest())
        try:
            row['python'] = dict(output=chematic.convert_format(case['input'], 'cjson', 'cjson'))
        except ValueError as error:
            row['python'] = dict(error=str(error))
        row['wasm'] = wasm_row
        for engine in ['python', 'wasm']:
            row[engine]['observation'] = observe(case['input'], row[engine].get('output'))
        for label, options in [('openbabel_default', []), ('openbabel_minified', ['-xm']), ('openbabel_verbose', ['-xv'])]:
            command = [str(args.obabel.resolve()), '-icjson', '-ocjson'] + options
            run = terminal(bounded_run(command, case['input'].encode()))
            run['command'] = command
            output = run['stdout'] if run['resource_status'] == 'completed' and run['exit_code'] == 0 else None
            run['observation'] = observe(case['input'], output)
            row[label] = run
        results.append(row)
    engines = ['python', 'wasm', 'openbabel_default', 'openbabel_minified', 'openbabel_verbose']
    summary = {e: {metric: sum(r[e]['observation'][metric] for r in results)
                    for metric in ['exact_source', 'source_fields_retained']} for e in engines}
    result = dict(schema='chematic.cjson-source-conversion-validation.v1', baseline=args.baseline,
        packet=identity(packet_path), python_package=identity(chematic.__file__), python_extension=identity(extension.__file__),
        wasm_module=identity(args.wasm_module), wasm_binary=identity(args.wasm_module.with_name('chematic_wasm_bg.wasm')),
        openbabel=identity(args.obabel), openbabel_version=terminal(bounded_run([str(args.obabel.resolve()), '-V'])),
        node_process=terminal(raw), case_count=len(results), summary=summary, cases=results)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, ensure_ascii=False, indent=2) + '\n')
    print(json.dumps(summary))
    if not args.baseline:
        require(all(r[e]['observation']['exact_source'] and r[e]['observation']['source_fields_retained']
                    for r in results for e in ['python', 'wasm']), 'candidate retention failed')


if __name__ == '__main__':
    main()
