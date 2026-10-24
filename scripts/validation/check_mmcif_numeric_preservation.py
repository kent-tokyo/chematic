#!/usr/bin/env python3
"""Retain paired mmCIF numeric outputs against an authored binary64 oracle."""
from __future__ import annotations
import argparse
import importlib
import json
from pathlib import Path
import subprocess
import sys

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
from scripts.checks.check_interchange_preservation import bounded_run, identity, require, strict_json
from scripts.validation.mmcif_numeric_observer import observe

NODE = r"""
const fs=require('node:fs'),w=require(process.argv[1]);
process.stdout.write(JSON.stringify(JSON.parse(fs.readFileSync(0,'utf8')).map(c=>{
 try{const d=JSON.parse(w.mmcif_to_json(c.input));
 return {id:c.id,output:w.write_mmcif_json(JSON.stringify(d.atoms),JSON.stringify(d.cell),d.space_group ?? '',c.id)}}
 catch(e){return {id:c.id,error:String(e)}}
})));
"""


def terminal(raw):
    return {key: value.decode('utf-8', errors='replace') if isinstance(value, bytes) else value for key, value in raw.items()}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--wasm-module', type=Path, required=True)
    parser.add_argument('--obabel', type=Path, default=Path('/opt/homebrew/bin/obabel'))
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--baseline', action='store_true')
    args = parser.parse_args()
    import chematic
    extension = importlib.import_module('chematic.chematic')
    fixture = ROOT / 'validation/mmcif_numeric_preservation/cases.json'
    cases = strict_json(fixture.read_text())['cases']
    for case in cases:
        require(observe(case['expected'], case['input'])['preserved'], 'source oracle mismatch: ' + case['id'])
    raw = terminal(bounded_run(['node', '-e', NODE, str(args.wasm_module.resolve())], json.dumps(cases).encode()))
    require(raw['resource_status'] == 'completed' and raw['exit_code'] == 0, 'Node runner failed')
    wasm = strict_json(raw['stdout'])
    require(len(wasm) == len(cases), 'Node row count mismatch')
    command = [str(args.obabel.resolve()), '-immcif', '-ommcif', '-ab']
    results = []
    for case, node in zip(cases, wasm, strict=True):
        require(node.pop('id') == case['id'], 'Node order mismatch')
        row = dict(id=case['id'], group=case['group'], input=case['input'], expected=case['expected'], wasm=node)
        try:
            parsed = chematic.parse_mmcif(case['input'])
            output = chematic.write_mmcif(parsed['atoms'], parsed['cell'], parsed['space_group'], case['id'])
            row['python'] = dict(output=output)
        except ValueError as error:
            row['python'] = dict(error=str(error))
        value = terminal(bounded_run(command, case['input'].encode(), timeout=60))
        if value['resource_status'] == 'completed' and value['exit_code'] == 0:
            value['output'] = value['stdout']
        row['openbabel'] = value
        for route in ('python', 'wasm', 'openbabel'):
            row[route]['observation'] = observe(case['expected'], row[route].get('output'))
        results.append(row)
    summary = {}
    for route in ('python', 'wasm', 'openbabel'):
        summary[route] = {}
        for group in ('control', 'reproduction', 'storage'):
            observations = [row[route]['observation'] for row in results if row['group'] == group]
            summary[route][group] = dict(count=len(observations), all_numeric=sum(o['preserved'] for o in observations),
                no_output=sum(o['outcome'] == 'no_output' for o in observations),
                channels={key: dict(eligible=sum(o['fields'][key] is not None for o in observations),
                                    preserved=sum(o['fields'][key] is True for o in observations))
                          for key in ('atom_identity', 'cartesian', 'occupancy', 'b_iso', 'cell')})
    baseline_commit = '342d1d5f9ad09aa6ade4e3079a902cba856209e2'
    comparison_base = 'fcbf6d0a8c901dee5220cca965cf119c0a3b7b46'
    paths = ['crates/chematic-mol/src/mmcif.rs', 'crates/chematic-py/src/formats/materials.rs',
             'crates/chematic-wasm/src/format_io.rs']
    if args.baseline:
        subprocess.run(['git', 'diff', '--exit-code', baseline_commit, comparison_base, '--', *paths], cwd=ROOT, check=True)
    record = dict(schema='chematic.mmcif-numeric-preservation-report.v1', baseline=args.baseline,
        baseline_runtime_commit=baseline_commit if args.baseline else None, reader_comparison_base=comparison_base,
        head=subprocess.check_output(['git', 'rev-parse', 'HEAD'], cwd=ROOT, text=True).strip(),
        scope='authored development packet; parsed finite binary64 numeric values, ordered atom id/element; no uncertainty/source bytes/full namespace/schema/bond/stereo/speed claim; storage group is arithmetic only',
        sources=[] if args.baseline else [identity(ROOT / path) for path in paths],
        baseline_sources=[identity(ROOT / path, subprocess.check_output(['git', 'show', baseline_commit + ':' + path], cwd=ROOT))
                          for path in paths] if args.baseline else [],
        runner=identity(__file__), observer=identity(ROOT / 'scripts/validation/mmcif_numeric_observer.py'),
        process_bounds=identity(ROOT / 'scripts/checks/check_interchange_preservation.py'), fixture=identity(fixture),
        python_version=sys.version, python_package=identity(chematic.__file__), python_extension=identity(extension.__file__),
        node_version=terminal(bounded_run(['node', '--version'])), wasm_module=identity(args.wasm_module),
        wasm_binary=identity(args.wasm_module.with_name('chematic_wasm_bg.wasm')),
        openbabel_binary=identity(args.obabel), openbabel_command=command,
        process_limits=dict(node_seconds=10, openbabel_seconds_per_case=60, combined_output_bytes=1048576),
        openbabel_version=terminal(bounded_run([str(args.obabel.resolve()), '-V'])), summary=summary, cases=results)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(record, indent=2, ensure_ascii=False, allow_nan=False) + '\n')
    print(json.dumps(summary))
    if not args.baseline:
        for route in ('python', 'wasm'):
            require(all(row[route]['observation']['preserved'] for row in results), route + ' failed')


if __name__ == '__main__':
    main()
