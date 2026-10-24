#!/usr/bin/env python3
"""Compare ordinary CML self-conversion on original source-retention fixtures."""
from __future__ import annotations
import argparse
import gzip
import importlib
import json
from pathlib import Path
import sys
import xml.etree.ElementTree as ET

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
from scripts.checks.check_interchange_preservation import bounded_run, identity, require, strict_json
from scripts.validation.check_cml_source_conversion import snapshot, compare

NODE = """
const fs=require('node:fs'),w=require(process.argv[1]);
process.stdout.write(JSON.stringify(JSON.parse(fs.readFileSync(0,'utf8')).map(c=>{
try{return {id:c.id,output:w.convert_common_format(c.input,'cml','cml')}}
catch(e){return {id:c.id,error:String(e)}}})));
"""


def terminal(result):
    return {k: v.decode('utf8', errors='replace') if isinstance(v, bytes) else v for k,v in result.items()}


def observe(source, result):
    output = result.get('output', '')
    result['exact_source_bytes'] = output == source
    if not output:
        result['outcome'] = 'no_output'
        return
    try:
        ET.fromstring(output)
    except ET.ParseError as error:
        result['outcome'] = 'unresolved_output'
        result['observer_error'] = str(error)
        return
    try:
        result['selected_fields'] = compare(snapshot(source), snapshot(output))
        result['outcome'] = 'selected_fields_retained' if result['selected_fields']['all_selected_fields_retained'] else 'selected_field_loss'
    except (ValueError, KeyError, TypeError, IndexError) as error:
        result['outcome'] = 'selected_graph_mismatch'
        result['observer_error'] = str(error)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--cli', type=Path, required=True)
    parser.add_argument('--wasm-module', type=Path, required=True)
    parser.add_argument('--obabel', type=Path, default=Path('/opt/homebrew/bin/obabel'))
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--baseline', action='store_true')
    parser.add_argument('--build-label', default='candidate working tree')
    args = parser.parse_args()
    import chematic
    extension = importlib.import_module('chematic.chematic')
    packet_path = ROOT/'validation/cml_source_preservation/cases.json'
    cases = []
    for case in strict_json(packet_path.read_text())['cases']:
        lf = case['source'].replace('\r\n', '\n')
        for label, text in [('lf', lf), ('crlf', lf.replace('\n', '\r\n')), ('no-final-newline', lf.rstrip('\n'))]:
            cases.append(dict(id=case['id']+'/'+label, input=text))
    raw = terminal(bounded_run(['node', '-e', NODE, str(args.wasm_module.resolve())], json.dumps(cases).encode()))
    require(raw['resource_status'] == 'completed' and raw['exit_code'] == 0, 'Node runner failed')
    node = strict_json(raw['stdout'])
    require(len(node) == len(cases), 'Node row count mismatch')
    rows = []
    for case, wasm in zip(cases, node):
        require(case['id'] == wasm['id'], 'Node row order mismatch')
        row = dict(case=case, wasm=wasm)
        try:
            row['python'] = dict(output=chematic.convert_format(case['input'], 'cml', 'cml'))
        except ValueError as error:
            row['python'] = dict(error=str(error))
        for label, command in [('cli', [str(args.cli.resolve()), 'convert', '--input-format', 'cml', '--output-format', 'cml']),
                               ('openbabel', [str(args.obabel.resolve()), '-icml', '-ocml']),
                               ('openbabel_metadata', [str(args.obabel.resolve()), '-icml', '-ocml', '-xp', '-xm'])]:
            result = terminal(bounded_run(command, case['input'].encode()))
            result['command'] = command
            if result['resource_status'] == 'completed' and result['exit_code'] == 0:
                result['output'] = result['stdout']
            row[label] = result
        for label in ['python', 'wasm', 'cli', 'openbabel', 'openbabel_metadata']:
            observe(case['input'], row[label])
        rows.append(row)
    labels = ['python', 'wasm', 'cli', 'openbabel', 'openbabel_metadata']
    summary = {label: dict(exact_source_bytes=sum(r[label]['exact_source_bytes'] for r in rows),
                          **{outcome: sum(r[label]['outcome'] == outcome for r in rows) for outcome in
                             ['selected_fields_retained', 'selected_field_loss', 'selected_graph_mismatch', 'no_output', 'unresolved_output']}) for label in labels}
    sources = ['crates/chematic-mol/src/conversion/cml_source.rs', 'crates/chematic-py/src/formats/cml.rs',
               'crates/chematic-py/python/chematic/__init__.py', 'crates/chematic-wasm/src/mol_io/common_formats.rs',
               'crates/chematic-cli/src/conversion_format.rs', 'crates/chematic-cli/src/main.rs']
    record = dict(schema='chematic.cml-common-retention.v1', baseline=args.baseline, build_label=args.build_label,
                  implementation_sources=[] if args.baseline else [identity(ROOT/p) for p in sources],
                  scope='CML -> CML common bridge; original six-document development packet times three newline forms; opaque retention, not chemical/stereo validation',
                  packet=identity(packet_path), runner=identity(__file__), observer=identity(ROOT/'scripts/validation/check_cml_source_conversion.py'),
                  python_package=identity(chematic.__file__), python_extension=identity(extension.__file__),
                  wasm_module=identity(args.wasm_module), wasm_binary=identity(args.wasm_module.with_name('chematic_wasm_bg.wasm')),
                  cli=identity(args.cli), openbabel=identity(args.obabel),
                  openbabel_version=terminal(bounded_run([str(args.obabel.resolve()), '-V'])),
                  node_terminal=raw, summary=summary, cases=rows)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_bytes(gzip.compress((json.dumps(record, indent=2, ensure_ascii=False)+'\n').encode(), mtime=0))
    print(json.dumps(summary))
    if not args.baseline:
        require(all(row[label]['exact_source_bytes'] and row[label]['outcome'] == 'selected_fields_retained'
                    for row in rows for label in labels[:3]), 'candidate retention failed')


if __name__ == '__main__':
    main()
