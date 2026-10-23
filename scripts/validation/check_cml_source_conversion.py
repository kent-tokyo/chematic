#!/usr/bin/env python3
"""Paired source-built CML retention and binding parity on original inputs."""
from __future__ import annotations

import argparse
from collections import Counter
import json
import math
from pathlib import Path
import platform
import sys
import xml.etree.ElementTree as ET

from scripts.checks.check_interchange_preservation import bounded_run, identity, require

ROOT = Path(__file__).resolve().parents[2]
NS = '{http://www.xml-cml.org/schema}'
AXES = {'x2', 'y2', 'x3', 'y3', 'z3'}
TOLERANCE = 1e-12
PYTHON_CALL = 'import chematic,sys; print(chematic.convert_with_report_json(sys.stdin.buffer.read().decode("utf-8"), "cml", "cml", sys.argv[1]))'
NODE_CALL = '''import {readFileSync} from 'node:fs';
import {pathToFileURL} from 'node:url';
const wasm = await import(pathToFileURL(process.argv[1]).href);
console.log(wasm.convert_with_report_json(readFileSync(0, 'utf8'), 'cml', 'cml', process.argv[2]));'''


def tree(node):
    # Expanded names remove prefix spelling; insignificant indentation is ignored.
    content = lambda value: value if value and value.strip() else ''
    return [node.tag, sorted(node.attrib.items()), content(node.text),
            content(node.tail), [tree(child) for child in node]]


def snapshot(text):
    root = ET.fromstring(text)
    require(root.tag == NS+'molecule', 'expected one CML molecule')
    atoms, ids, bonds = {}, {}, []
    for node in root.iter(NS+'atom'):
        attrs = node.attrib.copy()
        element, atom_id = attrs['elementType'], attrs['id']
        require(element in ('C', 'O') and element not in atoms and atom_id not in ids, 'unexpected atom identity')
        require(int(attrs.get('formalCharge', 0)) == 0, 'formal charge changed')
        require(int(attrs.get('isotopeNumber', attrs.get('isotope', 0))) == 0, 'isotope changed')
        coords = {k: float(attrs.pop(k)) for k in list(attrs) if k in AXES}
        require(all(math.isfinite(v) for v in coords.values()), 'nonfinite coordinate')
        atoms[element], ids[atom_id] = {'attributes': attrs, 'coordinates': coords}, element
    for node in root.iter(NS+'bond'):
        a, b = node.attrib['atomRefs2'].split()
        require(node.attrib['order'] in ('S', '1'), 'bond order changed')
        bonds.append(sorted([ids[a], ids[b]]))
    require(set(atoms) == {'C', 'O'} and bonds == [['C', 'O']], 'graph changed')
    return {'atoms': atoms, 'payload': [tree(n) for n in root if n.tag not in (NS+'atomArray', NS+'bondArray')]}


def compare(expected, actual):
    fields, deltas = {}, []
    for element, source in expected['atoms'].items():
        candidate = actual['atoms'][element]
        fields[element+'_source_attributes'] = all(candidate['attributes'].get(k) == v for k, v in source['attributes'].items())
        fields[element+'_coordinate_presence'] = source['coordinates'].keys() == candidate['coordinates'].keys()
        for axis, value in source['coordinates'].items():
            after = candidate['coordinates'].get(axis)
            delta = abs(value-after) if after is not None else None
            fields[element+'_'+axis] = delta is not None and delta <= TOLERANCE
            if delta is not None:
                deltas.append(delta)
    # Score source objects in their order. Added metadata outside those objects
    # must not turn an input with no property payload into a preservation loss.
    remaining = iter(actual['payload'])
    fields['property_extension_payload'] = all(
        any(candidate == source for candidate in remaining) for source in expected['payload'])
    return {'fields': fields, 'max_comparable_coordinate_error_angstrom': max(deltas, default=None),
            'all_selected_fields_retained': all(fields.values())}


def check_report(report, source, policy):
    require(report['schema'] == 'chematic.conversion-report.v1', 'wrong report schema')
    require(report['policy'] == policy and report['status'] == 'completed', 'unexpected policy/status')
    require(report['verification_scope'] == 'exact_source_bytes' and report['candidate_graph_verified'] is False, 'wrong verification scope')
    require(report['candidate_atom_mapping'] == report['candidate_bond_mapping'] == [], 'unexpected graph mapping')
    require({o['field']: o['status'] for o in report['observations']} == {
        'source_text': 'preserved', 'source_semantics': 'opaque_retained'}, 'wrong observations')
    require(report['output'] == source, 'source bytes changed')
    return report['output']


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--cli', type=Path, required=True)
    parser.add_argument('--python', default=sys.executable)
    parser.add_argument('--wasm-module', type=Path, required=True)
    parser.add_argument('--obabel', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    require(not args.output.exists(), 'output already exists')
    manifest = ROOT/'validation/cml_source_preservation/cases.json'
    packet = json.loads(manifest.read_text())
    version = bounded_run([str(args.obabel), '-V'])
    require(version['exit_code'] == 0 and b'Open Babel 3.2.1' in version['stdout'], 'expected Open Babel 3.2.1')
    rows, parity = [], []
    for case in packet['cases']:
        source, reports = case['source'], {}
        expected = snapshot(source)
        require(len(source.encode()) < 16384, 'packet input too large')
        commands = []
        for policy in ['strict', 'allow_loss', 'preserve_source']:
            commands += [
                ('cli', policy, [str(args.cli.resolve()), 'convert-with-report', '--input-format', 'cml', '--output-format', 'cml', '--policy', policy]),
                ('python', policy, [args.python, '-c', PYTHON_CALL, policy]),
                ('wasm', policy, ['node', '--input-type=module', '-e', NODE_CALL, str(args.wasm_module.resolve()), policy]),
            ]
        commands += [('openbabel', 'default', [str(args.obabel), '-icml', '-ocml']),
                     ('openbabel', 'properties_metadata', [str(args.obabel), '-icml', '-ocml', '-xp', '-xm'])]
        for engine, policy, command in commands:
            run = bounded_run(command, source.encode(), timeout=20, max_output=262144)
            row = {'case': case['id'], 'engine': engine, 'policy': policy, 'command': command,
                   **{k: v.decode('utf8', errors='replace') if isinstance(v, bytes) else v for k, v in run.items()},
                   'terminal': 'process_failure'}
            if run['resource_status'] != 'completed':
                row['terminal'] = run['resource_status']
            elif run['exit_code'] == 0:
                try:
                    if engine != 'openbabel':
                        report = json.loads(row['stdout'])
                        output = check_report(report, source, policy)
                        reports[engine, policy] = report
                    else:
                        output = row['stdout']
                    row['output'] = output
                    row['exact_source_bytes'] = output == source
                    row.update(compare(expected, snapshot(output)))
                    row['terminal'] = 'selected_fields_retained' if row['all_selected_fields_retained'] else 'selected_field_loss'
                except (ValueError, KeyError, TypeError, IndexError, ET.ParseError) as error:
                    row['terminal'], row['error'] = 'unresolved_output', str(error)
            rows.append(row)
        for policy in ['strict', 'allow_loss', 'preserve_source']:
            same = all((e, policy) in reports for e in ['cli', 'python', 'wasm'])
            same = same and reports['cli', policy] == reports['python', policy] == reports['wasm', policy]
            parity.append({'case': case['id'], 'policy': policy, 'identical_report': same})
    groups = [(e, p) for e in ['cli', 'python', 'wasm'] for p in ['strict', 'allow_loss', 'preserve_source']]
    groups += [('openbabel', 'default'), ('openbabel', 'properties_metadata')]
    result = {'schema': 'chematic.cml-source-conversion-development.v1',
              'scope': 'Original six-document development packet: explicit C/O single-bond graph, supplied 2D/3D component presence and numbers, source atom attributes and expanded property/extension payload. Complete source is retained opaquely; no verified full chemistry/stereo, edits, cross-format or independent holdout claim. Failures are not coverage.',
              'observer_revision': 2,
              'manifest': identity(manifest), 'fixture_metadata': {k: v for k, v in packet.items() if k != 'cases'},
              'oracle': {'name': 'Python ElementTree plus authored field expectations', 'version': platform.python_version(), 'script': identity(__file__),
                         'absolute_coordinate_tolerance_angstrom': TOLERANCE, 'relative_tolerance': 0, 'alignment': 'none',
                         'not_scored': 'Molecule root ID/title, inferred H, added payload outside original source subtrees, schema/chemical validity, stereo interpretation; comments/PIs/prefix and number spelling are exact-source diagnostics only.'},
              'tools': {'cli': identity(args.cli), 'wasm_js': identity(args.wasm_module), 'openbabel': identity(args.obabel)},
              'openbabel_version': version['stdout'].decode().strip(), 'cases': packet['cases'], 'rows': rows, 'binding_parity': parity,
              'counts': {e+'/'+p: dict(Counter(r['terminal'] for r in rows if (r['engine'], r['policy']) == (e, p))) for e, p in groups}}
    with args.output.open('x') as output:
        output.write(json.dumps(result, ensure_ascii=False, indent=2, allow_nan=False)+'\n')
    print(json.dumps({'counts': result['counts'], 'identical_reports': sum(r['identical_report'] for r in parity)}, sort_keys=True))
    require(all(p['identical_report'] for p in parity), 'binding reports differ or a binding failed')


if __name__ == '__main__':
    main()
