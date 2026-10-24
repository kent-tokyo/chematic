#!/usr/bin/env python3
"""Original unixyz packet, candidate bindings/CLI and installed Open Babel rewrite.

This measures fresh-process same-format workflows, not isolated parsers, speed,
all UniChem dialects or a released package. Comparator text is observed by a
small independent field splitter, never the candidate reader.
"""
from __future__ import annotations
import argparse
import hashlib
import importlib
import json
import math
from pathlib import Path
import subprocess
import sys
import tempfile

ROOT = Path(__file__).resolve().parents[2]
PACKET = ROOT / 'validation/unixyz_preservation/cases.json'
NODE = r'''
const fs=require('node:fs'),w=require(process.argv[1]);
const p=JSON.parse(fs.readFileSync(0,'utf8'));
const out={valid:[],refused:[]};
for(const c of p.valid){
 try{
  let r=JSON.parse(w.parse_unixyz_json(c.text));const states=[r];
  for(let i=0;i<3;i++){r=JSON.parse(w.parse_unixyz_json(w.write_unixyz_from_json(JSON.stringify(r))));states.push(r);}
  out.valid.push({id:c.id,states});
 }catch(e){out.valid.push({id:c.id,error:String(e)});}
}
for(const c of p.refused){
 try{out.refused.push({id:c.id,status:'accepted',output:w.parse_unixyz_json(c.text)});}
 catch(e){out.refused.push({id:c.id,status:!(e instanceof WebAssembly.RuntimeError)&&String(e).startsWith('UniChem XYZ line ')?'refused':'error',error:String(e)});}
}
process.stdout.write(JSON.stringify(out));
'''


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def run(args, text=None):
    result = subprocess.run(args, input=text, text=True, capture_output=True, timeout=30)
    return {'argv': [str(a) for a in args], 'returncode': result.returncode,
            'stdout': result.stdout, 'stderr': result.stderr}


def cli_refusal(result):
    return result['returncode'] == 2 and not result['stdout'] and result['stderr'].startswith('chematic: UniChem XYZ line ')


def observe(text):
    """Independent external-output field observation, not the candidate reader."""
    lines = text.splitlines()
    if len(lines) < 3:
        raise ValueError('missing external record')
    n = int(lines[1].strip())
    rows = [line.split() for line in lines[2:]]
    if n <= 0 or len(rows) != n or any(len(row) != 4 for row in rows):
        raise ValueError('incomplete/unexpected external atom section')
    return {'title': lines[0], 'atoms': [{'atomic_number': int(row[0]), 'position_angstrom': row[1:]} for row in rows]}


def nominal_bits(record):
    return [[float(v).hex() for v in atom['position_angstrom']] for atom in record['atoms']]


def comparisons(actual, expected):
    result = {'full_token_record_exact': actual == expected,
              'coordinate_bits_exact': nominal_bits(actual) == nominal_bits(expected),
              'title_exact': actual['title'] == expected['title'],
              'atomic_numbers_order_exact': [a['atomic_number'] for a in actual['atoms']] == [a['atomic_number'] for a in expected['atoms']]}
    result['declared_state_exact'] = all(result[key] for key in ['coordinate_bits_exact','title_exact','atomic_numbers_order_exact'])
    return result


def python_lane(packet, module):
    valid, refused = [], []
    for case in packet['valid']:
        record = module.parse_unixyz(case['text'])
        states = [record]
        for _ in range(3):
            record = module.parse_unixyz(module.write_unixyz(record))
            states.append(record)
        valid.append({'id': case['id'], 'states': states})
    for case in packet['refused']:
        try:
            result = module.parse_unixyz(case['text'])
            refused.append({'id': case['id'], 'status': 'accepted', 'output': result})
        except Exception as exc:
            status = 'refused' if isinstance(exc, ValueError) and str(exc).startswith('UniChem XYZ line ') else 'error'
            refused.append({'id': case['id'], 'status': status, 'error': repr(exc)})
    return {'valid': valid, 'refused': refused}


def cli_lane(packet, cli):
    valid, refused = [], []
    for case in packet['valid']:
        receipt = run([cli, 'unixyz', '--inspect'], case['text'])
        if receipt['returncode']:
            raise RuntimeError(receipt)
        record = json.loads(receipt['stdout'])
        states = [record]
        for _ in range(3):
            written = run([cli, 'unixyz', '--from-json'], json.dumps(record))
            parsed = run([cli, 'unixyz', '--inspect'], written['stdout'])
            if written['returncode'] or parsed['returncode']:
                raise RuntimeError([written, parsed])
            record = json.loads(parsed['stdout'])
            states.append(record)
        valid.append({'id': case['id'], 'states': states})
    for case in packet['refused']:
        receipt = run([cli, 'unixyz'], case['text'])
        refused.append({'id': case['id'], 'status': 'refused' if cli_refusal(receipt) else 'error', 'receipt': receipt})
    return {'valid': valid, 'refused': refused}


def openbabel_lane(packet, binary):
    rows = []
    with tempfile.TemporaryDirectory(prefix='chematic-unixyz-') as directory:
        for case in packet['valid']:
            source = Path(directory) / (case['id'] + '.unixyz')
            output = Path(directory) / (case['id'] + '.out.unixyz')
            source.write_text(case['text'], encoding='utf8')
            receipt = run([binary, '-iunixyz', source, '-ounixyz', '-O', output, '-ab'])
            row = {'id': case['id'], 'receipt': receipt}
            row['output'] = output.read_text() if output.exists() else ''
            try:
                if receipt['returncode'] != 0:
                    raise ValueError('nonzero Open Babel exit')
                actual = observe(row['output'])
                row.update(status='observed', state=actual, comparisons=comparisons(actual, case['expected']))
            except (ValueError, IndexError, KeyError) as exc:
                row.update(status='unobservable', error=str(exc))
            rows.append(row)
    return rows


def candidate_read_probe(packet, module, binary):
    """CJSON output avoids conflating unixyz reader acceptance with its writer's fixed columns."""
    rows = []
    with tempfile.TemporaryDirectory(prefix='chematic-unixyz-read-') as directory:
        for case in packet['valid']:
            source = Path(directory) / (case['id'] + '.unixyz')
            output = Path(directory) / (case['id'] + '.cjson')
            source.write_text(module.write_unixyz(case['expected']))
            receipt = run([binary, '-iunixyz', source, '-ocjson', '-O', output, '-ab'])
            text = output.read_text() if output.exists() else ''
            row = {'id': case['id'], 'receipt': receipt, 'output': text}
            try:
                data = json.loads(text)
                zs = data['atoms']['elements']['number']
                xyz = data['atoms']['coords']['3d']
                expected_zs = [atom['atomic_number'] for atom in case['expected']['atoms']]
                good = receipt['returncode'] == 0 and zs == expected_zs and len(xyz) == 3 * len(zs) and all(math.isfinite(v) for v in xyz)
                row.update(status='observed' if good else 'unexpected', atomic_numbers=zs, coordinates=xyz)
            except (ValueError, KeyError, TypeError) as exc:
                row.update(status='unobservable', error=str(exc))
            rows.append(row)
    return {'scope': 'Regenerated text read and CJSON output. Atom count/order and finite coordinate shape only; not isolated parser accuracy or coordinate-precision proof.',
            'observed': sum(r['status'] == 'observed' for r in rows), 'total': len(rows), 'receipts': rows}


def unit_probe(binary):
    source_text = 'units\n1\n6 .125 .25 .5\n'
    xyz_text = '1\nlayout title\nC .125 .25 .5\n'
    with tempfile.TemporaryDirectory(prefix='chematic-unixyz-unit-') as directory:
        source = Path(directory) / 'units.unixyz'
        source.write_text(source_text)
        receipt = run([binary, '-iunixyz', source, '-oxyz', '-ab'])
        layout_receipt = run([binary, '-ixyz', '-ounixyz', '-ab'], xyz_text)
        try:
            xyz = [float(v) for v in receipt['stdout'].splitlines()[2].split()[1:]]
            layout = observe(layout_receipt['stdout'])
        except (ValueError, IndexError):
            xyz, layout = None, None
        expected = {'title': 'layout title', 'atoms': [{'atomic_number': 6, 'position_angstrom': ['0.12500','0.25000','0.50000']}]}
        return {'source': source_text, 'receipt': receipt, 'expected_cartesian_angstrom': [.125,.25,.5],
                'observed': xyz, 'layout_source_xyz': xyz_text, 'layout_receipt': layout_receipt,
                'expected_layout': expected, 'layout_observed': layout,
                'exact': receipt['returncode'] == 0 and xyz == [.125,.25,.5] and layout_receipt['returncode'] == 0 and layout == expected}


def summarize_lane(lane, packet):
    expected = {c['id']: c['expected'] for c in packet['valid']}
    exact = sum('states' in row and len(row['states']) == 4 and
                all(state == expected[row['id']] for state in row['states']) for row in lane['valid'])
    return {'valid_total': len(packet['valid']), 'exact_all_four_states': exact,
            'refusal_total': len(packet['refused']),
            'typed_refused': sum(row['status'] == 'refused' for row in lane['refused'])}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--python-path', type=Path, required=True)
    parser.add_argument('--wasm-module', type=Path, required=True)
    parser.add_argument('--cli', type=Path, required=True)
    parser.add_argument('--openbabel', default='/opt/homebrew/bin/obabel')
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    packet = json.loads(PACKET.read_text())
    sys.path.insert(0, str(args.python_path.resolve()))
    module = importlib.import_module('chematic')
    node = run(['node', '-e', NODE, str(args.wasm_module.resolve())], json.dumps(packet))
    if node['returncode']:
        raise RuntimeError(node)
    lanes = {'python': python_lane(packet, module), 'wasm_node': json.loads(node['stdout']),
             'cli': cli_lane(packet, str(args.cli.resolve()))}
    summary = {name: summarize_lane(lane, packet) for name, lane in lanes.items()}
    external = openbabel_lane(packet, args.openbabel)
    emitted = {'valid': [dict(case, text=module.write_unixyz(case['expected'])) for case in packet['valid']]}
    external_emitted = openbabel_lane(emitted, args.openbabel)
    keys = ['declared_state_exact', 'full_token_record_exact', 'coordinate_bits_exact', 'title_exact', 'atomic_numbers_order_exact']
    observed = [r for r in external if r['status'] == 'observed']
    ob_summary = {'total': len(external), 'observed': len(observed)}
    for key in keys:
        values = [r['comparisons'][key] for r in observed if r['comparisons'][key] is not None]
        ob_summary[key] = {'exact': sum(values), 'total': len(values)}
    report = {
        'schema_version': 1, 'scope': 'Original single-record UniChem XYZ same-format rewrite; local source candidate. No timing, parser isolation or toolkit-wide claim.',
        'source_sha256': {str(path.relative_to(ROOT)): sha(path) for path in [
            ROOT / 'crates/chematic-mol/src/unixyz.rs', ROOT / 'crates/chematic-mol/src/unixyz/reader.rs',
            ROOT / 'crates/chematic-py/src/formats/unixyz.rs', ROOT / 'crates/chematic-wasm/src/unixyz.rs',
            ROOT / 'crates/chematic-cli/src/unixyz.rs', Path(__file__).resolve()]},
        'packet_sha256': sha(PACKET), 'candidate_base': '7704ee2030dda8f9ca54a07eb3fc22e94623fb49',
        'runtime_sha256': {'python_extension': sha(module.chematic.__file__), 'wasm': sha(args.wasm_module.with_name('chematic_wasm_bg.wasm')), 'cli': sha(args.cli)},
        'versions': {'python': sys.version, 'node': run(['node', '--version']), 'openbabel': run([args.openbabel, '-V'])},
        'summary': {'candidate': summary, 'openbabel': ob_summary},
        'candidate_receipts': lanes, 'openbabel_receipts': external,
        'openbabel_unixyz_inventory': run([args.openbabel, '-L', 'unixyz']),
        'unit_probe': unit_probe(args.openbabel),
        'candidate_read_probe': candidate_read_probe(packet, module, args.openbabel),
        'openbabel_candidate_text_rewrite': {
            'scope': 'Regenerated candidate text rewritten through OB unixyz. Counts describe external output observability, not isolated reader acceptance.',
            'observed': sum(row['status'] == 'observed' for row in external_emitted),
            'total': len(external_emitted), 'receipts': external_emitted},
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, ensure_ascii=False, indent=2) + '\n')
    print(json.dumps(report['summary'], ensure_ascii=False, indent=2))
    good = all(r['exact_all_four_states'] == r['valid_total'] and r['typed_refused'] == r['refusal_total'] for r in summary.values())
    return 0 if good and report['unit_probe']['exact'] and report['candidate_read_probe']['observed'] == len(packet['valid']) else 1


if __name__ == '__main__':
    raise SystemExit(main())
