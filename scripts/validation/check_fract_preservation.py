#!/usr/bin/env python3
"""Original fract packet, candidate bindings/CLI and installed Open Babel rewrite.

This measures fresh-process same-format workflows, not isolated parsers, speed,
all crystallography or a released package. Comparator text is observed by a
small independent field splitter, never the candidate reader.
"""
from __future__ import annotations
import argparse
import hashlib
import importlib
import json
from pathlib import Path
import subprocess
import sys
import tempfile

ROOT = Path(__file__).resolve().parents[2]
PACKET = ROOT / 'validation/fract_preservation/cases.json'
NODE = r'''
const fs=require('node:fs'),w=require(process.argv[1]);
const p=JSON.parse(fs.readFileSync(0,'utf8'));
const out={valid:[],refused:[]};
for(const c of p.valid){
 try{
  let r=JSON.parse(w.parse_fract_json(c.text));const states=[r];
  for(let i=0;i<3;i++){r=JSON.parse(w.parse_fract_json(w.write_fract_from_json(JSON.stringify(r))));states.push(r);}
  out.valid.push({id:c.id,states});
 }catch(e){out.valid.push({id:c.id,error:String(e)});}
}
for(const c of p.refused){
 try{out.refused.push({id:c.id,status:'accepted',output:w.parse_fract_json(c.text)});}
 catch(e){out.refused.push({id:c.id,status:!(e instanceof WebAssembly.RuntimeError)&&String(e).startsWith('fract line ')?'refused':'error',error:String(e)});}
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
    return result['returncode'] == 2 and not result['stdout'] and result['stderr'].startswith('chematic: fract line ')


def observe(text):
    """Whitespace-only observation, deliberately no candidate parser/validation."""
    lines = text.splitlines()
    if len(lines) < 3:
        raise ValueError('missing output record')
    cell = lines[1].replace(',', ' ').split()
    rows = [row.split() for row in lines[2:] if row.strip()]
    if len(cell) != 6 or not rows or any(len(row) != 4 for row in rows):
        raise ValueError('unexpected output shape')
    return {'title': lines[0], 'cell': cell,
            'atoms': [{'element': row[0], 'fractional': row[1:]} for row in rows]}


def numbers(record):
    return [*record['cell'], *[v for atom in record['atoms'] for v in atom['fractional']]]


def nominal_bits(record):
    return [float(token.split('(')[0]).hex() for token in numbers(record)]


def comparisons(actual, expected):
    source = numbers(expected)
    target = numbers(actual)
    uncertain = [i for i, token in enumerate(source) if '(' in token]
    result = {
        'full_token_record_exact': actual == expected,
        'nominal_numbers_exact': nominal_bits(actual) == nominal_bits(expected),
        'title_elements_order_exact': actual['title'] == expected['title'] and
            [a['element'] for a in actual['atoms']] == [a['element'] for a in expected['atoms']],
        'uncertainty_tokens_exact': all(i < len(target) and source[i] == target[i] for i in uncertain) if uncertain else None,
    }
    result['declared_state_exact'] = result['nominal_numbers_exact'] and result['title_elements_order_exact'] and result['uncertainty_tokens_exact'] is not False
    return result


def python_lane(packet, module):
    valid, refused = [], []
    for case in packet['valid']:
        record = module.parse_fract(case['text'])
        states = [record]
        for _ in range(3):
            record = module.parse_fract(module.write_fract(record))
            states.append(record)
        valid.append({'id': case['id'], 'states': states})
    for case in packet['refused']:
        try:
            result = module.parse_fract(case['text'])
            refused.append({'id': case['id'], 'status': 'accepted', 'output': result})
        except Exception as exc:
            status = 'refused' if isinstance(exc, ValueError) and str(exc).startswith('fract line ') else 'error'
            refused.append({'id': case['id'], 'status': status, 'error': repr(exc)})
    return {'valid': valid, 'refused': refused}


def cli_lane(packet, cli):
    valid, refused = [], []
    for case in packet['valid']:
        receipt = run([cli, 'fract', '--inspect'], case['text'])
        if receipt['returncode']:
            raise RuntimeError(receipt)
        record = json.loads(receipt['stdout'])
        states = [record]
        for _ in range(3):
            written = run([cli, 'fract', '--from-json'], json.dumps(record))
            parsed = run([cli, 'fract', '--inspect'], written['stdout'])
            if written['returncode'] or parsed['returncode']:
                raise RuntimeError([written, parsed])
            record = json.loads(parsed['stdout'])
            states.append(record)
        valid.append({'id': case['id'], 'states': states})
    for case in packet['refused']:
        receipt = run([cli, 'fract'], case['text'])
        refused.append({'id': case['id'], 'status': 'refused' if cli_refusal(receipt) else 'error', 'receipt': receipt})
    return {'valid': valid, 'refused': refused}


def openbabel_lane(packet, binary):
    rows = []
    with tempfile.TemporaryDirectory(prefix='chematic-fract-') as directory:
        for case in packet['valid']:
            source = Path(directory) / (case['id'] + '.fract')
            output = Path(directory) / (case['id'] + '.out.fract')
            source.write_text(case['text'], encoding='utf8')
            receipt = run([binary, '-ifract', source, '-ofract', '-O', output, '-ab'])
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


def unit_probe(binary):
    with tempfile.TemporaryDirectory(prefix='chematic-fract-unit-') as directory:
        source = Path(directory) / 'units.fract'
        source.write_text('units\n4 5 6 90 90 90\nC .25 .5 .75\n')
        receipt = run([binary, '-ifract', source, '-oxyz', '-ab'])
        try:
            xyz = [float(v) for v in receipt['stdout'].splitlines()[2].split()[1:]]
        except (ValueError, IndexError):
            xyz = None
        return {'source': source.read_text(), 'receipt': receipt,
                'expected_cartesian_angstrom': [1.0, 2.5, 4.5],
                'observed': xyz, 'exact': receipt['returncode'] == 0 and xyz == [1.0, 2.5, 4.5]}


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
    emitted = {'valid': [dict(case, text=module.write_fract(case['expected'])) for case in packet['valid']]}
    external_emitted = openbabel_lane(emitted, args.openbabel)
    keys = ['declared_state_exact', 'full_token_record_exact', 'nominal_numbers_exact', 'title_elements_order_exact', 'uncertainty_tokens_exact']
    observed = [r for r in external if r['status'] == 'observed']
    ob_summary = {'total': len(external), 'observed': len(observed)}
    for key in keys:
        values = [r['comparisons'][key] for r in observed if r['comparisons'][key] is not None]
        ob_summary[key] = {'exact': sum(values), 'total': len(values)}
    report = {
        'schema_version': 1, 'scope': 'Original single-cell fract same-format rewrite; local source candidate. No timing, parser isolation or toolkit-wide claim.',
        'source_sha256': {str(path.relative_to(ROOT)): sha(path) for path in [
            ROOT / 'crates/chematic-mol/src/fract.rs', ROOT / 'crates/chematic-mol/src/fract/reader.rs',
            ROOT / 'crates/chematic-py/src/formats/fract.rs', ROOT / 'crates/chematic-wasm/src/fract.rs',
            ROOT / 'crates/chematic-cli/src/fract.rs', Path(__file__).resolve()]},
        'packet_sha256': sha(PACKET), 'candidate_base': '7704ee2030dda8f9ca54a07eb3fc22e94623fb49',
        'runtime_sha256': {'python_extension': sha(module.chematic.__file__), 'wasm': sha(args.wasm_module.with_name('chematic_wasm_bg.wasm')), 'cli': sha(args.cli)},
        'versions': {'python': sys.version, 'node': run(['node', '--version']), 'openbabel': run([args.openbabel, '-V'])},
        'summary': {'candidate': summary, 'openbabel': ob_summary},
        'candidate_receipts': lanes, 'openbabel_receipts': external,
        'openbabel_fract_inventory': run([args.openbabel, '-L', 'fract']),
        'unit_probe': unit_probe(args.openbabel),
        'openbabel_reads_candidate_output': {
            'scope': 'Acceptance through OB fract read/write only; not isolated reader accuracy.',
            'observed': sum(row['status'] == 'observed' for row in external_emitted),
            'total': len(external_emitted), 'receipts': external_emitted},
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, ensure_ascii=False, indent=2) + '\n')
    print(json.dumps(report['summary'], ensure_ascii=False, indent=2))
    good = all(r['exact_all_four_states'] == r['valid_total'] and r['typed_refused'] == r['refusal_total'] for r in summary.values())
    return 0 if good and report['unit_probe']['exact'] and all(row['status'] == 'observed' for row in external_emitted) else 1


if __name__ == '__main__':
    raise SystemExit(main())
