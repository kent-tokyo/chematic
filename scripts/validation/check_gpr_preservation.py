#!/usr/bin/env python3
"""Original gpr packet, candidate bindings/CLI and installed Open Babel rewrite.

This measures fresh-process same-format workflows, not isolated parsers, speed,
all Ghemical dialects or a released package. Comparator text is observed by a
small independent section observer, never the candidate reader.
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
PACKET = ROOT / 'validation/gpr_preservation/cases.json'
NODE = r'''
const fs=require('node:fs'),w=require(process.argv[1]);
const p=JSON.parse(fs.readFileSync(0,'utf8'));
const out={valid:[],refused:[]};
for(const c of p.valid){
 try{
  let r=JSON.parse(w.parse_gpr_json(c.text));const states=[r];
  for(let i=0;i<3;i++){r=JSON.parse(w.parse_gpr_json(w.write_gpr_from_json(JSON.stringify(r))));states.push(r);}
  out.valid.push({id:c.id,states});
 }catch(e){out.valid.push({id:c.id,error:String(e)});}
}
for(const c of p.refused){
 try{out.refused.push({id:c.id,status:'accepted',output:w.parse_gpr_json(c.text)});}
 catch(e){out.refused.push({id:c.id,status:!(e instanceof WebAssembly.RuntimeError)&&String(e).startsWith('GPR line ')?'refused':'error',error:String(e)});}
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
    return result['returncode'] == 2 and not result['stdout'] and result['stderr'].startswith('chematic: GPR line ')


def observe(text):
    """Independent section/field observation of the external writer output."""
    lines = [row.split() for row in text.splitlines() if row.strip()]
    if len(lines) < 6 or lines[0] != ['!Header', 'gpr', '100'] or lines[1] != ['!Info', '1']:
        raise ValueError('missing/unexpected external header')
    index = 2
    def block(tag, width, counted=True):
        nonlocal index
        header = lines[index]
        if header[0] != tag or len(header) != (2 if counted else 1):
            raise ValueError('unexpected external section: ' + tag)
        count = int(header[1]) if counted else natoms
        index += 1
        rows = lines[index:index+count]
        if len(rows) != count or any(len(row) != width for row in rows):
            raise ValueError('incomplete external section: ' + tag)
        index += count
        return rows
    atoms = block('!Atoms', 2)
    natoms = len(atoms)
    bonds = block('!Bonds', 3)
    coords = block('!Coord', 4, False)
    charges = block('!Charges', 2, False) if lines[index][0] == '!Charges' else None
    if lines[index:] != [['!End']]:
        raise ValueError('unexpected external trailer')
    return {'atoms': [{'id': int(a), 'atomic_number': int(z)} for a,z in atoms],
            'bonds': [{'atoms': [int(a),int(b)], 'order': order} for a,b,order in bonds],
            'coordinates': [{'atom_id': int(row[0]), 'position_nm': row[1:]} for row in coords],
            'partial_charges': None if charges is None else [{'atom_id': int(i), 'value': q} for i,q in charges]}


def coordinate_bits(record):
    return [(row['atom_id'], [float(v).hex() for v in row['position_nm']]) for row in record['coordinates']]


def charge_bits(record):
    charges = record['partial_charges']
    return None if charges is None else [(row['atom_id'], float(row['value']).hex()) for row in charges]


def graph(record):
    return ([(a['id'], a['atomic_number']) for a in record['atoms']],
            sorted((min(b['atoms']), max(b['atoms']), b['order']) for b in record['bonds']))


def comparisons(actual, expected):
    result = {'full_token_record_exact': actual == expected,
              'graph_exact': graph(actual) == graph(expected),
              'coordinate_bits_exact': coordinate_bits(actual) == coordinate_bits(expected),
              'partial_charge_presence_and_bits_exact': charge_bits(actual) == charge_bits(expected)}
    result['declared_state_exact'] = all(result[key] for key in ['graph_exact', 'coordinate_bits_exact', 'partial_charge_presence_and_bits_exact'])
    return result


def python_lane(packet, module):
    valid, refused = [], []
    for case in packet['valid']:
        record = module.parse_gpr(case['text'])
        states = [record]
        for _ in range(3):
            record = module.parse_gpr(module.write_gpr(record))
            states.append(record)
        valid.append({'id': case['id'], 'states': states})
    for case in packet['refused']:
        try:
            result = module.parse_gpr(case['text'])
            refused.append({'id': case['id'], 'status': 'accepted', 'output': result})
        except Exception as exc:
            status = 'refused' if isinstance(exc, ValueError) and str(exc).startswith('GPR line ') else 'error'
            refused.append({'id': case['id'], 'status': status, 'error': repr(exc)})
    return {'valid': valid, 'refused': refused}


def cli_lane(packet, cli):
    valid, refused = [], []
    for case in packet['valid']:
        receipt = run([cli, 'gpr', '--inspect'], case['text'])
        if receipt['returncode']:
            raise RuntimeError(receipt)
        record = json.loads(receipt['stdout'])
        states = [record]
        for _ in range(3):
            written = run([cli, 'gpr', '--from-json'], json.dumps(record))
            parsed = run([cli, 'gpr', '--inspect'], written['stdout'])
            if written['returncode'] or parsed['returncode']:
                raise RuntimeError([written, parsed])
            record = json.loads(parsed['stdout'])
            states.append(record)
        valid.append({'id': case['id'], 'states': states})
    for case in packet['refused']:
        receipt = run([cli, 'gpr'], case['text'])
        refused.append({'id': case['id'], 'status': 'refused' if cli_refusal(receipt) else 'error', 'receipt': receipt})
    return {'valid': valid, 'refused': refused}


def openbabel_lane(packet, binary):
    rows = []
    with tempfile.TemporaryDirectory(prefix='chematic-gpr-') as directory:
        for case in packet['valid']:
            source = Path(directory) / (case['id'] + '.gpr')
            output = Path(directory) / (case['id'] + '.out.gpr')
            source.write_text(case['text'], encoding='utf8')
            receipt = run([binary, '-igpr', source, '-ogpr', '-O', output])
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
    source_text = '!Header gpr 100\n!Info 1\n!Atoms 1\n0 6\n!Bonds 0\n!Coord\n0 .125 .25 .5\n!Charges\n0 -.25\n!End\n'
    with tempfile.TemporaryDirectory(prefix='chematic-gpr-unit-') as directory:
        source = Path(directory) / 'units.gpr'
        source.write_text(source_text)
        receipt = run([binary, '-igpr', source, '-oxyz'])
        try:
            xyz = [float(v) for v in receipt['stdout'].splitlines()[2].split()[1:]]
        except (ValueError, IndexError):
            xyz = None
        return {'source': source_text, 'receipt': receipt,
                'expected_cartesian_angstrom': [1.25, 2.5, 5.0],
                'observed': xyz, 'exact': receipt['returncode'] == 0 and xyz == [1.25, 2.5, 5.0]}


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
    emitted = {'valid': [dict(case, text=module.write_gpr(case['expected'])) for case in packet['valid']]}
    external_emitted = openbabel_lane(emitted, args.openbabel)
    keys = ['declared_state_exact', 'full_token_record_exact', 'graph_exact', 'coordinate_bits_exact', 'partial_charge_presence_and_bits_exact']
    observed = [r for r in external if r['status'] == 'observed']
    ob_summary = {'total': len(external), 'observed': len(observed)}
    for key in keys:
        values = [r['comparisons'][key] for r in observed if r['comparisons'][key] is not None]
        ob_summary[key] = {'exact': sum(values), 'total': len(values)}
    report = {
        'schema_version': 1, 'scope': 'Original gpr 100 / one coordinate set same-format rewrite; local source candidate. No timing, parser isolation or toolkit-wide claim.',
        'source_sha256': {str(path.relative_to(ROOT)): sha(path) for path in [
            ROOT / 'crates/chematic-mol/src/gpr.rs', ROOT / 'crates/chematic-mol/src/gpr/reader.rs', ROOT / 'crates/chematic-mol/src/gpr/writer.rs',
            ROOT / 'crates/chematic-py/src/formats/gpr.rs', ROOT / 'crates/chematic-wasm/src/gpr.rs',
            ROOT / 'crates/chematic-cli/src/gpr.rs', Path(__file__).resolve()]},
        'packet_sha256': sha(PACKET), 'candidate_base': '7704ee2030dda8f9ca54a07eb3fc22e94623fb49',
        'runtime_sha256': {'python_extension': sha(module.chematic.__file__), 'wasm': sha(args.wasm_module.with_name('chematic_wasm_bg.wasm')), 'cli': sha(args.cli)},
        'versions': {'python': sys.version, 'node': run(['node', '--version']), 'openbabel': run([args.openbabel, '-V'])},
        'summary': {'candidate': summary, 'openbabel': ob_summary},
        'candidate_receipts': lanes, 'openbabel_receipts': external,
        'openbabel_gpr_inventory': run([args.openbabel, '-L', 'gpr']),
        'unit_probe': unit_probe(args.openbabel),
        'openbabel_reads_candidate_output': {
            'scope': 'Acceptance through OB gpr read/write only; not isolated reader accuracy.',
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
