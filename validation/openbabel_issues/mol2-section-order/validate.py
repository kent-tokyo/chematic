#!/usr/bin/env python3
"""Compare actual native/binding output and cross-read using installed Open Babel.

Run after the native dump example and Node test have written the result files.
Set PYTHONPATH to the freshly built Python package. This checks only the small
original fixture set; it does not run Antechamber or assert general MOL2 parity.
"""
import argparse
from collections import Counter
import hashlib
import json
from pathlib import Path
import subprocess

import chematic

FIXTURES = Path(__file__).resolve().parent


def sections(text):
    chunks = []
    for line in text.splitlines(keepends=True):
        if line.rstrip('\r\n').startswith('@<TRIPOS>'):
            chunks.append(line)
        else:
            assert chunks, 'unexpected section preamble'
            chunks[-1] += line
    return chunks


def external_read(obabel, output):
    proc = subprocess.run([obabel, '-imol2', '-omol'], input=output,
                          text=True, capture_output=True, check=True)
    lines = proc.stdout.splitlines()
    atom_count, bond_count = int(lines[3][:3]), int(lines[3][3:6])
    charge_codes = {0: 0, 1: 3, 2: 2, 3: 1, 4: 0, 5: -1, 6: -2, 7: -3}
    charges = [charge_codes[int(line[36:39])] for line in lines[4:4+atom_count]]
    for line in lines:
        if line.startswith('M  CHG'):
            fields = line.split()
            count = int(fields[2])
            assert len(fields) == 3 + count * 2
            for i in range(count):
                charges[int(fields[3 + 2*i])-1] = int(fields[4 + 2*i])
    return {'atoms': atom_count, 'bonds': bond_count, 'formal_charges': charges,
            'mol': proc.stdout, 'stderr': proc.stderr}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--results', type=Path, required=True)
    parser.add_argument('--obabel', default='obabel')
    args = parser.parse_args()
    fixture = json.loads((FIXTURES / 'cases.json').read_text())
    baseline = json.loads((args.results / 'baseline-native.json').read_text())['rows']
    candidate = json.loads((args.results / 'candidate-native.json').read_text())['rows']
    node = {row['id']: row for row in json.loads(
        (args.results / 'candidate-node.json').read_text())['rows']}
    assert len(baseline) == len(candidate) == len(fixture['cases']) == 9
    assert set(node) == {case['id'] for case in fixture['cases'] if 'smiles' in case}
    rows = []
    python_rows = []
    for case, old, new in zip(fixture['cases'], baseline, candidate, strict=True):
        assert case['id'] == old['id'] == new['id']
        assert Counter(sections(old['output'])) == Counter(sections(new['output']))
        assert new['default_output'] == old['output']
        assert {k: v for k, v in old.items() if k != 'output'} == {
            k: v for k, v in new.items() if k not in ('output', 'default_output')}
        assert new['formal_charges'] == case['formal_charges']
        markers = [chunk.splitlines()[0] for chunk in sections(new['output'])]
        assert markers[:3] == ['@<TRIPOS>MOLECULE', '@<TRIPOS>ATOM', '@<TRIPOS>BOND']
        if not any(case['formal_charges']):
            assert old['output'] == new['output']
        row = {'id': case['id'], 'section_content_identical': True,
               'bytes_changed': old['output'] != new['output']}
        if 'smiles' in case:
            mol = chematic.from_smiles(case['smiles'])
            output = mol.to_mol2(atom_bond_adjacent=True)
            assert mol.to_mol2() == old['output']
            assert output == new['output'] == node[case['id']]['output']
            assert node[case['id']]['default_output'] == old['output']
            python_rows.append({'id': case['id'], 'output': output, 'default_output': mol.to_mol2()})
            row['binding_bytes_identical'] = True
        # Both profiles are measured. Open Babel 3.2.1 drops UNITY after BOND,
        # so this is a disclosed diagnostic, not a general cross-read pass.
        external_default = external_read(args.obabel, new['default_output'])
        external_adjacent = external_read(args.obabel, new['output'])
        for external in [external_default, external_adjacent]:
            assert (external['atoms'], external['bonds']) == (new['atoms'], new['bonds'])
        assert external_default['formal_charges'] == case['formal_charges']
        row.update({'openbabel_default': external_default,
                    'openbabel_adjacent': external_adjacent,
                    'openbabel_adjacent_charge_match':
                        external_adjacent['formal_charges'] == case['formal_charges']})
        rows.append(row)
    source = (FIXTURES / 'charged-source.mol2').read_text()
    for text in [source, source.replace('\n', '\r\n')]:
        assert chematic.roundtrip_mol2_record(text) == text
        assert chematic.convert_format(text, 'mol2', 'mol2') == text
    version = subprocess.check_output([args.obabel, '-V'], text=True).strip()
    probe = subprocess.run([args.obabel, '-:[O-2]', '-omol2'], text=True,
                           capture_output=True, check=True)
    assert probe.stdout.index('@<TRIPOS>UNITY_ATOM_ATTR') < probe.stdout.index('@<TRIPOS>BOND')
    report = {
        'schema': 'chematic.mol2-section-order-validation.v1',
        'baseline_commit': '3e7f04b20593af4b07569d886b56afc1184468fe',
        'scope': '9 original native controls; 8 Python/Node byte parity and 9 two-profile Open Babel diagnostics',
        'limits': ['Antechamber not installed or executed',
                   'No full Amber workflow, broad chemistry parity, speed, or publication claim'],
        'native_section_content_checks': len(rows),
        'opt_in_changed_outputs': sum(row['bytes_changed'] for row in rows),
        'default_outputs_byte_identical': 9,
        'openbabel_default_charge_matches': 9,
        'openbabel_adjacent_charge_mismatches': sum(
            not row['openbabel_adjacent_charge_match'] for row in rows),
        'neutral_outputs_byte_identical': 2,
        'python_node_native_parity_cases': len(node),
        'source_exact_checks': 4,
        'openbabel_version': version,
        'openbabel_default_probe': {'output': probe.stdout, 'stderr': probe.stderr},
        'fixtures_sha256': {p.name: hashlib.sha256(p.read_bytes()).hexdigest()
                            for p in sorted(FIXTURES.glob('*')) if p.is_file()},
        'rows': rows,
    }
    (args.results / 'candidate-python.json').write_text(
        json.dumps({'rows': python_rows}, indent=2) + '\n')
    (args.results / 'cross-read-validation.json').write_text(json.dumps(report, indent=2) + '\n')
    print(f"{len(rows)} section comparisons, {len(node)} binding byte comparisons passed; Open Babel charge differences recorded")


if __name__ == '__main__':
    main()
