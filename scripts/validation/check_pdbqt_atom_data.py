#!/usr/bin/env python3
"""Record numeric PDBQT API/bridge outcomes with an optional Open Babel reader.

Run with a built package on PYTHONPATH. --baseline records old API absence and
conversion failure without treating those expected observations as gate failures.
The Open Babel probe is compiled from the existing original diagnostic; it reads
stored charges with automatic bonding disabled, not default ligand preparation.
"""
from __future__ import annotations

import argparse
import importlib
import json
from pathlib import Path
import sys
import struct

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
from scripts.checks.check_interchange_preservation import bounded_run, identity, require


def read_probe(stdout):
    lines = stdout.splitlines()
    require(bool(lines), 'probe returned no atom count')
    count = int(lines[0])
    require(len(lines) == count + 1, 'probe count mismatch')
    rows = [line.split() for line in lines[1:]]
    require(all(len(row) == 5 for row in rows), 'probe row shape mismatch')
    return dict(atomic_numbers=[int(row[0]) for row in rows],
                coords=[[float(x) for x in row[1:4]] for row in rows],
                partial_charges=[float(row[4]) for row in rows])


def numeric_data_equal(observed, expected):
    """Compare indexed data, including signed-zero and binary64 charge precision."""
    if observed != expected:
        return False
    def bits(data):
        values = [v for point in data['coords'] for v in point] + data['partial_charges']
        return [struct.pack('>d', v) for v in values]
    return bits(observed) == bits(expected)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', required=True, type=Path)
    parser.add_argument('--openbabel-probe', type=Path)
    parser.add_argument('--baseline', action='store_true')
    args = parser.parse_args()
    import chematic
    extension = importlib.import_module('chematic.chematic')
    packet_path = ROOT / 'validation/pdbqt_atom_data/cases.json'
    packet = json.loads(packet_path.read_text())
    result = dict(schema='chematic.pdbqt-atom-data-validation.v1',
                  baseline=args.baseline, python=sys.version,
                  packet=identity(packet_path), package=identity(chematic.__file__),
                  extension=identity(extension.__file__), cases=[])
    if args.openbabel_probe:
        result['openbabel_probe'] = identity(args.openbabel_probe)
    for case in packet['cases']:
        source, expected = case['input'], case['expected']
        row = dict(id=case['id'], profile=case['profile'], expected=expected)
        try:
            mol, coords, charges = chematic.from_pdbqt_with_data(source)
            numbers = json.loads(mol.to_cjson(coords))['atoms']['elements']['number']
            row['python_data'] = dict(atomic_numbers=numbers, coords=coords, partial_charges=charges)
        except AttributeError as error:
            row['python_data'] = dict(error=str(error))
        row['python_data_exact'] = numeric_data_equal(row['python_data'], expected)
        try:
            output = chematic.convert_format(source, 'pdbqt', 'pdbqt')
            atoms = [s for s in output.splitlines() if s.startswith('ATOM')]
            row['bridge'] = dict(output=output,
                coords=[[float(s[a:b]) for a, b in [(30, 38), (38, 46), (46, 54)]] for s in atoms],
                partial_charges=[float(s[70:76]) for s in atoms])
            row['bridge_preserved'] = (row['bridge']['coords'] == expected['coords']
                and len(atoms) == len(expected['partial_charges'])
                and all(abs(a-b) <= 0.000500000001 for a,b in zip(row['bridge']['partial_charges'], expected['partial_charges'])))
        except ValueError as error:
            row['bridge'] = dict(error=str(error))
            row['bridge_preserved'] = False
        # Existing API can be given coordinates, but previously defaulted charges to zero.
        output = chematic.convert_format(source, 'pdbqt', 'pdbqt', coords=expected['coords'])
        row['bridge_with_coords_charges'] = [float(s[70:76]) for s in output.splitlines() if s.startswith('ATOM')]
        if args.openbabel_probe:
            raw = bounded_run([str(args.openbabel_probe.resolve())], source.encode())
            probe = {k: v.decode('utf-8', errors='replace') if isinstance(v, bytes) else v for k,v in raw.items()}
            row['openbabel'] = probe
            if raw['resource_status'] == 'completed' and raw['exit_code'] == 0:
                probe['data'] = read_probe(probe['stdout'])
            row['openbabel_data_exact'] = numeric_data_equal(probe.get('data'), expected)
        result['cases'].append(row)
    result['summary'] = dict(cases=len(result['cases']),
        python_data_exact=sum(r['python_data_exact'] for r in result['cases']),
        bridge_preserved=sum(r['bridge_preserved'] for r in result['cases']))
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps(result['summary']))
    if not args.baseline:
        require(all(r['python_data_exact'] and r['bridge_preserved'] for r in result['cases']), 'candidate numeric contract failed')


if __name__ == '__main__':
    main()
