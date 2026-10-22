#!/usr/bin/env python3
"""Paired CML 2D coordinate rewrites on an original development packet.

Run as python3 -m scripts.validation.check_cml_coordinate_precision --probe
PROBE --obabel OBABEL --output NEW.json. No independent holdout or timing claim.
"""
from __future__ import annotations

import argparse
from collections import Counter
import hashlib
import json
import math
from pathlib import Path
import platform
import tempfile
import xml.etree.ElementTree as ET

from scripts.check_interchange_preservation import bounded_run, identity, require

ABS_TOLERANCE = 1e-10  # angstrom, no relative tolerance or alignment
COORDINATES = [
    ("dyadic-control", [(1.25, -2.5), (2.75, -2.5)]),
    ("decimal-control", [(0.1, -0.2), (1.5, -0.2)]),
    ("fractional-layout", [(0.123456789012345, 0.234567890123456),
                           (1.523456789012345, 0.334567890123456)]),
    ("opposite-signs", [(-1.23456789012345, 2.34567890123456),
                        (0.16543210987655, 2.44567890123456)]),
    ("near-origin", [(1.23456789012345e-8, -2.34567890123456e-8),
                     (1.400000012345679, 0.099999976543211)]),
    ("translated-layout", [(1000000.123456789, -1000000.234567891),
                           (1000001.523456789, -1000000.134567891)]),
    ("near-rounding-boundary", [(0.000049999999999, -0.000050000000001),
                                (1.400049999999999, 0.099949999999999)]),
    ("long-decimals", [(0.3333333333333333, -0.6666666666666666),
                       (1.7333333333333334, -0.5666666666666667)]),
]


def packet():
    cases = []
    for name, coords in COORDINATES:
        atoms = "\n".join(
            f'<atom id="a{i}" elementType="{element}" x2="{x!r}" y2="{y!r}"/>'
            for i, (element, (x, y)) in enumerate(zip(["C", "O"], coords), 1)
        )
        source = ('<molecule xmlns="http://www.xml-cml.org/schema">\n<atomArray>\n'
                  + atoms + '\n</atomArray>\n<bondArray>\n'
                  '<bond atomRefs2="a1 a2" order="1"/>\n</bondArray>\n</molecule>\n')
        cases.append({"id": name, "source": source, "expected_coordinates": coords,
                      "sha256": hashlib.sha256(source.encode()).hexdigest()})
    return cases


def coordinates(text):
    """Independent XML observer: map the unique C/O atoms, then check topology.

    Missing/3D-only coordinates do not count as supplied 2D zero coordinates.
    Atom ID renaming and atom/bond serialization order are immaterial.
    """
    root = ET.fromstring(text)
    atoms, ids, bonds = {}, {}, []
    for node in root.iter():
        if node.tag.split('}')[-1] == 'atom':
            a = node.attrib
            element = a['elementType']
            require(element in ('C', 'O') and element not in atoms, 'unexpected atom identity')
            require(a['id'] not in ids, 'duplicate atom ID')
            require(int(a.get('formalCharge', 0)) == 0, 'formal charge changed')
            require(int(a.get('isotopeNumber', a.get('isotope', 0))) == 0, 'isotope changed')
            xy = [float(a['x2']), float(a['y2'])]
            require(all(math.isfinite(v) for v in xy), 'nonfinite coordinate')
            atoms[element], ids[a['id']] = xy, element
    for node in root.iter():
        if node.tag.split('}')[-1] == 'bond':
            a, b = node.attrib['atomRefs2'].split()
            require(node.attrib['order'] in ('1', 'S'), 'bond order changed')
            bonds.append(sorted((ids[a], ids[b])))
    require(set(atoms) == {'C', 'O'} and bonds == [['C', 'O']], 'graph changed')
    return [atoms['C'], atoms['O']]


def compare(expected, actual):
    require(len(expected) == len(actual) == 2, 'coordinate count changed')
    require(all(len(row) == 2 for row in [*expected, *actual]), 'coordinate dimension changed')
    deltas = [abs(want - got) for a, b in zip(expected, actual) for want, got in zip(a, b)]
    require(all(math.isfinite(v) for v in deltas), 'nonfinite delta')
    return {"max_absolute_error_angstrom": max(deltas),
            "exact_binary64_values": all(delta == 0 for delta in deltas),
            "within_tolerance": max(deltas) <= ABS_TOLERANCE}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--probe', type=Path, required=True)
    parser.add_argument('--obabel', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    require(not args.output.exists(), 'output already exists')
    version = bounded_run([str(args.obabel), '-V'])
    require(version['exit_code'] == 0 and b'Open Babel 3.2.1' in version['stdout'], 'expected Open Babel 3.2.1')
    options = bounded_run([str(args.obabel), '-L', 'cml'])
    require(options['resource_status'] == 'completed' and options['exit_code'] == 0, 'CML inventory failed')
    cases, rows = packet(), []
    with tempfile.TemporaryDirectory(prefix='chematic-cml-precision-') as directory:
        path, output = Path(directory)/'input.cml', Path(directory)/'output.cml'
        for case in cases:
            require(len(case['source'].encode()) < 16384, 'packet input too large')
            require(coordinates(case['source']) == [list(xy) for xy in case['expected_coordinates']], 'bad source oracle')
            path.write_text(case['source'])
            for engine, command in [
                ('chematic', [str(args.probe), '--format', 'cml', '--path', str(path), '--rewrite-output', str(output)]),
                ('openbabel', [str(args.obabel), '-icml', '-ocml']),
            ]:
                output.unlink(missing_ok=True)
                run = bounded_run(command, case['source'].encode(), timeout=10, max_output=262144)
                row = {'case': case['id'], 'engine': engine, 'command': command,
                       **{k: v.decode('utf8', errors='replace') if isinstance(v, bytes) else v for k, v in run.items()}}
                row['terminal'] = 'process_failure'
                if run['resource_status'] != 'completed':
                    row['terminal'] = run['resource_status']
                elif run['exit_code'] == 0:
                    try:
                        if engine == 'chematic':
                            require(output.stat().st_size <= 262144, 'output too large')
                        row['output'] = output.read_text() if engine == 'chematic' else row['stdout']
                        row['coordinates'] = coordinates(row['output'])
                        row.update(compare(case['expected_coordinates'], row['coordinates']))
                        row['terminal'] = 'preserved' if row['within_tolerance'] else 'coordinate_loss'
                    except (ValueError, KeyError, IndexError, OSError, ET.ParseError) as error:
                        row['terminal'], row['error'] = 'unresolved_output', str(error)
                rows.append(row)
    result = {
        'schema': 'chematic.cml-coordinate-precision-development.v1',
        'scope': 'CML→CML supplied 2D coordinates and the two-atom methanol graph only. Default Open Babel CLI options. No stereo/metadata/3D, independent holdout, timing or toolkit-wide claim; failures/refusals are not preservation.',
        'provenance': 'Original synthetic layouts, 2026-10-22; packet() recipe, no third-party or sealed input.',
        'license': 'MIT OR Apache-2.0', 'redistribution': 'allowed: original inputs',
        'minimization': 'Two distinct-element atoms provide unambiguous correspondence and one bond; layouts isolate decimal rounding, signs, translation and proximity to the origin. Not optimized geometry.',
        'oracle': {'name': 'Python ElementTree plus authored binary64 coordinates',
                   'version': platform.python_version(), 'script': identity(__file__),
                   'absolute_tolerance_angstrom': ABS_TOLERANCE, 'relative_tolerance': 0,
                   'alignment': 'none', 'exact_value_metric': 'separate diagnostic'},
        'engines': {'chematic': identity(args.probe), 'openbabel': identity(args.obabel)},
        'openbabel_version': version['stdout'].decode().strip(),
        'openbabel_cml_options': {k: v.decode('utf8', errors='replace') if isinstance(v, bytes) else v for k, v in options.items()},
        'cases': cases, 'rows': rows,
        'counts': {e: dict(Counter(r['terminal'] for r in rows if r['engine'] == e)) for e in ['chematic', 'openbabel']},
    }
    with args.output.open('x') as output:
        output.write(json.dumps(result, indent=2, allow_nan=False)+'\n')
    print(json.dumps(result['counts'], sort_keys=True))


if __name__ == '__main__':
    main()
