#!/usr/bin/env python3
"""Independent Decimal geometry checks and pinned Open Babel #2822 diagnostic."""
import argparse
from decimal import Decimal, localcontext
import hashlib
import itertools
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import chematic

ROOT = Path(__file__).resolve().parents[2]
FIXTURE = ROOT / 'validation/openbabel_issues/turbomole/lattice.json'
BOHR = Decimal('0.529177210544')
TOL = Decimal('1e-12')


def digest(data):
    return hashlib.sha256(data).hexdigest()


def near(actual, expected):
    assert abs(Decimal(str(actual)) - expected) <= TOL, (actual, expected)


def check_written(text, matrix, positions, elements, flags):
    """Decode only the fixed canonical writer layout, independently of Rust."""
    rows = text.splitlines()
    assert rows[0:2] == ['$periodic 3', '$lattice']
    for row, expected in zip(rows[2:5], matrix):
        values = row.split()
        assert len(values) == 3
        for value, wanted in zip(values, expected):
            near(Decimal(value) * BOHR, wanted)
    assert rows[5] == '$coord'
    assert len(rows) == 7 + len(elements) and rows[-1] == '$end'
    for row, expected, element, frozen in zip(rows[6:-1], positions, elements, flags):
        values = row.split()
        assert values[3] == element.lower()
        assert values[4:] == (['f'] if frozen else [])
        for value, wanted in zip(values[:3], expected):
            near(Decimal(value) * BOHR, wanted)


def generated_cases():
    elements = ['Si', 'Cl', 'C', 'H']
    fractions = [[Decimal(x) for x in row] for row in [('.25', '-.5', '1.5'), ('0', '0', '0'), ('-.125', '1.25', '.5'), ('2', '.25', '-.5')]]
    checked = []
    for order, group_order, mode, handedness in itertools.product(itertools.permutations(range(4)), itertools.permutations(range(3)), ['bohr', 'angs', 'frac'], [1, -1]):
        raw_matrix = [[Decimal(x) for x in row] for row in [[4, 1, 0], [-1, 5, 2], [1, 0, handedness * 6]]]
        factor = Decimal(1) if mode == 'angs' else BOHR
        matrix = [[v * factor for v in row] for row in raw_matrix]
        expected = [[sum(fractions[i][j] * matrix[j][k] for j in range(3)) for k in range(3)] for i in order]
        raw_positions = [fractions[i] if mode == 'frac' else [sum(fractions[i][j] * raw_matrix[j][k] for j in range(3)) for k in range(3)] for i in order]
        flags = [i % 2 == 0 for i in order]
        names = [elements[i] for i in order]
        coord = '$coord' + (' ' + mode if mode != 'bohr' else '') + '\n'
        coord += ''.join(' '.join(map(str, p)) + ' ' + e.lower() + (' f' if frozen else '') + '\n' for p, e, frozen in zip(raw_positions, names, flags))
        lattice = '$lattice' + (' angs' if mode == 'angs' else '') + '\n'
        lattice += ''.join(' '.join(map(str, row)) + '\n' for row in raw_matrix)
        groups = ['$periodic 3\n', lattice, coord]
        source = ''.join(groups[i] for i in group_order) + '$end\n'
        g = chematic.parse_turbomole(source)
        assert g['source_text'] == source and g['periodicity'] == 3
        assert g['elements'] == names and g['frozen'] == flags
        for got, wanted in zip(g['positions_angstrom'], expected):
            for a, b in zip(got, wanted): near(a, b)
        for got, wanted in zip(g['lattice_angstrom'], matrix):
            for a, b in zip(got, wanted): near(a, b)
        written = chematic.write_turbomole(g['elements'], g['positions_angstrom'], lattice_angstrom=g['lattice_angstrom'], frozen=g['frozen'])
        check_written(written, matrix, expected, names, flags)
        back = chematic.parse_turbomole(written)
        assert back['elements'] == names and back['frozen'] == flags
        for got, wanted in zip(back['positions_angstrom'], expected):
            for a, b in zip(got, wanted): near(a, b)
        for got, wanted in zip(back['lattice_angstrom'], matrix):
            for a, b in zip(got, wanted): near(a, b)
        checked.append({'order':order, 'group_order':group_order, 'mode':mode, 'handedness':handedness, 'input_sha256':digest(source.encode()), 'output_sha256':digest(written.encode())})
    return checked


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--obabel', default='/opt/homebrew/bin/obabel')
    parser.add_argument('--output', type=Path, default=ROOT/'validation/results/openbabel-turbomole-lattice-2026-10-15.json')
    args = parser.parse_args()
    version = subprocess.run([args.obabel, '-V'], capture_output=True, text=True, check=True, timeout=10).stdout.strip()
    if not version.startswith('Open Babel 3.2.1 '): raise RuntimeError(f'expected Open Babel 3.2.1, got {version}')
    manifest = json.loads(FIXTURE.read_text())
    payload = FIXTURE.parent / manifest['payload']
    text = payload.read_text()
    assert digest(payload.read_bytes()) == manifest['sha256']
    g = chematic.parse_turbomole(text)
    assert g['elements'] == manifest['expected']['elements'] and g['frozen'] == manifest['expected']['frozen']
    for actual, wanted in zip(g['positions_angstrom'][0], manifest['expected']['position_angstrom']):
        near(actual, Decimal(str(wanted)))
    for i, wanted in enumerate(manifest['expected']['lattice_diagonal_angstrom']):
        near(g['lattice_angstrom'][i][i], Decimal(str(wanted)))
    process = subprocess.run([args.obabel, '-itmol', '-ovasp'], input=text, capture_output=True, text=True, check=True, timeout=15)
    # POSCAR first line can be empty (unnamed source), so do not strip it.
    rows = process.stdout.splitlines()
    babel_lattice = [[float(v) for v in row.split()] for row in rows[2:5]]
    assert babel_lattice == [[0., 0., 0.]] * 3, 'upstream reader behavior changed; inspect before rebaseline'
    babel_position = list(map(float, rows[8].split()[:3]))
    assert all(abs(a-b) < 2e-9 for a,b in zip(babel_position, g['positions_angstrom'][0]))
    # Writer direction is a separate control: 3.2.1 already emits lattice.
    poscar = 'synthetic control\n1.0\n4 0 0\n0 5 0\n0 0 6\nSi\n1\nDirect\n.25 .5 .75\n'
    with tempfile.TemporaryDirectory(prefix='chematic-tmol-control-') as directory:
        input_path = Path(directory) / 'POSCAR'
        input_path.write_text(poscar)
        control = subprocess.run([args.obabel, '-ivasp', str(input_path), '-otmol'], capture_output=True, text=True, check=True, timeout=15)
    control_g = chematic.parse_turbomole(control.stdout)
    assert control_g['periodicity'] == 3
    for row, expected in zip(control_g['lattice_angstrom'], [[4,0,0],[0,5,0],[0,0,6]]):
        assert all(abs(a-b) < 1e-8 for a,b in zip(row, expected))
    with localcontext() as context:
        context.prec = 40
        checked = generated_cases()
    paths = ['crates/chematic-mol/src/turbomole.rs', 'crates/chematic-mol/src/turbomole/reader.rs', 'crates/chematic-py/src/formats/turbomole.rs', 'crates/chematic-mol/tests/turbomole_contract.rs', 'crates/chematic-py/tests/test_turbomole.py', 'crates/chematic-py/python/chematic/__init__.pyi', 'scripts/validation/validate_turbomole_lattice.py']
    record = {
        'schema_version':1, 'status':'local-source-candidate', 'issue':'https://github.com/openbabel/openbabel/issues/2822',
        'base':'cac480fb source snapshot plus preserved earlier RDKit/Open Babel candidates through stereo-group inversion; patch isolates TURBOMOLE geometry I/O',
        'scope':'standalone nonperiodic and 3D geometry documents; original lattice orientation, atom order, coordinates and frozen flags; no general graph-conversion registration',
        'runtime':{'python':sys.version.split()[0], 'openbabel':version, 'chematic_extension_sha256':digest(Path(chematic.__file__).resolve().read_bytes())},
        'source_sha256':{p:digest((ROOT/p).read_bytes()) for p in paths},
        'fixture':{'path':str(FIXTURE.relative_to(ROOT)), 'manifest_sha256':digest(FIXTURE.read_bytes()), **manifest},
        'results':{'comparisons':len(checked), 'coordinate_mismatches':0, 'lattice_mismatches':0, 'fixed_flag_mismatches':0, 'source_text_mismatches':0, 'atom_permutations':24, 'group_orders':6, 'coordinate_modes':['bohr','angs','frac'], 'handedness':[1,-1], 'absolute_tolerance_angstrom':str(TOL), 'case_sha256':digest(json.dumps(checked,sort_keys=True,separators=(',',':')).encode())},
        'openbabel_reader_diagnostic':{'command':['obabel','-itmol','-ovasp'], 'output_lattice':babel_lattice, 'expected_lattice_angstrom':g['lattice_angstrom'], 'positions_angstrom':babel_position, 'output_sha256':digest(process.stdout.encode()), 'stderr':process.stderr.strip(), 'coordinate_tolerance_angstrom':2e-9},
        'openbabel_writer_control':{'command':['obabel','-ivasp','<absolute-POSCAR-path>','-otmol'], 'periodicity':control_g['periodicity'], 'retains_lattice':True, 'input_sha256':digest(poscar.encode()), 'output_sha256':digest(control.stdout.encode()), 'numeric_tolerance_angstrom':1e-8, 'note':'Writer direction works in this pinned version; only reader-direction cell loss was reproduced.'},
        'limitations':'TURBOMOLE executable was not run; independent arithmetic and documented layout are the expected-geometry oracle. 1D/2D, control files, external file options, unknown groups, charge/spin/isotope/ghost labels are rejected. $coord angs/bohr are explicitly supported MCTC-style extensions. Comments/formatting survive only in source_text; canonical writer emits fresh bohr data. No topology, bond perception, wrapping, cell rotation or POSCAR adapter added.',
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(record,indent=2)+'\n')
    print(json.dumps(record['results']))


if __name__ == '__main__':
    main()
