#!/usr/bin/env python3
"""Validate original EXYZ geometry against pinned Open Babel, without source reuse."""
import argparse
import hashlib
import itertools
import json
import math
from pathlib import Path
import subprocess
import tempfile

import chematic

ROOT = Path(__file__).resolve().parents[1]


def digest(payload):
    return hashlib.sha256(payload).hexdigest()


def run_obabel(executable, source, output_format):
    command = [executable, '-iexyz', str(source), '-o' + output_format, '-ab']
    run = subprocess.run(command, text=True, capture_output=True, check=True, timeout=10)
    assert '1 molecule converted' in run.stderr, run.stderr
    return run.stdout


def emit(g):
    return chematic.write_exyz(g['elements'], g['positions_angstrom'], header=g['header'],
        lattice_angstrom=g['lattice_angstrom'], offset_angstrom=g['offset_angstrom'],
        virtual_atoms=g['virtual_atoms'])


def assert_close(got, want, tolerance):
    assert len(got) == len(want)
    errors = [abs(a-b) for a,b in zip(got,want)]
    assert all(math.isfinite(v) and v <= tolerance for v in errors), (got, want)
    return max(errors, default=0.0)


def flat(rows):
    return list(itertools.chain.from_iterable(rows))


def fixture(seed, periodic, virtual_style, crlf):
    # Small original geometries; skewed cells include a rotated row basis and
    # both handednesses. No chemical graph or real compound dataset is used.
    symbols = ['C','H','N','O','Si','Cl'][:seed % 6 + 1]
    points = [[math.sin(seed+i) + 0.12345678901234567,
               -0.75*i + seed/7, (-1)**i * (seed+i)/11] for i in range(len(symbols))]
    flags = [virtual_style == 2 or (virtual_style == 1 and i % 2 == 0)
             for i in range(len(symbols))]
    theta = seed / 9
    c, s = math.cos(theta), math.sin(theta)
    rotate = lambda row: [c*row[0]-s*row[1], s*row[0]+c*row[1], row[2]]
    lattice = [rotate([2+seed/8,0.25,0]), rotate([-0.5,3,0.75]),
               rotate([0.5,-0.25,(-1)**seed*4])] if periodic else None
    offset = [seed/3-1, -0.12345678901234567, 1.75] if periodic else None
    header = f'original synthetic {seed}' + (' %VIRTUAL' if virtual_style == 3 else '')
    text = chematic.write_exyz(symbols, points, header=header, lattice_angstrom=lattice,
        offset_angstrom=offset, virtual_atoms=flags)
    if crlf:
        text = text.replace('\n','\r\n')
    expected = {'elements':symbols, 'positions_angstrom':points, 'virtual_atoms':flags,
                'lattice_angstrom':lattice, 'offset_angstrom':offset}
    return text, expected


def check_geometry(text):
    g = chematic.parse_exyz(text)
    back = chematic.parse_exyz(emit(g))
    for key in g.keys() - {'source_text'}:
        assert back[key] == g[key], key
    assert g['source_text'] == text
    return g


def compare_periodic(executable, source, g):
    doc = json.loads(run_obabel(executable, source, 'cjson'))
    numbers = dict(zip(['C','H','N','O','Si','Cl'], [6,1,7,8,14,17]))
    assert doc['atoms']['elements']['number'] == [numbers[s] for s in g['elements']]
    cjson_error = assert_close(doc['atoms']['coords']['3d'], flat(g['positions_angstrom']), 1e-12)
    cell_error = assert_close(doc['unitCell']['cellVectors'], flat(g['lattice_angstrom']), 2e-12)
    output = run_obabel(executable, source, 'exyz')
    back = chematic.parse_exyz(output)
    assert back['elements'] == g['elements']
    rounded_error = assert_close(flat(back['positions_angstrom']), flat(g['positions_angstrom']), 5.00001e-6)
    rounded_cell = assert_close(flat(back['lattice_angstrom']), flat(g['lattice_angstrom']), 5.00001e-6)
    rounded_offset = assert_close(back['offset_angstrom'], g['offset_angstrom'], 5.00001e-6)
    # Open Babel 3.2.1 emits no per-atom VIRTUAL marker for these inputs.
    assert back['virtual_atoms'] == [False] * len(g['elements'])
    return cjson_error, cell_error, max(rounded_error, rounded_cell, rounded_offset), sum(g['virtual_atoms'])


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--obabel', default='/opt/homebrew/bin/obabel')
    parser.add_argument('--output', type=Path, default=ROOT/'validation/results/openbabel-exyz-2026-10-15.json')
    args = parser.parse_args()
    version = subprocess.check_output([args.obabel,'-V'], text=True).strip()
    assert version.split()[2] == '3.2.1', version
    records = []
    maxima = [0.0, 0.0, 0.0]
    lost_virtual_atoms = 0
    with tempfile.TemporaryDirectory(prefix='chematic-exyz-oracle-') as directory:
        source = Path(directory)/'geometry.exyz'
        for seed, periodic, virtual_style, crlf in itertools.product(range(12), [False,True], range(4), [False,True]):
            identifier = f's{seed}-p{int(periodic)}-v{virtual_style}-crlf{int(crlf)}'
            text, expected = fixture(seed,periodic,virtual_style,crlf)
            g = check_geometry(text)
            for key, value in expected.items():
                assert g[key] == value, (identifier, key)
            if periodic:
                source.write_bytes(text.encode())
                *errors, dropped = compare_periodic(args.obabel, source, g)
                maxima = [max(a,b) for a,b in zip(maxima, errors)]
                lost_virtual_atoms += dropped
            records.append({
                'id': identifier, 'sha256': digest(text.replace('\r\n','\n').encode()),
                'provenance': 'Original geometry generated by fixture() in scripts/exyz_openbabel_gate.py; 2026-10-15',
                'license': 'MIT OR Apache-2.0', 'redistribution': 'allowed; project-authored synthetic data',
                'minimization': '1-6 atoms; one skewed 3D cell only when periodic; flags and CRLF varied independently',
                'oracle': 'Open Babel 3.2.1 -iexyz FILE -ocjson/-oexyz -ab' if periodic else 'exact chematic structured roundtrip only',
                'expected': {'structured_roundtrip': 'exact', 'units': 'angstrom',
                    'cjson_tolerance': 2e-12 if periodic else None,
                    'openbabel_write_tolerance': 5.00001e-6 if periodic else None},
                'status': 'passed',
            })
        # Record the two observed Open Babel deltas without declaring parity.
        plain = '1\nplain\nC 0 0 0\n'
        source.write_text(plain)
        run = subprocess.run([args.obabel,'-iexyz',str(source),'-ocjson','-ab'], capture_output=True,text=True,check=True,timeout=10)
        assert '0 molecules converted' in run.stderr and 'Expecting unitcell' in run.stderr
        exyz_reader_plain = run.stderr.strip()
        run = subprocess.run([args.obabel,'-ixyz',str(source),'-oexyz','-ab'], capture_output=True,text=True,check=True,timeout=10)
        assert '1 molecule converted' in run.stderr
        synthesized = chematic.parse_exyz(run.stdout)
        assert synthesized['lattice_angstrom'] == [[1,0,0],[0,1,0],[0,0,1]]
    payload = json.dumps(records,sort_keys=True,separators=(',',':')).encode()
    report = {
        'schema_version':1, 'status':'local-source-candidate; uncommitted and unpushed',
        'oracle':version, 'scope':'192 original documents; 96 periodic Open Babel cross-reads; nonperiodic exact internal roundtrip',
        'records':records, 'normalized_records_sha256':digest(payload),
        'results': {'structured_roundtrips':len(records), 'periodic_openbabel_cases':96,
            'max_cjson_coordinate_error_angstrom':maxima[0], 'max_cjson_lattice_error_angstrom':maxima[1],
            'max_openbabel_written_error_angstrom':maxima[2], 'unexpected_mismatches':0,
            'virtual_atoms_omitted_by_openbabel_writer':lost_virtual_atoms},
        'observed_openbabel_deltas': {
            'virtual_markers':'Per-atom markers omitted by Open Babel 3.2.1 EXYZ writer; not part of compatibility claim',
            'nonperiodic_exyz_reader':exyz_reader_plain,
            'xyz_to_exyz_writer':'Adds %PBC and identity cell to nonperiodic XYZ input'},
        'limitations':['Single document only','No chemical bond/charge/spin/isotope inference',
            'No WASM/CLI/common-conversion bridge bindings','Normalized lattice determinant magnitude must exceed 1e-6'],
    }
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(report['results'],indent=2))


if __name__ == '__main__':
    main()
