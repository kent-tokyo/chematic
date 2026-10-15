#!/usr/bin/env python3
"""Original synthetic DMol3 .incoor against pinned Open Babel and ASE."""
import argparse
import hashlib
import itertools
import json
from pathlib import Path
import platform
import subprocess
import tempfile

import ase
from ase import Atoms
from ase.io import read, write
from ase.units import Bohr
import chematic

ROOT = Path(__file__).resolve().parents[1]
NATIVE_BOHR = 0.529177210544
OB_READ_BOHR = 0.529177249
OB_WRITE_INVERSE_BOHR = 1.889725989
ASE_BOHR = 0.5291772105638411


def digest(data):
    return hashlib.sha256(data).hexdigest()


def flatten(rows):
    return [v for row in rows for v in row]


def err(actual, expected):
    assert len(actual) == len(expected)
    return max((abs(a-b) for a, b in zip(actual, expected)), default=0.)


def emit(g):
    return chematic.write_dmol3(g['elements'], g['positions_bohr'], cell_bohr=g['cell_bohr'])


def recipe(seed, periodic, fortran, crlf):
    symbols = ['Si', 'Cl', 'Fe', 'H'][:seed % 4+1]
    positions = [[seed/7+i/3, -i/5+.125, (-1)**i*(seed+1)/11] for i in range(len(symbols))]
    cell = [[4+seed/8, 1, 0], [-1, 5, 2], [1, 0, (-1)**seed*6]] if periodic else None
    expected = dict(elements=symbols, positions_bohr=positions, cell_bohr=cell, periodicity=3 if periodic else 0,
                    positions_angstrom=[[v*NATIVE_BOHR for v in p] for p in positions],
                    cell_angstrom=None if cell is None else [[v*NATIVE_BOHR for v in p] for p in cell], bohr_angstrom=NATIVE_BOHR)
    # Independently authored input; the production writer is not the input oracle.
    num = lambda v: format(v, '.17e').replace('e', 'D' if fortran else 'e')
    lines = []
    if cell is not None:
        lines = ['$cell vectors', *(' '.join(map(num, row)) for row in cell)]
    lines += ['$coordinates', *(s+' '+' '.join(map(num, p)) for s, p in zip(symbols, positions)), '$end']
    text = '\n'.join(lines)+'\n'
    return expected, text.replace('\n', '\r\n') if crlf else text


def native_check(expected, text):
    g = chematic.parse_dmol3(text)
    for key, value in expected.items():
        assert g[key] == value, key
    assert g['source_text'] == text
    back = chematic.parse_dmol3(emit(g))
    for key in g.keys()-{'source_text'}:
        assert back[key] == g[key], key
    return g


def obabel(executable, path, output):
    r = subprocess.run([executable, '-idmol', str(path), '-o'+output, '-ab'], capture_output=True, text=True, check=True, timeout=10)
    assert '1 molecule converted' in r.stderr, r.stderr
    return r.stdout


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--obabel', default='/opt/homebrew/bin/obabel')
    parser.add_argument('--output', type=Path, default=ROOT/'validation/results/openbabel-dmol3-2026-10-15.json')
    args = parser.parse_args()
    version = subprocess.check_output([args.obabel, '-V'], text=True).strip()
    assert version.split()[2] == '3.2.1' and ase.__version__ == '3.27.0'
    assert Bohr == ASE_BOHR
    manifest_path = ROOT/'validation/openbabel_issues/dmol3/contract-fixtures.json'
    manifest = json.loads(manifest_path.read_text())
    for fixture in manifest['records']:
        assert digest(fixture['payload'].encode()) == fixture['sha256']
        native_check({k: v for k, v in fixture['expected'].items() if k not in ('units', 'tolerance')}, fixture['payload'])
    maxima = dict(ob_read_position_angstrom=0., ob_read_cell_angstrom=0., ase_read_position_angstrom=0., ase_read_cell_angstrom=0.,
                  ob_write_scaled_bohr=0., ase_write_raw_bohr=0., ob_raw_roundtrip_drift_bohr=0., ob_native_unit_delta_angstrom=0.)
    records = []
    common_count = 0
    with tempfile.TemporaryDirectory(prefix='chematic-dmol3-oracle-') as directory:
        path = Path(directory)/'input.incoor'
        ase_out = Path(directory)/'ase.incoor'
        for seed, periodic, fortran, crlf in itertools.product(range(16), (False, True), (False, True), (False, True)):
            expected, text = recipe(seed, periodic, fortran, crlf)
            g = native_check(expected, text)
            record = dict(id=f's{seed}-p{int(periodic)}-D{int(fortran)}-crlf{int(crlf)}', seed=seed, periodic=periodic, fortran=fortran, crlf=crlf,
                          provenance='Original synthetic recipe in scripts/validate_dmol3.py, authored 2026-10-15', license='MIT OR Apache-2.0', redistribution='allowed; project-authored data',
                          minimization='One to four atoms and optional skewed three-vector cell isolate raw units/precision/orientation',
                          oracle='Exact analytic raw Bohr; Open Babel 3.2.1 and ASE 3.27.0 with separately pinned unit constants; common E-exponent lane only',
                          expected=dict(status='exact native structured roundtrip and declared engine-normalized comparison', units='Bohr; explicit per-engine Angstrom projection', tolerance=0 if fortran else 2e-12),
                          sha256=digest(text.replace('\r\n', '\n').encode()), native_status='pass')
            if not fortran:
                path.write_bytes(text.encode())
                ob = json.loads(obabel(args.obabel, path, 'cjson'))
                assert ob['atoms']['elements']['number'] == [dict(Si=14, Cl=17, Fe=26, H=1)[s] for s in g['elements']]
                p = flatten(g['positions_bohr'])
                e = err(ob['atoms']['coords']['3d'], [v*OB_READ_BOHR for v in p]); assert e <= 2e-12
                maxima['ob_read_position_angstrom'] = max(maxima['ob_read_position_angstrom'], e)
                maxima['ob_native_unit_delta_angstrom'] = max(maxima['ob_native_unit_delta_angstrom'], err(ob['atoms']['coords']['3d'], flatten(g['positions_angstrom'])))
                if periodic:
                    e = err(ob['unitCell']['cellVectors'], [v*OB_READ_BOHR for v in flatten(g['cell_bohr'])]); assert e <= 2e-12
                    maxima['ob_read_cell_angstrom'] = max(maxima['ob_read_cell_angstrom'], e)
                else:
                    assert 'unitCell' not in ob
                cell = flatten(g['cell_bohr']) if periodic else [0.]*9
                raw = p+(cell if periodic else [])
                if periodic:
                    a = read(str(path), format='dmol-incoor')
                    assert a.get_chemical_symbols() == g['elements'] and a.pbc.tolist() == [True]*3
                    e = err(a.positions.ravel().tolist(), [v*ASE_BOHR for v in p]); assert e <= 2e-12
                    maxima['ase_read_position_angstrom'] = max(maxima['ase_read_position_angstrom'], e)
                    e = err(a.cell.array.ravel().tolist(), [v*ASE_BOHR for v in cell]); assert e <= 2e-12
                    maxima['ase_read_cell_angstrom'] = max(maxima['ase_read_cell_angstrom'], e)
                    write(str(ase_out), a, format='dmol-incoor')
                    back = chematic.parse_dmol3(ase_out.read_text())
                    assert back['elements'] == g['elements'] and back['periodicity'] == g['periodicity']
                    e = err(flatten(back['positions_bohr'])+flatten(back['cell_bohr']), raw); assert e <= 2e-12
                    maxima['ase_write_raw_bohr'] = max(maxima['ase_write_raw_bohr'], e)
                back = chematic.parse_dmol3(obabel(args.obabel, path, 'dmol'))
                assert back['elements'] == g['elements'] and back['periodicity'] == g['periodicity']
                actual = flatten(back['positions_bohr'])+(flatten(back['cell_bohr']) if periodic else [])
                e = err(actual, [v*OB_READ_BOHR*OB_WRITE_INVERSE_BOHR for v in raw]); assert e <= 2e-12
                maxima['ob_write_scaled_bohr'] = max(maxima['ob_write_scaled_bohr'], e)
                maxima['ob_raw_roundtrip_drift_bohr'] = max(maxima['ob_raw_roundtrip_drift_bohr'], err(actual, raw))
                record['openbabel_status'] = 'pass within declared normalization/tolerance'
                record['ase_status'] = 'pass' if periodic else 'outside scope: pinned ASE .incoor supports 3D periodicity only'
                common_count += 1
            records.append(record)
        text = '$cell vectors\n1 0 0\n0 1 0\n0 0 1\n$coordinates\nC 1.25D-1 0 0\n$end\n'
        assert chematic.parse_dmol3(text)['positions_bohr'][0][0] == 0.125
        path.write_text(text)
        ob_x = json.loads(obabel(args.obabel, path, 'cjson'))['atoms']['coords']['3d'][0]
        assert ob_x == 1.25*OB_READ_BOHR
        try:
            read(str(path), format='dmol-incoor')
        except ValueError:
            ase_d = 'rejected'
        else:
            raise AssertionError('pinned ASE D-exponent behavior changed')
        path.write_text('$coordinates\nC 0 0 0\n$end\n')
        try:
            read(str(path), format='dmol-incoor')
        except UnboundLocalError:
            ase_nonperiodic = 'raises UnboundLocalError without cell'
        else:
            raise AssertionError('pinned ASE nonperiodic-reader behavior changed')
        try:
            write(str(ase_out), Atoms('C', positions=[[0, 0, 0]], pbc=False), format='dmol-incoor')
        except ValueError as exc:
            assert 'PBC must be all true' in str(exc)
            ase_nonperiodic_write = 'rejects nonperiodic Atoms'
        else:
            raise AssertionError('pinned ASE nonperiodic-writer behavior changed')
    report = dict(schema_version=1, date='2026-10-15', feature='checked standalone DMol3 .incoor input/output',
                  runtime=dict(python=platform.python_version(), openbabel=version, ase=ase.__version__),
                  units=dict(native_bohr_angstrom=NATIVE_BOHR, native_codata='2022 central value', openbabel_read_bohr_angstrom=OB_READ_BOHR, openbabel_write_angstrom_bohr=OB_WRITE_INVERSE_BOHR, openbabel_raw_roundtrip_scale=OB_READ_BOHR*OB_WRITE_INVERSE_BOHR, ase_bohr_angstrom=ASE_BOHR),
                  oracle_settings=dict(openbabel_read='obabel -idmol INPUT -ocjson -ab', openbabel_write='obabel -idmol INPUT -odmol -ab', ase='ase.io.read/write(format="dmol-incoor"), default bohr=True; ASE read/write 3D only', read_tolerance_angstrom=2e-12, write_tolerance_bohr=2e-12, common_scope='canonical elements, E exponent, 0D/3D, LF/CRLF; raw-Bohr comparison adjusts only the explicitly pinned Open Babel read/write scale'),
                  provenance='Original synthetic MIT OR Apache-2.0 fixtures/recipes; public format description only, no Open Babel/ASE implementation copied',
                  fixture_manifest_sha256=digest(manifest_path.read_bytes()),
                  results=dict(native_generated_cases=len(records), native_manifest_cases=len(manifest['records']), openbabel_read_cases=common_count, openbabel_write_cases=common_count, ase_read_cases=common_count//2, ase_write_cases=common_count//2, unexpected_mismatches=0, maximum_absolute_errors=maxima),
                  diagnostics=dict(fortran_D=dict(native_raw_x_bohr=0.125, openbabel_x_angstrom=ob_x, ase_status=ase_d, input_sha256=digest(text.encode())), ase_nonperiodic_read=ase_nonperiodic, ase_nonperiodic_write=ase_nonperiodic_write), records=records,
                  scope='Standalone Bohr coordinate documents only; no output logs, electronic-state/charge inference, bond perception, calculation execution, physical validity or whole-library parity claim')
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2)+'\n')
    print(json.dumps(report['results'], indent=2))


if __name__ == '__main__':
    main()
