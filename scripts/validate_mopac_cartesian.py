#!/usr/bin/env python3
"""Original synthetic MOPAC Cartesian contracts against Open Babel 3.2.1."""
import argparse
import hashlib
import itertools
import json
from pathlib import Path
import platform
import subprocess
import tempfile

import chematic

ROOT = Path(__file__).resolve().parents[1]
FIELDS = ('keywords', 'title', 'comment', 'labels', 'optimization_flags', 'translation_vectors_angstrom', 'translation_flags')
NUMBERS = {'Si': 14, 'Cl': 17, 'Fe': 26, 'H': 1}


def digest(data):
    return hashlib.sha256(data).hexdigest()


def emit(g):
    return chematic.write_mopac_cartesian(g['elements'], g['positions_angstrom'], **{k: g[k] for k in FIELDS})


def recipe(seed, dimension, flagged, labeled, crlf):
    symbols = ['Si', 'Cl', 'Fe', 'H'][:seed % 4+1]
    positions = [[seed/7+i/3, -i/5+0.125, (-1)**i*(seed+1)/11] for i in range(len(symbols))]
    vectors = [[4+seed/8, 1, 0], [-1, 5, 2], [1, 0, (-1)**seed*6]][:dimension]
    codes = ['-2', '-1', '0', '1', '2', 'T']
    atom_flags = [[codes[(seed+3*i+j) % 6] for j in range(3)] if flagged else None for i in range(len(symbols))]
    vector_flags = [['0', '1', '0'] if flagged else None for _ in vectors]
    labels = [f'site {i} alpha' if labeled else None for i in range(len(symbols))]
    expected = dict(keywords='PM7 XYZ CHARGE=1 MINI', title=f'original {seed} ', comment='Cartesian contract',
                    elements=symbols, positions_angstrom=positions, labels=labels, optimization_flags=atom_flags,
                    translation_vectors_angstrom=vectors, translation_flags=vector_flags, periodicity=dimension)
    # Input is authored independently of the production writer.
    lines = [expected['keywords'], expected['title'], expected['comment']]
    identities = [s+('('+label+')' if label is not None else '') for s, label in zip(symbols, labels)]
    for identity, point, flags in zip(identities+['Tv']*dimension, positions+vectors, atom_flags+vector_flags):
        fields = [identity]
        for j, value in enumerate(point):
            fields.append(format(value, '.17g'))
            if flags is not None:
                fields.append(flags[j])
        lines.append(' '.join(fields))
    text = '\n'.join(lines)+'\n\n'
    return expected, text.replace('\n', '\r\n') if crlf else text


def native_check(expected, text):
    g = chematic.parse_mopac_cartesian(text)
    assert g['source_text'] == text
    for key, value in expected.items():
        assert g[key] == value, key
    back = chematic.parse_mopac_cartesian(emit(g))
    for key in g.keys()-{'source_text'}:
        assert back[key] == g[key], key
    return g


def obabel(executable, path, fmt, *extra):
    r = subprocess.run([executable, '-imop', str(path), '-o'+fmt, '-ab', *extra], capture_output=True, text=True, check=True, timeout=10)
    assert '1 molecule converted' in r.stderr, r.stderr
    return r.stdout


def error(actual, expected):
    assert len(actual) == len(expected)
    return max((abs(a-b) for a, b in zip(actual, expected)), default=0.)


def flatten(rows):
    return [x for row in rows for x in row]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--obabel', default='/opt/homebrew/bin/obabel')
    parser.add_argument('--output', type=Path, default=ROOT/'validation/results/openbabel-mopac-cartesian-2026-10-15.json')
    args = parser.parse_args()
    version = subprocess.check_output([args.obabel, '-V'], text=True).strip()
    assert version.split()[2] == '3.2.1', version
    manifest_path = ROOT/'validation/openbabel_issues/mopac/contract-fixtures.json'
    manifest = json.loads(manifest_path.read_text())
    for fixture in manifest['records']:
        text = fixture['payload']
        assert digest(text.replace('\r\n', '\n').encode()) == fixture['sha256']
        expected = {k: v for k, v in fixture['expected'].items() if k not in ('units', 'tolerance')}
        native_check(expected, text)
    records = []
    ob_count = 0
    maxima = dict(read_positions_angstrom=0., read_cell_angstrom=0., write_positions_angstrom=0., write_cell_angstrom=0.)
    with tempfile.TemporaryDirectory(prefix='chematic-mopac-oracle-') as directory:
        path = Path(directory)/'input.mop'
        for seed, dimension, flagged, labeled, crlf in itertools.product(range(8), range(4), (False, True), (False, True), (False, True)):
            expected, text = recipe(seed, dimension, flagged, labeled, crlf)
            g = native_check(expected, text)
            common = flagged and dimension in (0, 3)
            record = dict(id=f's{seed}-d{dimension}-f{int(flagged)}-l{int(labeled)}-crlf{int(crlf)}',
                          normalized_input_sha256=digest(text.replace('\r\n', '\n').encode()),
                          provenance='Original synthetic recipe in scripts/validate_mopac_cartesian.py, authored 2026-10-15',
                          license='MIT OR Apache-2.0', redistribution='allowed; project-authored data',
                          minimization='One to four atoms and zero to three skewed vectors isolate coordinates/flags/labels/line endings',
                          oracle='Exact analytic expected values; common-lane Open Babel 3.2.1 commands and tolerances in oracle_settings',
                          expected=expected, native_status='pass', openbabel_common_geometry=common)
            if common:
                path.write_bytes(text.encode())
                ob = json.loads(obabel(args.obabel, path, 'cjson'))
                assert ob['atoms']['elements']['number'] == [NUMBERS[s] for s in g['elements']]
                e = error(ob['atoms']['coords']['3d'], flatten(g['positions_angstrom']))
                assert e <= 2e-12, e
                maxima['read_positions_angstrom'] = max(maxima['read_positions_angstrom'], e)
                if dimension == 3:
                    e = error(ob['unitCell']['cellVectors'], flatten(g['translation_vectors_angstrom']))
                    assert e <= 2e-12, e
                    maxima['read_cell_angstrom'] = max(maxima['read_cell_angstrom'], e)
                else:
                    assert 'unitCell' not in ob
                # Explicit unit-cell option; defaults deliberately replace headers and flags.
                output = obabel(args.obabel, path, 'mop', '-xu')
                back = chematic.parse_mopac_cartesian(output)
                assert back['elements'] == g['elements'] and back['periodicity'] == dimension
                e = error(flatten(back['positions_angstrom']), flatten(g['positions_angstrom']))
                assert e <= 5.00001e-6, e
                maxima['write_positions_angstrom'] = max(maxima['write_positions_angstrom'], e)
                e = error(flatten(back['translation_vectors_angstrom']), flatten(g['translation_vectors_angstrom']))
                assert e <= 5.00001e-6, e
                maxima['write_cell_angstrom'] = max(maxima['write_cell_angstrom'], e)
                assert back['optimization_flags'] == [['1']*3 for _ in g['elements']]
                assert back['translation_flags'] == [['1']*3 for _ in g['translation_vectors_angstrom']]
                assert back['labels'] == [None]*len(g['elements'])
                assert back['keywords'] == 'PUT KEYWORDS HERE' and back['title'] == str(path) and back['comment'] == ''
                assert ob['properties']['totalCharge'] == 0
                record['openbabel_geometry_status'] = 'pass'
                record['openbabel_known_deltas'] = ['five-decimal output', 'flags reset to 1', 'labels omitted', 'headers replaced', 'CHARGE keyword not projected']
                ob_count += 1
            records.append(record)
        diagnostics = []
        for name, body in [('unflagged', 'Si .125 -1.25 2.5'), ('empty-label', 'Si() .125 0 -1.25 1 2.5 0')]:
            text = 'PM7 XYZ\noriginal\n\n'+body+'\n\n'
            assert chematic.parse_mopac_cartesian(text)['elements'] == ['Si']
            path.write_text(text)
            r = subprocess.run([args.obabel, '-imop', str(path), '-ocjson', '-ab'], capture_output=True, text=True, timeout=10, check=True)
            assert '0 molecules converted' in r.stderr and not r.stdout
            diagnostics.append(dict(id=name, native_status='accepted-preserved', openbabel_status='rejected', stderr=r.stderr.strip(), normalized_input_sha256=digest(text.encode())))
        # Isolate Fortran exponent behavior instead of treating its interpretation as parity.
        text = 'PM7 XYZ\noriginal\n\nSi 1.25D-1 0 -1.25 1 2.5 0\n\n'
        path.write_text(text)
        ob = json.loads(obabel(args.obabel, path, 'cjson'))
        native = chematic.parse_mopac_cartesian(text)['positions_angstrom'][0][0]
        assert native == 0.125 and ob['atoms']['coords']['3d'][0] == 1.25
        diagnostics.append(dict(id='fortran-D-exponent', native_x_angstrom=native, openbabel_x_angstrom=ob['atoms']['coords']['3d'][0], normalized_input_sha256=digest(text.encode())))
    report = dict(schema_version=1, date='2026-10-15', feature='MOPAC Cartesian geometry input/output',
                  runtime=dict(python=platform.python_version(), openbabel=version),
                  scope='Original synthetic Cartesian geometry; no bond perception, charge/spin inference, calculation validity, quantum chemistry or whole-library comparison',
                  provenance='All fixtures/recipes project-authored MIT OR Apache-2.0; public manuals define grammar only; no upstream implementation copied',
                  oracle_settings=dict(read='obabel -imop INPUT -ocjson -ab', write='obabel -imop INPUT -omop -ab -xu', units='angstrom', read_tolerance=2e-12, write_tolerance=5.00001e-6, comparison_scope='flagged canonical-symbol input, 0D/3D, optional nonempty labels; all metadata/0D-3D exact native roundtrip'),
                  results=dict(native_generated_cases=len(records), native_manifest_cases=len(manifest['records']), openbabel_read_cases=ob_count, openbabel_write_cases=ob_count, unexpected_mismatches=0, maximum_absolute_errors=maxima),
                  manifest_sha256=digest(manifest_path.read_bytes()), diagnostics=diagnostics, records=records)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2)+'\n')
    print(json.dumps(report['results'], indent=2))


if __name__ == '__main__':
    main()
