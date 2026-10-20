#!/usr/bin/env python3
"""Reproduce Open Babel #2845 on independent fixtures; pin the external engine."""
import argparse
import hashlib
import json
from pathlib import Path
import subprocess
import sys

import chematic

ROOT = Path(__file__).resolve().parents[2]
TYPES = {'charged-ions.pqr': ['Zn', 'Cl'], 'charged-records.pqr': ['N', 'C', 'HD', 'Zn', 'Cl']}


def sha256(data):
    return hashlib.sha256(data).hexdigest()


def rows(text):
    return [s for s in text.splitlines() if s.startswith(('ATOM  ', 'HETATM'))]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--obabel', default='/opt/homebrew/bin/obabel')
    parser.add_argument('--output', type=Path, default=ROOT / 'validation/results/openbabel-pqr-pdbqt-2026-10-14.json')
    args = parser.parse_args()
    version = subprocess.run([args.obabel, '-V'], capture_output=True, text=True, check=True, timeout=10).stdout.strip()
    if not version.startswith('Open Babel 3.2.1 '):
        raise RuntimeError(f'expected Open Babel 3.2.1, got {version}')
    fixtures = []
    for filename, atom_types in TYPES.items():
        path = ROOT / 'validation/openbabel_issues/pqr' / filename
        text = path.read_text()
        source = chematic.parse_pqr(text)
        report = chematic.pqr_to_pdbqt(text, atom_types, allow_radius_loss=True)
        output = rows(report['pdbqt'])
        expected = [a['charge'] for a in source['atoms']]
        actual = [float(r[70:76]) for r in output]
        assert len(output) == len(expected)
        assert max(abs(a-b) for a,b in zip(actual,expected)) <= 0.0005
        assert report['source_charges'] == expected
        assert report['source_radii'] == [a['radius'] for a in source['atoms']]
        entry = {
            'id': filename.removesuffix('.pqr'),
            'path': str(path.relative_to(ROOT)), 'sha256': sha256(text.encode()),
            'provenance': {'source': 'independently authored from APBS PQR field documentation', 'retrieved': '2026-10-14', 'issue': 'https://github.com/openbabel/openbabel/issues/2845', 'transformations': 'none; synthetic source records, no upstream attachment copied'},
            'license': 'MIT OR Apache-2.0',
            'redistribution': {'status': 'allowed', 'rationale': 'original project-authored fixture'},
            'minimization': 'two separated ions isolate charge retention; five records additionally exercise H, insertion codes, missing chain and noncontiguous serials',
            'oracle': {'name': 'APBS PQR documented fields and independent fixed-column PDBQT decoder', 'version': 'APBS 3.4.1 docs / AutoDock 4.2.6 manual', 'settings': 'source-order identity; columns 71-76 charges; no inferred bonds/types'},
            'expected': {'charges': expected, 'charge_unit': 'elementary charge', 'max_absolute_rounding_tolerance': 0.0005, 'radius_sidecar_unit': 'angstrom'},
            'chematic': {'output_charges': actual, 'source_radii': report['source_radii'], 'radii_omitted': report['radii_omitted'], 'max_charge_rounding_error': report['max_charge_rounding_error'], 'max_coordinate_rounding_error_angstrom': report['max_coordinate_rounding_error'], 'atom_types': atom_types, 'output_sha256': sha256(report['pdbqt'].encode())},
        }
        # The strict-column ion fixture is the upstream comparison control.
        # Do not confuse OB's failure on loose PQR fields with charge loss.
        if filename == 'charged-ions.pqr':
            comparisons = []
            for options in [[], ['-xr', '-xh', '-xp', '-xn']]:
                command = [args.obabel, '-ipqr', str(path), '-opdbqt', *options]
                proc = subprocess.run(command, capture_output=True, text=True, timeout=10, check=True)
                ob_rows = rows(proc.stdout)
                ob_charges = [float(r[70:76]) for r in ob_rows]
                assert len(ob_rows) == len(expected), proc.stdout
                assert [r[77:79].strip() for r in ob_rows] == atom_types
                assert all(float(r[30:38]) == a['x'] for r,a in zip(ob_rows,source['atoms']))
                assert ob_charges == [0.0, 0.0], 'upstream behavior changed; investigate before rebaseline'
                comparisons.append({'options': options, 'output_charges': ob_charges, 'stderr': proc.stderr.strip(), 'output_sha256': sha256(proc.stdout.encode())})
            entry['openbabel'] = comparisons
        fixtures.append(entry)
    module = Path(chematic.__file__).resolve()
    sources = ['crates/chematic-mol/src/pqr_pdbqt.rs', 'crates/chematic-py/src/formats/pqr_pdbqt.rs']
    record = {'schema_version': 1, 'status': 'local-source-candidate', 'issue': 'https://github.com/openbabel/openbabel/issues/2845', 'base': 'cac480fb23cd8adeec6e5e9bed1f18677f87b392 + preserved RDKit #9675 candidate; this patch contains only PQR changes', 'scope': 'rigid record adaptation with caller-supplied types; no docking quality, large protein, pH or flexible-ligand claim', 'runtime': {'python': sys.version.split()[0], 'openbabel': version, 'chematic_extension_sha256': sha256(module.read_bytes())}, 'source_sha256': {s:sha256((ROOT/s).read_bytes()) for s in sources}, 'fixtures': fixtures}
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(record, indent=2) + '\n')
    print(json.dumps({'fixtures':len(fixtures), 'upstream_charge_loss_reproduced':True, 'source_charge_retention':True, 'output':str(args.output)}))


if __name__ == '__main__':
    main()
