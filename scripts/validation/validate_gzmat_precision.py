#!/usr/bin/env python3
"""Check GZMAT scalar precision and geometry against RDKit and Open Babel."""
import argparse
import hashlib
import itertools
import json
from pathlib import Path
import subprocess
import sys

import chematic
import numpy as np
from rdkit import Chem, rdBase
from rdkit.Chem import rdMolTransforms as transform

ROOT = Path(__file__).resolve().parents[2]
FIXTURE = ROOT / 'validation/openbabel_issues/gzmat/precision-geometries.json'


def digest(data):
    return hashlib.sha256(data).hexdigest()


def reference(elements, positions, refs):
    c = Chem.Conformer(len(elements))
    for i, p in enumerate(positions):
        c.SetAtomPosition(i, p)
    result = []
    for i, r in enumerate(refs):
        values = []
        if r:
            values.append(transform.GetBondLength(c, i, r[0]))
        if len(r) > 1:
            values.append(transform.GetAngleDeg(c, i, r[0], r[1]))
        if len(r) > 2:
            values.append(transform.GetDihedralDeg(c, i, *r))
        result.append(values)
    return result


def readback(obabel, text):
    run = subprocess.run([obabel, '-igzmat', '-ocjson', '-ab'], input=text,
                         text=True, capture_output=True, timeout=10, check=True)
    assert '1 molecule converted' in run.stderr, run.stderr
    doc = json.loads(run.stdout)
    p = np.array(doc['atoms']['coords']['3d']).reshape(-1, 3)
    assert np.isfinite(p).all()
    return doc, p


def distance_matrix(p):
    return np.linalg.norm(p[:, None, :] - p[None, :, :], axis=2)


def signed_volume(p):
    return float(np.linalg.det(p[1:] - p[0]))


def scalar_errors(native, expected):
    errors = []
    for row, values in zip(native, expected):
        for k, (actual, want) in enumerate(zip(row['values'], values)):
            error = abs(actual - want) if k < 2 else abs((actual-want+180) % 360-180)
            errors.append(error)
    return errors


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--obabel', default='/opt/homebrew/bin/obabel')
    parser.add_argument('--output', type=Path, default=ROOT/'validation/results/openbabel-gzmat-precision-2026-10-15.json')
    args = parser.parse_args()
    assert rdBase.rdkitVersion == '2026.09.1'
    version = subprocess.check_output([args.obabel, '-V'], text=True).strip()
    assert version.split()[2] == '3.2.1', version
    manifest = json.loads(FIXTURE.read_text())
    assert digest(json.dumps(manifest['cases'], sort_keys=True, separators=(',', ':')).encode()) == manifest['sha256']
    # A proper rotation and translation must not change internal values/handedness.
    az, ax = np.deg2rad([43., 27.])
    rz = np.array([[np.cos(az), -np.sin(az), 0], [np.sin(az), np.cos(az), 0], [0, 0, 1]])
    rx = np.array([[1, 0, 0], [0, np.cos(ax), -np.sin(ax)], [0, np.sin(ax), np.cos(ax)]])
    rotation = rz @ rx
    checked = []
    max_error = {str(d): 0. for d in (8, 12, 17)}
    max_scalar_reference_error = 0.
    volumes_checked = 0
    exact_values = 0
    for case in manifest['cases']:
        elements = case['elements']
        n = len(elements)
        orders = list(itertools.islice(itertools.permutations(range(n)), 24))
        styles = [list([list(range(i-1, max(-1, i-4), -1)) for i in range(n)])]
        star = [list(range(min(i, 3))) for i in range(n)]
        if star != styles[0]:
            styles.append(star)
        for order in orders:
            e = [elements[i] for i in order]
            base = np.array(case['positions_angstrom'])[list(order)]
            for rotated in (False, True):
                p = base @ rotation.T + [2.3, -4.1, 1.7] if rotated else base
                for refs in styles:
                    before = repr((e, p.tolist(), refs))
                    native = chematic.cartesian_to_zmatrix(e, p.tolist(), refs)
                    expected = reference(e, p, refs)
                    error = max(scalar_errors(native, expected), default=0.)
                    max_scalar_reference_error = max(max_scalar_reference_error, error)
                    assert error < 1e-10, (case['id'], order, refs, native, expected)
                    for digits in (8, 12, 17):
                        text = chematic.write_gzmat(e, p.tolist(), refs, charge=0, multiplicity=1, significant_digits=digits)
                        scalars = [float(line.split('=')[1]) for line in text.splitlines() if '=' in line]
                        flat = [v for row in native for v in row['values']]
                        assert len(scalars) == len(flat)
                        if digits == 17:
                            assert scalars == flat
                            exact_values += len(scalars)
                        for value, original in zip(scalars, flat):
                            assert abs(value-original) <= max(abs(original), 1.) * 5 * 10**(-digits)
                        doc, q = readback(args.obabel, text)
                        assert doc['atoms']['elements']['number'] == [Chem.GetPeriodicTable().GetAtomicNumber(s) for s in e]
                        assert doc['properties']['totalCharge'] == 0 and doc['properties']['totalSpinMultiplicity'] == 1
                        error = float(np.max(np.abs(distance_matrix(p)-distance_matrix(q))))
                        max_error[str(digits)] = max(max_error[str(digits)], error)
                        assert error < manifest['expected']['geometry_tolerance_angstrom'][str(digits)], (case['id'], order, digits, error)
                        for indices in itertools.combinations(range(n), 4):
                            vp, vq = signed_volume(p[list(indices)]), signed_volume(q[list(indices)])
                            if abs(vp) > 1e-8:
                                assert vp * vq > 0, (case['id'], order, digits, vp, vq)
                                volumes_checked += 1
                        assert repr((e, p.tolist(), refs)) == before
                        checked.append([case['id'], order, rotated, refs, digits, digest(text.encode())])
    # Independently authored probe reproduces precision-option insensitivity.
    probe = manifest['cases'][3]
    xyz = str(len(probe['elements']))+'\noriginal precision probe\n'+''.join(s+' '+ ' '.join(format(v,'.17g') for v in p)+'\n' for s,p in zip(probe['elements'],probe['positions_angstrom']))
    def upstream(extra):
        return subprocess.run([args.obabel, '-ixyz', '-ogzmat', *extra], input=xyz, text=True, capture_output=True, check=True).stdout
    default, requested = upstream([]), upstream(['-x8'])
    assert default == requested
    upstream_r2 = float(next(line.split('=')[1] for line in requested.splitlines() if line.startswith('r2=')))
    assert abs(upstream_r2 - probe['positions_angstrom'][1][0]) > 1e-5
    paths = ['crates/chematic-mol/src/gzmat.rs', 'crates/chematic-mol/src/gzmat/geometry.rs', 'crates/chematic-mol/src/gzmat/writer.rs', 'crates/chematic-py/src/formats/gzmat.rs', 'crates/chematic-mol/tests/gzmat_precision.rs', 'crates/chematic-py/tests/test_gzmat_precision.py', 'scripts/validation/validate_gzmat_precision.py']
    result = {'schema_version': 1, 'status': 'local-source-candidate',
        'issue': 'https://github.com/openbabel/openbabel/issues/2794',
        'base': 'Preserved source candidates through RDKit meso-group reduction; no main integration',
        'scope': 'Caller-referenced geometry conversion and Gaussian-style Z-matrix serialization; no general Gaussian job compatibility claim',
        'runtime': {'python':sys.version.split()[0], 'numpy':np.__version__, 'rdkit':rdBase.rdkitVersion, 'openbabel':version, 'openbabel_sha256':digest(Path(args.obabel).read_bytes()), 'openbabel_readback':'-igzmat -ocjson -ab; preserve indexed elements, disable bonding', 'chematic_extension_sha256':digest(Path(chematic.__file__).resolve().read_bytes())},
        'source_sha256': {p:digest((ROOT/p).read_bytes()) for p in paths},
        'fixture': {'path':str(FIXTURE.relative_to(ROOT)), 'manifest_sha256':digest(FIXTURE.read_bytes()), **manifest},
        'results': {'cases':len(checked), 'fixtures':len(manifest['cases']), 'significant_digits':[8,12,17], 'computed_scalars_roundtrip_exact_at_17_digits':exact_values, 'signed_tetrahedra_verified':volumes_checked, 'max_scalar_reference_error_angstrom_or_degree':max_scalar_reference_error, 'max_pairwise_distance_error_angstrom':max_error, 'geometry_mismatches':0, 'element_mismatches':0, 'handedness_mismatches':0, 'input_mutations':0, 'case_sha256':digest(json.dumps(checked, separators=(',', ':')).encode())},
        'upstream_diagnostic': {'default_output_equals_x8_output':True, 'requested_r2':upstream_r2, 'source_r2':probe['positions_angstrom'][1][0], 'gzmat_output':requested},
        'limitations': 'Geometry-only, not topology. No isotope/ghost/dummy/periodic inputs, automatic reference choice, Gaussian input reader, ORCA internal writer, energy or job execution. Charge/spin are explicit, not inferred or checked for physical electron parity. Absolute Cartesian frame is discarded. Near-linear angles and reference planes (sine <=1e-10), coincident references, overflow and rounding-induced invalid frames fail. 17 digits preserves computed f64 scalars; geometry conversion itself is floating point, not bit-exact Cartesian preservation. Open Babel readback is measured; Gaussian/ORCA execution is untested.'}
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps(result['results']))


if __name__ == '__main__':
    main()
