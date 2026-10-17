"""Independent ASCII STL geometry checks using trimesh and an ASE CIF lane."""
import ast
import hashlib
import io
import json
import math
from pathlib import Path

import numpy as np
import trimesh
import ase
from ase.io import read as ase_read
import chematic as ch

HERE = Path(__file__).parent
ROOT = HERE.parents[2]


def fixtures():
    tree = ast.parse((ROOT / 'crates/chematic-py/tests/test_crystal.py').read_text())
    return {n.targets[0].id: ast.literal_eval(n.value) for n in tree.body
            if isinstance(n, ast.Assign) and isinstance(n.targets[0], ast.Name)
            and n.targets[0].id in ('NACL_P1_CIF', 'SYMMETRIC_C2C_CIF')}


def main():
    # Compare byte hashes captured before extracting shared mesh preparation.
    baseline = json.loads((HERE / 'obj-baseline.json').read_text())
    for row in baseline:
        source = ch.PeriodicStructure(ch.Lattice.cubic(4.), [ch.Site([('C', 1.)], tuple(p)) for p in row['positions']])
        obj = ch.write_crystal_obj_checked(source, row['radii'], subdivisions=row['level'])
        assert hashlib.sha256(obj.text.encode('ascii')).hexdigest() == row['sha256']
    sources = []
    for count in (1, 2, 5):
        source = ch.PeriodicStructure(ch.Lattice.cubic(10.), [
            ch.Site([('C', 1.)], (i * 1.5 - .3, .2 * i - .2, i * .4)) for i in range(count)
        ])
        sources.append((f'authored-{count}', source))
    original = fixtures()
    cif_image_shifts = {}
    for name, text in original.items():
        source = ch.PeriodicStructure.from_cif(text)
        independent = ase_read(io.StringIO(text), format='cif')
        assert len(independent) == source.site_count()
        shifts = np.linalg.solve(independent.cell.array.T,
            (independent.positions - source.cartesian_positions()).T).T
        assert np.allclose(shifts, np.rint(shifts), atol=1e-12)
        cif_image_shifts[name] = np.rint(shifts).astype(int).tolist()
        sources.append((name, source))
    rows, totals = [], dict(spheres=0, triangles=0, vertices=0)
    for name, source in sources:
        radii = [.2 + .04 * i for i in range(source.site_count())]
        for level in range(5):
            for scale in (.25, 3., 12.5):
                result = ch.write_crystal_stl_checked(source, radii,
                    millimeters_per_angstrom=scale, subdivisions=level)
                mesh = trimesh.load_mesh(io.BytesIO(result.text.encode('ascii')), file_type='stl', process=True)
                assert len(mesh.faces) == result.triangle_count
                assert mesh.is_watertight and mesh.is_winding_consistent
                assert mesh.euler_number == 2 * source.site_count()
                assert np.isfinite(mesh.vertices).all() and mesh.volume > 0
                assert np.all(mesh.area_faces > 0)
                assert np.allclose(np.linalg.norm(mesh.face_normals, axis=1), 1., atol=1e-12)
                parts = mesh.split(only_watertight=False)
                assert len(parts) == source.site_count()
                centers = np.asarray(source.cartesian_positions()) * scale
                taken, errors = set(), []
                for part in parts:
                    assert part.euler_number == 2 and part.volume > 0
                    center = part.vertices.mean(axis=0)
                    site = int(np.argmin(np.linalg.norm(centers - center, axis=1)))
                    assert site not in taken
                    taken.add(site)
                    radius = radii[site] * scale
                    error = float(np.max(np.abs(np.linalg.norm(part.vertices - centers[site], axis=1) - radius)))
                    assert np.allclose(center, centers[site], atol=1e-10)
                    assert error < 1e-10
                    assert len(part.vertices) == 4 ** (level + 1) + 2
                    assert len(part.faces) == 8 * 4**level
                    assert 0 < part.volume <= 4 / 3 * math.pi * radius**3 * (1 + 1e-10)
                    expected_bounds = np.array([centers[site] - radius, centers[site] + radius])
                    assert np.allclose(part.bounds, expected_bounds, atol=1e-10)
                    errors.append(error)
                raw_normals = np.asarray([[float(x) for x in line.split()[2:]]
                    for line in result.text.splitlines() if line.lstrip().startswith('facet normal ')])
                # process=True may renumber faces; independently load unprocessed
                # triangles in serialization order before comparing raw normals.
                raw = trimesh.load_mesh(io.BytesIO(result.text.encode('ascii')), file_type='stl', process=False)
                cross = np.cross(raw.triangles[:, 1] - raw.triangles[:, 0], raw.triangles[:, 2] - raw.triangles[:, 0])
                computed = cross / np.linalg.norm(cross, axis=1)[:, None]
                assert np.allclose(raw_normals, computed, atol=1e-10)
                totals['spheres'] += result.sphere_count
                totals['triangles'] += result.triangle_count
                totals['vertices'] += len(mesh.vertices)
                rows.append(dict(source=name, level=level, mm_per_angstrom=scale,
                                 spheres=result.sphere_count, triangles=result.triangle_count,
                                 max_radial_error_mm=max(errors), bytes=len(result.text)))
                if name == 'NACL_P1_CIF' and level == 2 and scale == 3.:
                    (HERE / 'nacl-p1.cif').write_text(original[name])
                    (HERE / 'nacl-atom-spheres.stl').write_text(result.text, encoding='ascii')
    report = dict(trimesh=trimesh.__version__, ase=ase.__version__, cases=len(rows),
                  obj_byte_identical_cases=len(baseline), totals=totals,
                  ase_minus_native_lattice_image_shifts=cif_image_shifts,
                  failures=0, rows=rows)
    (HERE / 'oracle.json').write_text(json.dumps(report, indent=2) + '\n')
    print(json.dumps({k: v for k, v in report.items() if k != 'rows'}))


if __name__ == '__main__':
    main()
