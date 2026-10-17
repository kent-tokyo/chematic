#!/usr/bin/env python3
"""Independent trimesh OBJ readback and ASE CIF-to-mesh geometry checks."""
from __future__ import annotations

import argparse
import hashlib
import io
import json
from pathlib import Path

import ase
from ase import Atoms
from ase.io import read, write
import numpy as np
import trimesh
import chematic as ch


def inspect(source, radii, level, expected_centers):
    result = ch.write_crystal_obj_checked(source, radii, subdivisions=level)
    scene = trimesh.load(
        io.StringIO(result.text),
        file_type="obj",
        force="scene",
        process=False,
        group_material=False,
        skip_materials=True,
    )
    assert len(scene.geometry) == source.site_count()
    max_center_error = max_radius_error = 0.0
    vertices = faces = 0
    for i, radius in enumerate(radii):
        name = f"site_{i + 1}_{source.sites[i].species[0][0]}"
        mesh = scene.geometry[name]
        assert mesh.is_watertight and mesh.is_winding_consistent
        assert mesh.euler_number == 2
        assert mesh.volume > 0
        assert mesh.faces.shape == (8 * 4**level, 3)
        assert mesh.vertices.shape == (4 ** (level + 1) + 2, 3)
        center_error = float(np.max(np.abs(mesh.center_mass - expected_centers[i])))
        distances = np.linalg.norm(mesh.vertices - expected_centers[i], axis=1)
        radius_error = float(np.max(np.abs(distances - radius)))
        assert center_error < 1e-10 and radius_error < 1e-10
        unit_directions = (mesh.vertices - expected_centers[i]) / radius
        assert np.allclose(mesh.vertex_normals, unit_directions, atol=1e-10)
        assert np.all(
            np.einsum(
                "ij,ij->i",
                mesh.face_normals,
                mesh.triangles_center - expected_centers[i],
            )
            > 0
        )
        max_center_error = max(max_center_error, center_error)
        max_radius_error = max(max_radius_error, radius_error)
        vertices += len(mesh.vertices)
        faces += len(mesh.faces)
    assert vertices == result.vertex_count and faces == result.face_count
    return {
        "spheres": len(radii),
        "subdivisions": level,
        "vertices": vertices,
        "faces": faces,
        "max_center_error_angstrom": max_center_error,
        "max_radius_error_angstrom": max_radius_error,
    }, result


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    rng = np.random.default_rng(2680)
    cells = [
        np.diag([3.0, 5.0, 7.0]),
        np.array([[4.0, 0.0, 0.0], [1.0, 5.0, 0.0], [-0.7, 0.4, 6.0]]),
        np.array([[2.5, 0.2, 0.1], [0.1, 4.0, 0.3], [0.3, -0.5, 6.5]]),
    ]
    records = []
    symbols = ["C", "O", "Si", "Na", "Cl", "Fe", "H", "S", "Au", "F", "B", "N", "Xe"]
    for cell_index, cell in enumerate(cells):
        for count in [1, 2, 7, 13]:
            fractional = rng.uniform(-0.5, 1.5, size=(count, 3))
            radii = rng.uniform(0.15, 0.85, size=count).tolist()
            source = ch.PeriodicStructure(
                ch.Lattice.from_matrix(cell.tolist()),
                [
                    ch.Site([(symbol, 1.0)], tuple(point))
                    for symbol, point in zip(symbols[:count], fractional)
                ],
            )
            expected = fractional @ cell
            for level in range(5):
                record, _ = inspect(source, radii, level, expected)
                records.append({"kind": "explicit_cell", "cell": cell_index, **record})
            # Independently author/read a P1 CIF with ASE, then pass that text
            # through chematic's native reader before generating OBJ geometry.
            atoms = Atoms(
                symbols[:count], scaled_positions=fractional, cell=cell, pbc=True
            )
            buffer = io.BytesIO()
            write(buffer, atoms, format="cif")
            cif = buffer.getvalue().decode("utf-8")
            expected_atoms = read(io.StringIO(cif), format="cif")
            parsed = ch.PeriodicStructure.from_cif(cif)
            record, _ = inspect(parsed, radii, 2, expected_atoms.positions)
            records.append({"kind": "ase_cif", "cell": cell_index, **record})
    all_symbols = [ase.data.chemical_symbols[i] for i in range(1, 119)]
    source = ch.PeriodicStructure(
        ch.Lattice.cubic(20.0),
        [
            ch.Site([(symbol, 1.0)], (i / 118.0, 0.25, 0.75))
            for i, symbol in enumerate(all_symbols)
        ],
    )
    record, _ = inspect(source, [0.2] * 118, 2, source.cartesian_positions())
    records.append({"kind": "118_elements", **record})
    report = {
        "trimesh_version": trimesh.__version__,
        "ase_version": ase.__version__,
        "conditions": len(records),
        "spheres": sum(r["spheres"] for r in records),
        "vertices": sum(r["vertices"] for r in records),
        "faces": sum(r["faces"] for r in records),
        "max_center_error_angstrom": max(
            r["max_center_error_angstrom"] for r in records
        ),
        "max_radius_error_angstrom": max(
            r["max_radius_error_angstrom"] for r in records
        ),
        "extension_sha256": hashlib.sha256(
            Path(ch.chematic.__file__).read_bytes()
        ).hexdigest(),
        "records": records,
    }
    args.output.write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps({k: v for k, v in report.items() if k != "records"}, indent=2))


if __name__ == "__main__":
    main()
