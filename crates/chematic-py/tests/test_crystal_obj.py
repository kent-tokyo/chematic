"""Display geometry invariants and typed limits for crystal atom OBJ meshes."""

from collections import Counter
from dataclasses import FrozenInstanceError
import math

import numpy as np
import pytest
import chematic as ch
import chematic._crystal_obj as obj
from .test_crystal import NACL_P1_CIF, SYMMETRIC_C2C_CIF


def structure(count=1, fractional=(0.0, 0.0, 0.0)):
    return ch.PeriodicStructure(
        ch.Lattice.cubic(1.0), [ch.Site([("C", 1.0)], fractional)] * count
    )


def parse(text):
    vertices, normals, faces, objects = [], [], [], []
    for line in text.splitlines():
        fields = line.split()
        if fields[0] == "v":
            vertices.append([float(x) for x in fields[1:]])
        elif fields[0] == "vn":
            normals.append([float(x) for x in fields[1:]])
        elif fields[0] == "f":
            faces.append([tuple(int(x) - 1 for x in f.split("//")) for f in fields[1:]])
        elif fields[0] == "o":
            objects.append(fields[1])
        else:
            assert fields[0] == "#"
    return np.array(vertices), np.array(normals), faces, objects


@pytest.mark.parametrize("level", range(5))
def test_closed_outward_spheres_and_normals(level):
    source = ch.PeriodicStructure.from_cif(NACL_P1_CIF)
    centers = source.cartesian_positions()
    result = ch.write_crystal_obj_checked(source, [0.4, 0.7], subdivisions=level)
    vertices, normals, faces, objects = parse(result.text)
    per_site = 4 ** (level + 1) + 2
    assert result.sphere_count == 2
    assert result.vertex_count == len(vertices) == 2 * per_site
    assert result.face_count == len(faces) == 16 * 4**level
    assert len(normals) == per_site
    assert objects == ["site_1_Na", "site_2_Cl"]
    assert np.allclose(np.linalg.norm(normals, axis=1), 1.0, atol=1e-14)
    for site, radius in enumerate([0.4, 0.7]):
        selected = vertices[site * per_site : (site + 1) * per_site]
        assert np.allclose(np.mean(selected, axis=0), centers[site], atol=1e-14)
        assert np.allclose(
            np.linalg.norm(selected - centers[site], axis=1), radius, atol=1e-14
        )
        assert np.allclose((selected - centers[site]) / radius, normals, atol=1e-14)
    edges = Counter()
    for face in faces:
        indexes, normal_indexes = zip(*face)
        assert all(0 <= i < len(vertices) for i in indexes)
        assert all(0 <= i < len(normals) for i in normal_indexes)
        site = indexes[0] // per_site
        assert all(i // per_site == site for i in indexes)
        a, b, c = vertices[list(indexes)]
        assert np.dot(np.cross(b - a, c - a), (a + b + c) / 3 - centers[site]) > 0
        for i, j in zip(indexes, indexes[1:] + indexes[:1]):
            edges[tuple(sorted((i, j)))] += 1
    assert set(edges.values()) == {2}
    assert len(vertices) - len(edges) + len(faces) == 4


def test_no_wrap_and_no_mutation():
    source = structure(fractional=(2.5, -1.25, 0.125))
    before = source.to_poscar()
    result = ch.write_crystal_obj_checked(source, [0.25])
    vertices, _, _, _ = parse(result.text)
    assert np.allclose(vertices.mean(axis=0), [2.5, -1.25, 0.125])
    assert source.to_poscar() == before
    assert result == ch.write_crystal_obj_checked(source, [0.25])
    assert result.units == "angstrom"
    assert "periodic_boundary_conditions" in result.omitted_features
    with pytest.raises(FrozenInstanceError):
        result.face_count = 0


def test_expanded_cif_and_unexpanded_rejection():
    expanded = ch.PeriodicStructure.from_cif(SYMMETRIC_C2C_CIF)
    result = ch.write_crystal_obj_checked(expanded, [0.3] * expanded.site_count())
    assert result.sphere_count == expanded.site_count() == 4
    asymmetric = ch.PeriodicStructure.from_cif(SYMMETRIC_C2C_CIF, expand_symmetry=False)
    with pytest.raises(ValueError, match="unexpanded"):
        ch.write_crystal_obj_checked(asymmetric, [0.3])


@pytest.mark.parametrize(
    "species", [[("C", 0.5)], [("C", 0.5), ("N", 0.5)], [("C", 0.9999999999)]]
)
def test_partial_and_mixed_occupancy(species):
    source = ch.PeriodicStructure(
        ch.Lattice.cubic(3.0), [ch.Site(species, (0.0, 0.0, 0.0))]
    )
    with pytest.raises(ValueError, match="fully occupied"):
        ch.write_crystal_obj_checked(source, [0.3])


@pytest.mark.parametrize("value", [True, False, 1.0, "1", None])
def test_subdivision_type(value):
    with pytest.raises(TypeError):
        ch.write_crystal_obj_checked(structure(), [0.5], subdivisions=value)


@pytest.mark.parametrize("value", [-1, 5, 10**100])
def test_subdivision_range(value):
    with pytest.raises(ValueError):
        ch.write_crystal_obj_checked(structure(), [0.5], subdivisions=value)


@pytest.mark.parametrize("value", [None, "CIF", object()])
def test_structure_type(value):
    with pytest.raises(TypeError):
        ch.write_crystal_obj_checked(value, [0.5])


@pytest.mark.parametrize("value", [None, 0.5, "a", b"a", (x for x in [0.5])])
def test_radii_container_type(value):
    with pytest.raises(TypeError):
        ch.write_crystal_obj_checked(structure(), value)


@pytest.mark.parametrize("radii", [[], [0.5, 0.5]])
def test_radii_count(radii):
    with pytest.raises(ValueError, match="one value"):
        ch.write_crystal_obj_checked(structure(), radii)


@pytest.mark.parametrize(
    "value", [True, False, "0.5", None, 1j, object(), np.bool_(True)]
)
def test_radius_type(value):
    with pytest.raises(TypeError):
        ch.write_crystal_obj_checked(structure(), [value])


@pytest.mark.parametrize("value", [0.0, -1.0, math.nan, math.inf, -math.inf, 10**1000])
def test_radius_value(value):
    with pytest.raises(ValueError):
        ch.write_crystal_obj_checked(structure(), [value])


@pytest.mark.parametrize("values", [[1], (0.5,), [np.float64(0.5)], range(1, 2)])
def test_real_sequence_types(values):
    assert ch.write_crystal_obj_checked(structure(), values).sphere_count == 1


def test_numeric_collapse_is_refused():
    with pytest.raises(ValueError, match="coordinate scale"):
        ch.write_crystal_obj_checked(structure(fractional=(1e16, 0.0, 0.0)), [1.0])
    with pytest.raises(ValueError, match="overflow"):
        ch.write_crystal_obj_checked(structure(fractional=(1e308, 0.0, 0.0)), [1e308])


def test_site_and_resolution_caps():
    for count in [0, 1001]:
        with pytest.raises(ValueError, match="1..1000"):
            ch.write_crystal_obj_checked(structure(count), [0.5] * count)
    with pytest.raises(ValueError, match="budget"):
        ch.write_crystal_obj_checked(structure(98), [0.5] * 98, subdivisions=4)


def test_output_cap(monkeypatch):
    monkeypatch.setattr(obj, "_MAX_OUTPUT", 100)
    with pytest.raises(ValueError, match="output exceeds"):
        ch.write_crystal_obj_checked(structure(), [0.5])


def test_labels_cannot_inject_obj_commands():
    source = ch.PeriodicStructure(
        ch.Lattice.cubic(1.0),
        [ch.Site([("C", 1.0)], (0.0, 0.0, 0.0), "x\nf 1 2 3\nmtllib evil")],
    )
    text = ch.write_crystal_obj_checked(source, [0.5]).text
    assert "mtllib" not in text
    assert parse(text)[3] == ["site_1_C"]


def test_maximum_site_count_at_coarse_resolution():
    result = ch.write_crystal_obj_checked(
        structure(1000), [0.25] * 1000, subdivisions=0
    )
    assert result.sphere_count == 1000
    assert result.vertex_count == 6000 and result.face_count == 8000


def test_coincident_sites_remain_separate_objects():
    result = ch.write_crystal_obj_checked(structure(2), [0.25, 0.5], subdivisions=0)
    vertices, _, faces, objects = parse(result.text)
    assert objects == ["site_1_C", "site_2_C"]
    assert len(vertices) == 12 and len(faces) == 16
    assert all(all(index < 6 for index, _ in face) for face in faces[:8])
    assert all(all(index >= 6 for index, _ in face) for face in faces[8:])
