"""ASCII facet grammar, physical display scaling and checked mesh limits."""
from dataclasses import FrozenInstanceError
import math

import numpy as np
import pytest
import chematic as ch
import chematic._crystal_stl as stl
from .test_crystal import NACL_P1_CIF, SYMMETRIC_C2C_CIF
from .test_crystal_obj import structure


def parse(text):
    lines = text.splitlines()
    assert lines[0] == "solid chematic_atom_spheres_mm"
    assert lines[-1] == "endsolid chematic_atom_spheres_mm"
    assert (len(lines) - 2) % 7 == 0
    normals, triangles = [], []
    for offset in range(1, len(lines) - 1, 7):
        head = lines[offset].split()
        assert head[:2] == ["facet", "normal"] and len(head) == 5
        normals.append([float(v) for v in head[2:]])
        assert lines[offset + 1].split() == ["outer", "loop"]
        face = []
        for row in lines[offset + 2:offset + 5]:
            fields = row.split()
            assert fields[0] == "vertex" and len(fields) == 4
            face.append([float(x) for x in fields[1:]])
        triangles.append(face)
        assert lines[offset + 5].split() == ["endloop"]
        assert lines[offset + 6].split() == ["endfacet"]
    return np.asarray(normals), np.asarray(triangles)


@pytest.mark.parametrize("level", range(5))
def test_scaled_spheres_have_outward_unit_normals_and_expected_facets(level):
    source = ch.PeriodicStructure.from_cif(NACL_P1_CIF)
    before = source.to_poscar()
    result = ch.write_crystal_stl_checked(source, [.4, .7], millimeters_per_angstrom=3., subdivisions=level)
    normals, triangles = parse(result.text)
    assert result.sphere_count == 2
    assert result.triangle_count == len(triangles) == 16 * 4**level
    assert np.isfinite(triangles).all() and np.isfinite(normals).all()
    assert np.allclose(np.linalg.norm(normals, axis=1), 1., atol=1e-14)
    for site, radius in enumerate([.4, .7]):
        faces = triangles[site * 8 * 4**level:(site + 1) * 8 * 4**level]
        center = source.cartesian_positions()[site] * 3.
        assert np.allclose(np.linalg.norm(faces - center, axis=2), radius * 3., atol=1e-13)
        cross = np.cross(faces[:, 1] - faces[:, 0], faces[:, 2] - faces[:, 0])
        assert np.all(np.einsum('ij,ij->i', cross, faces.mean(axis=1) - center) > 0)
        expected = cross / np.linalg.norm(cross, axis=1)[:, None]
        assert np.allclose(normals[site * len(faces):(site + 1) * len(faces)], expected, atol=1e-12)
    assert source.to_poscar() == before
    assert result.units == "millimeter" and result.millimeters_per_angstrom == 3.
    assert "site_and_sphere_component_identity" in result.omitted_features
    with pytest.raises(FrozenInstanceError):
        result.triangle_count = 1


def test_scale_is_required_and_is_uniform_without_wrapping():
    source = structure(fractional=(2.5, -1.25, .125))
    with pytest.raises(TypeError):
        ch.write_crystal_stl_checked(source, [.25])
    result = ch.write_crystal_stl_checked(source, [.25], millimeters_per_angstrom=10., subdivisions=0)
    _, triangles = parse(result.text)
    points = np.unique(triangles.reshape(-1, 3), axis=0)
    assert np.allclose(points.mean(axis=0), [25., -12.5, 1.25])
    assert np.allclose(np.linalg.norm(points - points.mean(axis=0), axis=1), 2.5)
    assert result == ch.write_crystal_stl_checked(source, [.25], millimeters_per_angstrom=10., subdivisions=0)


@pytest.mark.parametrize("value", [True, False, None, "1", object()])
def test_invalid_scale_type(value):
    with pytest.raises(TypeError):
        ch.write_crystal_stl_checked(structure(), [.5], millimeters_per_angstrom=value)


@pytest.mark.parametrize("value", [0., -1., math.nan, math.inf, -math.inf, 10**1000])
def test_invalid_scale_value(value):
    with pytest.raises(ValueError):
        ch.write_crystal_stl_checked(structure(), [.5], millimeters_per_angstrom=value)


@pytest.mark.parametrize("value", [True, 2., "2", None])
def test_subdivision_type(value):
    with pytest.raises(TypeError):
        ch.write_crystal_stl_checked(structure(), [.5], millimeters_per_angstrom=1., subdivisions=value)


@pytest.mark.parametrize("value", [-1, 5])
def test_subdivision_range(value):
    with pytest.raises(ValueError):
        ch.write_crystal_stl_checked(structure(), [.5], millimeters_per_angstrom=1., subdivisions=value)


def test_site_face_and_output_caps(monkeypatch):
    for count in [0, 1001]:
        with pytest.raises(ValueError, match="1..1000"):
            ch.write_crystal_stl_checked(structure(count), [.5] * count, millimeters_per_angstrom=1.)
    with pytest.raises(ValueError, match="triangle budget"):
        ch.write_crystal_stl_checked(structure(25), [.5] * 25, millimeters_per_angstrom=1., subdivisions=4)
    monkeypatch.setattr(stl, '_MAX_OUTPUT', 100)
    with pytest.raises(ValueError, match="output exceeds"):
        ch.write_crystal_stl_checked(structure(), [.5], millimeters_per_angstrom=1.)


@pytest.mark.parametrize("species", [[('C', .5)], [('C', .5), ('N', .5)]])
def test_unsupported_occupancy(species):
    source = ch.PeriodicStructure(ch.Lattice.cubic(1.), [ch.Site(species, (0., 0., 0.))])
    with pytest.raises(ValueError, match="fully occupied"):
        ch.write_crystal_stl_checked(source, [.5], millimeters_per_angstrom=1.)


def test_unexpanded_symmetry_is_refused():
    source = ch.PeriodicStructure.from_cif(SYMMETRIC_C2C_CIF, expand_symmetry=False)
    with pytest.raises(ValueError, match="unexpanded"):
        ch.write_crystal_stl_checked(source, [.5], millimeters_per_angstrom=1.)


@pytest.mark.parametrize("radii,error", [([], ValueError), ([0.], ValueError), ([True], TypeError), (iter([.5]), TypeError)])
def test_shared_radius_contract(radii, error):
    with pytest.raises(error):
        ch.write_crystal_stl_checked(structure(), radii, millimeters_per_angstrom=1.)


def test_scaled_overflow_underflow_and_coordinate_collapse_are_refused():
    for source, radius, scale in [(structure(), 1e200, 1e200), (structure(), 1e-200, 1e-200),
                                  (structure(fractional=(1e16, 0., 0.)), 1., 1.)]:
        with pytest.raises(ValueError):
            ch.write_crystal_stl_checked(source, [radius], millimeters_per_angstrom=scale)


def test_wrong_structure_and_hostile_labels():
    with pytest.raises(TypeError):
        ch.write_crystal_stl_checked(None, [.5], millimeters_per_angstrom=1.)
    source = ch.PeriodicStructure(ch.Lattice.cubic(1.), [ch.Site([('C', 1.)], (0., 0., 0.), 'bad\nendsolid evil')])
    result = ch.write_crystal_stl_checked(source, [.5], millimeters_per_angstrom=1.)
    assert 'evil' not in result.text
    assert parse(result.text)[1].shape == (128, 3, 3)


def test_maximum_sites_at_coarse_resolution():
    result = ch.write_crystal_stl_checked(structure(1000), [.25] * 1000, millimeters_per_angstrom=1., subdivisions=0)
    assert result.sphere_count == 1000 and result.triangle_count == 8000
