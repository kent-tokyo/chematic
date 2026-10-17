"""Explicit geometry/radius contracts inspired by RDKit #9601."""
import itertools
import math

import pytest
import chematic as c


def test_single_sphere_and_radius_override():
    mol = c.from_smiles('C')
    a = mol.sasa_with_radii_checked([[0,0,0]], [1], probe_radius=0)
    b = mol.sasa_with_radii_checked([[0,0,0]], [2], probe_radius=0.5)
    assert a['total'] == pytest.approx(4*math.pi)
    assert b['per_atom'] == pytest.approx([25*math.pi])
    assert b['total'] == pytest.approx(sum(b['per_atom']))


@pytest.mark.parametrize('coords,radii,field', [([], [1], 'coordinates'), ([[0,0,0]], [], 'radii'), ([[0,0,0]]*2, [1], 'coordinates'), ([[0,0,0]], [1,2], 'radii')])
def test_counts(coords, radii, field):
    with pytest.raises(ValueError, match=field):
        c.from_smiles('C').sasa_with_radii_checked(coords, radii)


@pytest.mark.parametrize('radius', [0, -1, math.nan, math.inf, -math.inf])
def test_bad_radius(radius):
    with pytest.raises(ValueError, match='radius at atom 0'):
        c.from_smiles('C').sasa_with_radii_checked([[0,0,0]], [radius])


@pytest.mark.parametrize('probe', [-1, math.nan, math.inf, -math.inf])
def test_bad_probe(probe):
    with pytest.raises(ValueError, match='probe_radius'):
        c.from_smiles('C').sasa_with_radii_checked([[0,0,0]], [1], probe_radius=probe)


@pytest.mark.parametrize('points', [-1, 0, 31, 65537])
def test_points(points):
    with pytest.raises(ValueError, match='sphere_points'):
        c.from_smiles('C').sasa_with_radii_checked([[0,0,0]], [1], sphere_points=points)


@pytest.mark.parametrize('point', [32, 65536])
def test_point_boundaries(point):
    assert c.from_smiles('C').sasa_with_radii_checked([[0,0,0]], [1], probe_radius=0, sphere_points=point)['total'] == pytest.approx(4*math.pi)


@pytest.mark.parametrize('coordinate', [math.nan, math.inf, -math.inf])
def test_coordinate_finite(coordinate):
    with pytest.raises(ValueError, match='non-finite coordinate'):
        c.from_smiles('C').sasa_with_radii_checked([[0,coordinate,0]], [1])


@pytest.mark.parametrize('row', [[0,0], [0,0,0,0], ['x',0,0], [None,0,0]])
def test_coordinate_types(row):
    with pytest.raises((TypeError, ValueError)):
        c.from_smiles('C').sasa_with_radii_checked([row], [1])


@pytest.mark.parametrize('radii', [[None], ['1'], None, 1])
def test_radius_types(radii):
    with pytest.raises(TypeError):
        c.from_smiles('C').sasa_with_radii_checked([[0,0,0]], radii)


def test_empty_molecule():
    assert c.RWMol().GetMol().sasa_with_radii_checked([], []) == {'total':0.0, 'per_atom':[]}


def test_explicit_h_atoms_count():
    mol = c.from_smiles('[H]O[H]')
    with pytest.raises(ValueError, match='radii'):
        mol.sasa_with_radii_checked([[0,0,0],[20,0,0],[40,0,0]], [1])
    out = mol.sasa_with_radii_checked([[0,0,0],[20,0,0],[40,0,0]], [1,2,3], probe_radius=0)
    assert out['per_atom'] == pytest.approx([4*math.pi,16*math.pi,36*math.pi])


def test_duplicate_and_nested_spheres():
    mol = c.from_smiles('CC')
    with pytest.raises(ValueError, match='coincident'):
        mol.sasa_with_radii_checked([[0,0,0]]*2, [1,1])
    assert mol.sasa_with_radii_checked([[0,0,0]]*2, [1,2], probe_radius=0)['per_atom'] == pytest.approx([0,16*math.pi])


def test_translation_and_permutation():
    mol = c.from_smiles('CO')
    coords = [[0,0,0],[0,1,0]]
    out = mol.sasa_with_radii_checked(coords, [1,1.5], probe_radius=0, sphere_points=4096)
    translated = mol.sasa_with_radii_checked([[v+1e12 for v in row] for row in coords], [1,1.5], probe_radius=0, sphere_points=4096)
    assert translated == out
    permuted = mol.renumber_atoms([1,0]).sasa_with_radii_checked(coords[::-1], [1.5,1], probe_radius=0, sphere_points=4096)
    assert permuted['per_atom'] == out['per_atom'][::-1]


def test_multiple_occluders():
    out = c.from_smiles('CCC').sasa_with_radii_checked(
        [[0,0,0], [0,1,0], [0,2,0]], [1,1,1],
        probe_radius=0, sphere_points=4096,
    )
    assert out['per_atom'] == pytest.approx([3*math.pi, 2*math.pi, 3*math.pi])
    assert out['total'] == pytest.approx(8*math.pi)


def test_bounded_generators_and_original_exceptions():
    mol = c.from_smiles('C')
    consumed=[]
    def radii():
        for i in itertools.count():
            consumed.append(i)
            yield 1.0
    with pytest.raises(ValueError, match='radii'):
        mol.sasa_with_radii_checked([[0,0,0]], radii())
    assert consumed == [0,1]
    def broken():
        yield [0,0,0]
        raise RuntimeError('coordinate source')
    with pytest.raises(RuntimeError, match='coordinate source'):
        mol.sasa_with_radii_checked(broken(), [1])
    good = mol.sasa_with_radii_checked((row for row in [[0,0,0]]), iter([1]), probe_radius=0)
    assert good['total']==pytest.approx(4*math.pi)


@pytest.mark.parametrize('radius', [1e-300,1e300])
def test_area_range(radius):
    with pytest.raises(ValueError, match='range'):
        c.from_smiles('C').sasa_with_radii_checked([[0,0,0]], [radius], probe_radius=0)


def test_far_finite_coordinates():
    out=c.from_smiles('CC').sasa_with_radii_checked([[1e308,0,0],[-1e308,0,0]], [1,1], probe_radius=0)
    assert out['total']==pytest.approx(8*math.pi)


def test_total_overflow():
    with pytest.raises(ValueError, match='range'):
        c.from_smiles('CC').sasa_with_radii_checked([[0,0,0],[1e155,0,0]], [3e153,3e153], probe_radius=0)


def test_work_and_atom_caps():
    with pytest.raises(ValueError, match='50000000'):
        c.from_smiles('C'*40).sasa_with_radii_checked([[i*.01,0,0] for i in range(40)], [1]*40, probe_radius=0, sphere_points=65536)
    with pytest.raises(ValueError, match='2048 atoms'):
        c.from_smiles('C'*2049).sasa_with_radii_checked([], [])


def test_integer_and_keyword_contract():
    mol=c.from_smiles('C')
    with pytest.raises(TypeError):
        mol.sasa_with_radii_checked([[0,0,0]], [1], sphere_points=32.5)
    with pytest.raises(OverflowError):
        mol.sasa_with_radii_checked([[0,0,0]], [1], sphere_points=1<<100)
    with pytest.raises(TypeError):
        mol.sasa_with_radii_checked([[0,0,0]], [1], typo=True)
