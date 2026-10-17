"""Analytical and boundary contracts for explicit-charge dielectric fields."""

from concurrent.futures import ThreadPoolExecutor
from dataclasses import FrozenInstanceError
import math

import numpy as np
import pytest
from chematic.rdkit_compat import rdMIF

K = 1389.048632870505
DEPTHS = [0.0] * 7 + [0.4, 0.9, 1.4, 1.9, 2.6, 4.0, 4.0, 4.0, 4.0]


def test_sparse_model_and_strict_squared_threshold():
    field = rdMIF.CoulombDielectric([1.0], [[0, 0, 0]])
    assert repr(field) == "CoulombDielectric(sites=1)" and field.site_count == 1
    assert field(2, 0, 0, 4) == 0
    assert field(2, 0, 0, math.nextafter(4.0, math.inf)) == pytest.approx(K / 84)
    assert field(0, 0, 0, 1) == pytest.approx(K / 42)
    assert field(2, 0, 0, 0) == 0


@pytest.mark.parametrize("size", range(1, 16))
@pytest.mark.parametrize("alpha", [0.0, 0.4])
def test_every_neighbor_depth_and_soft_core_analytically(size, alpha):
    # Zero-charge neighbors contribute to the cavity, but not directly to energy.
    q = [1.0] + [0.0] * (size - 1)
    xyz = [[0.0, 0.0, 0.0]] * size
    field = rdMIF.CoulombDielectric(q, xyz, softcoreParam=alpha)
    d2 = alpha if alpha else 1.0
    contrast = (4 - 80) / (4 + 80)
    expected = (
        K
        / 4
        * (
            1 / math.sqrt(d2)
            + contrast / math.sqrt(d2 + 4 * DEPTHS[size - 1] * DEPTHS[size])
        )
    )
    assert field(0, 0, 0, 1) == pytest.approx(expected, rel=1e-14)


def test_all_sites_count_for_cavity_even_outside_energy_threshold():
    xyz = [[0.0, 0.0, 0.0]] + [[2.0, 0.0, 0.0]] * 7
    field = rdMIF.CoulombDielectric([1.0] + [0.0] * 7, xyz)
    expected = K / 4 * (1 + (-76 / 84) / math.sqrt(1 + 4 * 0.4 * 0.9))
    assert field(0, 0, 0, 1) == pytest.approx(expected)
    assert field(0, 0, 0, 100) == field(0, 0, 0, 1)
    assert field(0, 0, 0, 1) != pytest.approx(
        rdMIF.CoulombDielectric([1.0], [[0, 0, 0]])(0, 0, 0, 1)
    )


def test_neighbor_distance_four_is_strict():
    near = rdMIF.CoulombDielectric(
        [1.0] + [0.0] * 7, [[0, 0, 0]] + [[math.nextafter(4.0, 0), 0, 0]] * 7
    )
    edge = rdMIF.CoulombDielectric([1.0] + [0.0] * 7, [[0, 0, 0]] + [[4.0, 0, 0]] * 7)
    assert near(0, 0, 0, 1) > edge(0, 0, 0, 1)
    assert edge(0, 0, 0, 1) == pytest.approx(K / 42)


def test_equal_permittivities_reduce_to_scaled_vacuum_and_sign_total_abs():
    q = [1.0, -0.2, 0.4]
    xyz = [[0, 0, 0], [1, 0, 0], [0, 1, 0]]
    vacuum = rdMIF.Coulomb(q, xyz)(2, 3, 0, 100)
    for permittivity in (1.0, 4.0, 80.0, 1e6):
        field = rdMIF.CoulombDielectric(q, xyz, epsilon=permittivity, xi=permittivity)
        assert field(2, 3, 0, 100) == pytest.approx(vacuum / permittivity)
    signed = rdMIF.CoulombDielectric(q, xyz, probeCharge=-2)(2, 3, 0, 100)
    normal = rdMIF.CoulombDielectric(q, xyz)(2, 3, 0, 100)
    assert signed == pytest.approx(-2 * normal)
    assert rdMIF.CoulombDielectric(q, xyz, absVal=True)(2, 3, 0, 100) == -abs(normal)
    assert (
        rdMIF.CoulombDielectric([1.0, -1.0], [[0, 0, 0], [0, 0, 0]], absVal=True)(
            1, 0, 0, 2
        )
        == 0
    )


def test_source_copy_frozen_and_repeat_parallel_queries():
    q = np.array([1.0] + [0.2] * 12)
    xyz = np.zeros((13, 3))
    field = rdMIF.CoulombDielectric(q, xyz, epsilon=np.float64(80), xi=np.int64(4))
    points = [(0.0, 0.0, 0.0), (2.0, 0.0, 0.0), (8.0, 0.0, 0.0)] * 5
    expected = [field(*point, 100) for point in points]
    q[:] = 100
    xyz[:] = 100
    with ThreadPoolExecutor(max_workers=3) as pool:
        actual = list(pool.map(lambda p: field(*p, 100), points))
    assert actual == expected
    with pytest.raises(FrozenInstanceError):
        field._contrast = 0


@pytest.mark.parametrize("key", ["epsilon", "xi"])
@pytest.mark.parametrize("value", [True, "80", 1j, None])
def test_permittivity_types(key, value):
    with pytest.raises(TypeError):
        rdMIF.CoulombDielectric([1.0], [[0, 0, 0]], **{key: value})


@pytest.mark.parametrize("key", ["epsilon", "xi"])
@pytest.mark.parametrize("value", [0.0, -1.0, 0.5, math.inf, math.nan, 1e7, 10**1000])
def test_permittivity_bounds(key, value):
    with pytest.raises(ValueError):
        rdMIF.CoulombDielectric([1.0], [[0, 0, 0]], **{key: value})


def test_quadratic_cap_before_indexing_and_maximum_count():
    class Huge:
        def __len__(self):
            return 1025

        def __getitem__(self, i):
            raise AssertionError("must not index")

    with pytest.raises(ValueError, match="1024"):
        rdMIF.CoulombDielectric(Huge(), Huge())
    with pytest.raises(ValueError, match="1024"):
        rdMIF.CoulombDielectric([1.0], Huge())
    field = rdMIF.CoulombDielectric([0.0] * 1024, [[0, 0, 0]] * 1024)
    assert field.site_count == 1024 and field(0, 0, 0, 1) == 0


def test_shared_input_and_query_rejections_empty_and_extreme_values():
    for q, xyz, error in [
        ([True], [[0, 0, 0]], TypeError),
        (iter([1.0]), [[0, 0, 0]], TypeError),
        ([1.0], [[0, 0]], ValueError),
        ([1.0], [], ValueError),
        ([1.0], [[math.nan, 0, 0]], ValueError),
    ]:
        with pytest.raises(error):
            rdMIF.CoulombDielectric(q, xyz)
    for key, value in [
        ("probeCharge", True),
        ("absVal", 1),
        ("softcoreParam", 1e-7),
        ("cutoff", 0),
    ]:
        with pytest.raises((TypeError, ValueError)):
            rdMIF.CoulombDielectric([1.0], [[0, 0, 0]], **{key: value})
    field = rdMIF.CoulombDielectric(
        [1e6], [[-1e100, 0, 0]], probeCharge=1e6, epsilon=1e6, xi=1
    )
    assert field(1e100, 0, 0, 1e200) == 0
    assert math.isfinite(field(-1e100, 0, 0, 1))
    with pytest.raises(TypeError):
        field(True, 0, 0, 1)
    with pytest.raises(ValueError):
        field(0, 0, 0, -1)
    with pytest.raises(ValueError):
        field(math.inf, 0, 0, 1)
    empty = rdMIF.CoulombDielectric([], [], absVal=True)
    assert empty(0, 0, 0, 1) == 0 and math.copysign(1, empty(0, 0, 0, 1)) == -1
