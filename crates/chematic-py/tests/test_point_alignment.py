"""Point input/transform contracts motivated by RDKit #9601."""

from dataclasses import dataclass
import math

import numpy as np
import pytest
from chematic.rdkit_compat import rdAlignment

REFERENCE = np.array(
    [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 2.0, 0.0], [0.0, 0.0, 3.0]]
)


@dataclass(frozen=True)
class Point:
    x: float
    y: float
    z: float


def transformed(points, matrix):
    return np.asarray(points) @ matrix[:3, :3].T + matrix[:3, 3]


@pytest.mark.parametrize(
    "kind", ["list", "tuple", "array", "point", "strided", "readonly", "endian"]
)
def test_container_inputs_and_probe_to_reference_direction(kind):
    probe = REFERENCE @ np.array(
        [[0.0, -1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, 1.0]]
    ) + [3.0, -7.0, 2.0]
    before = probe.copy()
    if kind == "list":
        points = probe.tolist()
    elif kind == "tuple":
        points = tuple(map(tuple, probe))
    elif kind == "point":
        points = [Point(*row) for row in probe]
    elif kind == "strided":
        backing = np.zeros((4, 6))
        backing[:, ::2] = probe
        points = backing[:, ::2]
    elif kind == "readonly":
        points = probe.copy()
        points.flags.writeable = False
    elif kind == "endian":
        points = probe.astype(">f8")
    else:
        points = probe
    ssd, matrix = rdAlignment.GetAlignmentTransform(REFERENCE, points)
    assert ssd == pytest.approx(0.0, abs=1e-25)
    np.testing.assert_allclose(transformed(probe, matrix), REFERENCE, atol=1e-14)
    np.testing.assert_array_equal(probe, before)
    assert matrix.shape == (4, 4) and matrix.dtype == np.float64
    np.testing.assert_array_equal(matrix[3], [0.0, 0.0, 0.0, 1.0])
    assert np.linalg.det(matrix[:3, :3]) == pytest.approx(1.0)


def test_weighted_ssd_is_not_rmsd_and_matches_transformed_residuals():
    probe = REFERENCE.copy()
    probe[-1] += [0.3, -0.4, 1.0]
    weights = np.array([1.0, 2.0, 5.0, 9.0])
    ssd, matrix = rdAlignment.GetAlignmentTransform(REFERENCE, probe, weights)
    errors = transformed(probe, matrix) - REFERENCE
    assert ssd == pytest.approx(float(np.sum(weights[:, None] * errors**2)))
    ssd_scaled, same = rdAlignment.GetAlignmentTransform(REFERENCE, probe, weights * 7)
    assert ssd_scaled == pytest.approx(ssd * 7)
    np.testing.assert_allclose(same, matrix, atol=1e-13)


@pytest.mark.parametrize(
    "weights", [None, [], (), np.array([], dtype=float), [1.0, 1.0, 1.0, 1.0]]
)
def test_default_or_empty_weights(weights):
    ssd, matrix = rdAlignment.GetAlignmentTransform(REFERENCE, REFERENCE + 2, weights)
    assert ssd == pytest.approx(0.0, abs=1e-25)
    np.testing.assert_allclose(
        transformed(REFERENCE + 2, matrix), REFERENCE, atol=1e-14
    )


def test_reflect_requires_an_improper_transform():
    mirror = REFERENCE * [-1.0, 1.0, 1.0]
    proper, _ = rdAlignment.GetAlignmentTransform(REFERENCE, mirror)
    reflected, matrix = rdAlignment.GetAlignmentTransform(
        REFERENCE, mirror, reflect=True
    )
    assert proper > 0.1
    assert reflected == pytest.approx(0.0, abs=1e-25)
    assert np.linalg.det(matrix[:3, :3]) == pytest.approx(-1.0)
    np.testing.assert_allclose(transformed(mirror, matrix), REFERENCE, atol=1e-14)


@pytest.mark.parametrize("n", [1, 2, 3, 4])
@pytest.mark.parametrize("reflect", [False, True])
def test_degenerate_point_clouds(n, reflect):
    reference = np.array([[float(i), 0.0, 0.0] for i in range(n)])
    probe = reference + [1.0, 2.0, 3.0]
    ssd, matrix = rdAlignment.GetAlignmentTransform(reference, probe, reflect=reflect)
    assert ssd == pytest.approx(0.0, abs=1e-25)
    np.testing.assert_allclose(transformed(probe, matrix), reference, atol=1e-14)
    assert np.linalg.det(matrix[:3, :3]) == pytest.approx(-1.0 if reflect else 1.0)


@pytest.mark.parametrize(
    "bad,exc",
    [
        ([], ValueError),
        ([[0.0, 0.0]], ValueError),
        ([[0.0, 0.0, 0.0, 0.0]], ValueError),
        (np.empty((1, 0)), ValueError),
        (np.zeros((1, 3, 1)), ValueError),
        ([[math.nan, 0.0, 0.0]], ValueError),
        ([[math.inf, 0.0, 0.0]], ValueError),
        ([[1e101, 0.0, 0.0]], ValueError),
        ([[10**400, 0.0, 0.0]], ValueError),
        (None, TypeError),
        ("abc", TypeError),
        (iter([[0.0, 0.0, 0.0]]), TypeError),
        ([[True, 0.0, 0.0]], TypeError),
        ([["0", 0.0, 0.0]], TypeError),
        ([[1j, 0.0, 0.0]], TypeError),
        ([object()], TypeError),
        (np.zeros((1, 3), dtype=object), TypeError),
        (np.zeros((1, 3), dtype=bool), TypeError),
        (np.zeros((1, 3), dtype=complex), TypeError),
        (np.ma.array(np.zeros((1, 3))), TypeError),
    ],
)
def test_invalid_points(bad, exc):
    with pytest.raises(exc):
        rdAlignment.GetAlignmentTransform(bad, [[0.0, 0.0, 0.0]])


@pytest.mark.parametrize(
    "weights,exc",
    [
        ([1.0], ValueError),
        ([-1.0] * 4, ValueError),
        ([0.0] * 4, ValueError),
        ([float("nan")] * 4, ValueError),
        ([float("inf")] * 4, ValueError),
        ([1e101] * 4, ValueError),
        ([1e100, 1e-300, 1.0, 1.0], ValueError),
        (np.zeros((4, 1)), ValueError),
        ([True] * 4, TypeError),
        ("1234", TypeError),
        (iter([1.0] * 4), TypeError),
        (np.ones(4, dtype=complex), TypeError),
        (np.ones(4, dtype=object), TypeError),
    ],
)
def test_invalid_weights(weights, exc):
    with pytest.raises(exc):
        rdAlignment.GetAlignmentTransform(REFERENCE, REFERENCE, weights)


def test_mismatched_counts_and_preflight_cap():
    with pytest.raises(ValueError, match="counts"):
        rdAlignment.GetAlignmentTransform(REFERENCE, REFERENCE[:3])
    too_many = np.broadcast_to(np.zeros(3), (65_537, 3))
    with pytest.raises(ValueError, match="65536"):
        rdAlignment.GetAlignmentTransform(too_many, too_many)


@pytest.mark.parametrize(
    "kwargs,exc",
    [
        ({"reflect": 1}, TypeError),
        ({"maxIterations": True}, TypeError),
        ({"maxIterations": 20}, ValueError),
        ({"maxIterations": 0}, ValueError),
    ],
)
def test_options_are_not_silently_ignored(kwargs, exc):
    with pytest.raises(exc):
        rdAlignment.GetAlignmentTransform(REFERENCE, REFERENCE, **kwargs)


def test_bounded_large_values_and_tiny_clouds_stay_finite():
    for scale in [1e-200, 1e90]:
        points = REFERENCE * scale
        ssd, matrix = rdAlignment.GetAlignmentTransform(points, points, [1e90] * 4)
        assert math.isfinite(ssd) and np.isfinite(matrix).all()
        np.testing.assert_allclose(
            transformed(points, matrix) / scale, REFERENCE, atol=1e-14
        )


def test_point_cap_boundary_is_supported():
    points = np.zeros((65_536, 3))
    ssd, matrix = rdAlignment.GetAlignmentTransform(points, points)
    assert ssd == 0.0
    np.testing.assert_array_equal(matrix, np.eye(4))
