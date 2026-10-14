"""Checked plain-CIF output: Open Babel #2818 and coordinate contracts."""
import math

import pytest
import chematic


CELL = dict(a=4.0, b=6.0, c=8.0, alpha=90.0, beta=90.0, gamma=90.0)


def silicon():
    return chematic.from_smiles("[Si]")


def test_singular_issue_cell_rejected():
    cell = dict(CELL, b=0.0)
    with pytest.raises(ValueError, match="length b"):
        chematic.write_cif_checked(silicon(), [[0, 0, 0]], cell)


@pytest.mark.parametrize("key", list(CELL))
@pytest.mark.parametrize("value", [0.0, -1.0, math.nan, math.inf])
def test_invalid_cell_fields_rejected(key, value):
    with pytest.raises(ValueError, match="invalid cell parameters"):
        chematic.write_cif_checked(silicon(), [[0, 0, 0]], dict(CELL, **{key: value}))


@pytest.mark.parametrize("angles", [(60, 60, 120), (10, 10, 150), (90, 90, 1e-6)])
def test_degenerate_angles_rejected(angles):
    cell = dict(CELL, **dict(zip(["alpha", "beta", "gamma"], angles)))
    with pytest.raises(ValueError, match="singular"):
        chematic.write_cif_checked(silicon(), [[0, 0, 0]], cell)


@pytest.mark.parametrize("cell", [None, CELL, dict(CELL, alpha=75, beta=83, gamma=64), dict(CELL, a=1e-6, b=2e-6, c=3e-6)])
def test_positions_round_trip(cell):
    coords = [[0.125123456789, -0.375987654321, 1.234567891234]]
    result = chematic.parse_cif(chematic.write_cif_checked(silicon(), coords, cell))
    assert result["cell"] == cell
    assert result["coords"][0] == pytest.approx(coords[0], abs=1e-12)


@pytest.mark.parametrize("coords", [[], [[0, 0, 0], [1, 2, 3]]])
def test_coordinate_count_rejected(coords):
    with pytest.raises(ValueError, match="requires 1 coordinates"):
        chematic.write_cif_checked(silicon(), coords)


@pytest.mark.parametrize("axis", [0, 1, 2])
@pytest.mark.parametrize("value", [math.nan, math.inf, -math.inf])
def test_nonfinite_coordinate_rejected(axis, value):
    coords = [[0.0, 0.0, 0.0]]
    coords[0][axis] = value
    with pytest.raises(ValueError, match="non-finite Cartesian"):
        chematic.write_cif_checked(silicon(), coords)


def test_missing_cell_key_rejected():
    with pytest.raises(ValueError, match="missing required key 'gamma'"):
        chematic.write_cif_checked(silicon(), [[0, 0, 0]], dict(a=4, b=6, c=8, alpha=90, beta=90))


def test_conversion_overflow_rejected():
    cell = dict(CELL, a=1e-100, b=1e-100, c=1e-100)
    with pytest.raises(ValueError, match="overflows"):
        chematic.write_cif_checked(silicon(), [[1e308, 0, 0]], cell)
