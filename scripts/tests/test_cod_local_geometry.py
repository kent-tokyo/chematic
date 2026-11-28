import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "benchmarks"))
from cod_local_geometry import local_geometry


def test_local_geometry_ignores_neighbor_order():
    coordinates = [[0, 0, 0], [1, 0, 0], [0, 2, 0], [-3, 0, 0]]
    lengths, angles = local_geometry(coordinates, 0, [1, 2, 3])
    assert lengths == [1, 2, 3]
    assert angles == pytest.approx([90, 90, 180])
    assert local_geometry(coordinates, 0, [3, 1, 2]) == (lengths, angles)


def test_local_geometry_rejects_collapsed_or_nonfinite_bonds():
    with pytest.raises(ValueError):
        local_geometry([[0, 0, 0], [0, 0, 0]], 0, [1])
    with pytest.raises(ValueError):
        local_geometry([[0, 0, 0], [float("nan"), 0, 0]], 0, [1])
