"""Independent analytical and typed input contracts for explicit Coulomb fields."""

from dataclasses import FrozenInstanceError
import itertools
import math

import numpy as np
import pytest
from chematic.rdkit_compat import rdMIF

K = 1389.048632870505


def test_single_charge_and_squared_strict_threshold():
    field = rdMIF.Coulomb([1.0], [[0.0, 0.0, 0.0]])
    assert field.site_count == 1 and repr(field) == "Coulomb(sites=1)"
    assert field(1.0, 0.0, 0.0, 2.0) == K
    assert field(2.0, 0.0, 0.0, 4.0) == 0.0
    assert field(2.0, 0.0, 0.0, math.nextafter(4.0, math.inf)) == K / 2
    assert field(0.0, 0.0, 0.0, 0.0) == 0.0
    assert field(0.0, 0.0, 0.0, 1.0) == K


def test_probe_sign_dipole_and_negative_absolute_total():
    q = [1.0, -1.0]
    xyz = [[-1.0, 0.0, 0.0], [1.0, 0.0, 0.0]]
    assert rdMIF.Coulomb(q, xyz)(0.0, 0.0, 0.0, 2.0) == 0.0
    value = rdMIF.Coulomb(q, xyz)(-2.0, 0.0, 0.0, 10.0)
    assert value == pytest.approx(K * (1 - 1 / 3))
    assert rdMIF.Coulomb(q, xyz, probeCharge=-2.0)(
        -2.0, 0.0, 0.0, 10.0
    ) == pytest.approx(-2 * value)
    assert rdMIF.Coulomb(q, xyz, absVal=True)(-2.0, 0.0, 0.0, 10.0) == -abs(value)
    assert (
        rdMIF.Coulomb([1.0], [[0.0, 0.0, 0.0]], probeCharge=0.0)(1.0, 0.0, 0.0, 2.0)
        == 0.0
    )


@pytest.mark.parametrize("alpha", [0.0, 1e-6, 0.1, 2.0, 1e200])
def test_softcore_and_hard_floor(alpha):
    field = rdMIF.Coulomb([2.0], [[0.0, 0.0, 0.0]], softcoreParam=alpha, cutoff=0.001)
    expected = 2 * K / math.sqrt(alpha if alpha else 1e-6)
    assert field(0.0, 0.0, 0.0, 1.0) == pytest.approx(expected, rel=1e-15)
    assert field(0.1, 0.0, 0.0, 1.0) == pytest.approx(
        2 * K / math.sqrt(alpha + 0.1 * 0.1 if alpha else 0.1 * 0.1), rel=1e-15
    )


@pytest.mark.parametrize("array", [False, True])
def test_inputs_are_copied_and_field_is_frozen(array):
    q = np.array([1.0, 2.0]) if array else [1.0, 2.0]
    xyz = (
        np.array([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]])
        if array
        else [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]]
    )
    field = rdMIF.Coulomb(q, xyz)
    before = field(2.0, 0.0, 0.0, 10.0)
    q[0] = 100.0
    xyz[0][0] = 100.0
    assert field(2.0, 0.0, 0.0, 10.0) == before
    with pytest.raises(FrozenInstanceError):
        field._alpha = 1.0


def test_integer_sequences_noncontiguous_numpy_and_numpy_scalars():
    q = np.array([1.0, 99.0, -1.0, 99.0])[::2]
    xyz = np.array([[-1.0, 0.0, 0.0], [1.0, 0.0, 0.0]])[:, ::-1]
    field = rdMIF.Coulomb(q, xyz, probeCharge=np.float64(2.0), cutoff=np.float32(1.0))
    assert field(0, 0, 0, np.float64(2.0)) == 0.0
    assert rdMIF.Coulomb((1,), ((0, 0, 0),))(1, 0, 0, 2) == K


@pytest.mark.parametrize(
    "charges,positions",
    [
        (1.0, [[0, 0, 0]]),
        ("1", [[0, 0, 0]]),
        ([1.0], "000"),
        (iter([1.0]), [[0, 0, 0]]),
        ([1.0], iter([[0, 0, 0]])),
        ({0: 1.0}, [[0, 0, 0]]),
        ([True], [[0, 0, 0]]),
        ([np.bool_(True)], [[0, 0, 0]]),
        ([1j], [[0, 0, 0]]),
        (["1"], [[0, 0, 0]]),
        ([1.0], [["0", 0, 0]]),
        ([1.0], [[False, 0, 0]]),
        ([1.0], [[0j, 0, 0]]),
        (np.ones((1, 1)), [[0, 0, 0]]),
        ([1.0], [iter([0, 0, 0])]),
    ],
)
def test_type_errors(charges, positions):
    with pytest.raises(TypeError):
        rdMIF.Coulomb(charges, positions)


@pytest.mark.parametrize(
    "charges,positions",
    [
        ([1.0], []),
        ([1.0], [[0, 0]]),
        ([1.0], [[0, 0, 0, 0]]),
        ([math.nan], [[0, 0, 0]]),
        ([math.inf], [[0, 0, 0]]),
        ([1e6 + 1], [[0, 0, 0]]),
        ([1.0], [[math.inf, 0, 0]]),
        ([1.0], [[math.nan, 0, 0]]),
        ([1.0], [[1e101, 0, 0]]),
    ],
)
def test_value_errors(charges, positions):
    with pytest.raises(ValueError):
        rdMIF.Coulomb(charges, positions)


@pytest.mark.parametrize(
    "key,value",
    [
        ("probeCharge", True),
        ("probeCharge", "1"),
        ("probeCharge", 1j),
        ("absVal", 1),
        ("softcoreParam", True),
        ("cutoff", False),
    ],
)
def test_parameter_types(key, value):
    with pytest.raises(TypeError):
        rdMIF.Coulomb([1.0], [[0, 0, 0]], **{key: value})


@pytest.mark.parametrize(
    "key,value",
    [
        ("probeCharge", 1e7),
        ("probeCharge", math.nan),
        ("softcoreParam", -1.0),
        ("softcoreParam", 1e-7),
        ("softcoreParam", math.inf),
        ("softcoreParam", 1e201),
        ("cutoff", 0.0),
        ("cutoff", 0.0001),
        ("cutoff", math.inf),
        ("cutoff", 1e101),
    ],
)
def test_parameter_ranges(key, value):
    with pytest.raises(ValueError):
        rdMIF.Coulomb([1.0], [[0, 0, 0]], **{key: value})


@pytest.mark.parametrize(
    "args",
    [
        (True, 0.0, 0.0, 1.0),
        ("0", 0.0, 0.0, 1.0),
        (0.0, 0.0, 0.0, "1"),
        (0.0, 0.0, 0.0, True),
    ],
)
def test_call_type_errors(args):
    with pytest.raises(TypeError):
        rdMIF.Coulomb([1.0], [[0, 0, 0]])(*args)


@pytest.mark.parametrize(
    "args",
    [
        (math.nan, 0.0, 0.0, 1.0),
        (0.0, math.inf, 0.0, 1.0),
        (0.0, 0.0, 1e101, 1.0),
        (0.0, 0.0, 0.0, -1.0),
        (0.0, 0.0, 0.0, math.inf),
        (0.0, 0.0, 0.0, 1e201),
    ],
)
def test_call_value_errors(args):
    with pytest.raises(ValueError):
        rdMIF.Coulomb([1.0], [[0, 0, 0]])(*args)


def test_empty_field_and_large_finite_ranges():
    assert rdMIF.Coulomb([], [])(0.0, 0.0, 0.0, 1.0) == 0.0
    assert (
        math.copysign(1.0, rdMIF.Coulomb([], [], absVal=True)(0.0, 0.0, 0.0, 1.0)) == -1
    )
    field = rdMIF.Coulomb([1e6], [[-1e100, 0.0, 0.0]], probeCharge=1e6)
    assert field(1e100, 0.0, 0.0, 1e200) == 0.0
    assert math.isfinite(field(-1e100, 0.0, 0.0, 1e200))


def test_site_cap_before_indexing_and_maximum_supported_count():
    class Huge:
        def __len__(self):
            return 100001

        def __getitem__(self, index):
            raise AssertionError("must not read huge input")

    with pytest.raises(ValueError, match="100000"):
        rdMIF.Coulomb(Huge(), Huge())
    assert (
        rdMIF.Coulomb([0.0] * 100000, [[0.0, 0.0, 0.0]] * 100000).site_count == 100000
    )
    with pytest.raises(TypeError):
        rdMIF.Coulomb(itertools.repeat(0.0), itertools.repeat([0.0, 0.0, 0.0]))


def test_indexing_errors_propagate_without_changing_inputs():
    class Broken:
        def __len__(self):
            return 2

        def __getitem__(self, i):
            if i:
                raise RuntimeError("source failed")
            return 1.0

    xyz = [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]]
    with pytest.raises(RuntimeError, match="source failed"):
        rdMIF.Coulomb(Broken(), xyz)
    assert xyz == [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]]
