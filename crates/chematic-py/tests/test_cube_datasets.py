"""Synthetic CUBE dataset identity, voxel order and explicit scalar exports."""

from dataclasses import FrozenInstanceError
import itertools
import math

import numpy as np
import pytest

from chematic import (
    CubeDatasetDocument,
    CubeDatasetSelection,
    VolumetricGrid,
    read_cube_datasets,
)


def cube(natoms=-1, nval=None, ids=(7, 9), values="1 10 2 20", shape=(1, 1, 2)):
    header = f"{natoms} .125 -.25 .375"
    if nval is not None:
        header += f" {nval}"
    atoms = "6 4.0 .1 .2 .3\n" * abs(natoms)
    identifiers = f"{len(ids)} " + " ".join(map(str, ids)) + "\n" if natoms < 0 else ""
    return (
        "synthetic title\nscalar field\n"
        + header
        + "\n"
        + f"{shape[0]} .5 .1 0\n{shape[1]} 0 .75 .2\n{shape[2]} .3 0 1\n"
        + atoms
        + identifiers
        + values
        + "\n"
    )


def test_negative_atom_identifiers_are_demultiplexed_without_losing_source():
    source = cube()
    with pytest.raises(ValueError, match="multi-dataset"):
        VolumetricGrid.from_cube(source)
    doc = read_cube_datasets(source)
    assert isinstance(doc, CubeDatasetDocument)
    assert doc.source_text == source
    assert doc.shape == (1, 1, 2) and doc.dataset_count == 2
    assert doc.dataset_ids == (7, 9)
    assert doc.datasets == ((1.0, 2.0), (10.0, 20.0))
    selected = doc.select(1)
    assert isinstance(selected, CubeDatasetSelection)
    assert selected.source_text == source
    assert selected.dataset_index == 1 and selected.dataset_id == 9
    assert selected.omitted_indices == (0,)
    assert selected.text.splitlines()[2].split() == ["1", ".125", "-.25", ".375"]
    assert selected.text.splitlines()[:2] == source.splitlines()[:2]
    assert selected.text.splitlines()[3:7] == source.splitlines()[3:7]
    grid = selected.to_grid()
    assert grid.units == "bohr" and grid.shape == (1, 1, 2)
    assert grid.origin == (0.125, -0.25, 0.375)
    np.testing.assert_array_equal(grid.values, [10, 20])
    assert grid.atoms == [("C", 4.0, (0.1, 0.2, 0.3))]  # ECP nuclear charge


@pytest.mark.parametrize("negative,nval", [(False, 2), (True, None), (True, 1)])
def test_voxel_dataset_order_for_nonorthogonal_asymmetric_grid(negative, nval):
    # Deliberately asymmetric: dataset index fastest, then z, then y, then x.
    source_values = [
        1000 * x + 100 * y + 10 * z + d
        for x, y, z, d in itertools.product(range(2), range(3), range(4), range(2))
    ]
    text = cube(
        -1 if negative else 1,
        nval,
        (-17, 42),
        " ".join(map(str, source_values)),
        (2, 3, 4),
    )
    doc = read_cube_datasets(text)
    assert doc.dataset_ids == ((-17, 42) if negative else None)
    for dataset in range(2):
        grid = doc.select(dataset).to_grid()
        for x, y, z in itertools.product(range(2), range(3), range(4)):
            assert grid.get(x, y, z) == 1000 * x + 100 * y + 10 * z + dataset
        np.testing.assert_array_equal(
            grid.axes, [[0.5, 0.1, 0], [0, 0.75, 0.2], [0.3, 0, 1]]
        )


def test_wrapped_ids_crlf_and_values_in_same_line():
    text = cube(ids=(-3, 21)).replace("2 -3 21\n", "2\n -3\n21 ").replace("\n", "\r\n")
    doc = read_cube_datasets(text)
    assert doc.source_text == text
    assert doc.dataset_ids == (-3, 21)
    assert doc.datasets == ((1, 2), (10, 20))
    assert "\r" not in doc.select(0).text


@pytest.mark.parametrize("nval", [None, 1])
def test_ordinary_single_dataset_and_no_atom_case(nval):
    doc = read_cube_datasets(cube(0, nval, values="-.0 1.2345678901234567e-120"))
    assert doc.dataset_ids is None and doc.dataset_count == 1
    selected = doc.select(0)
    assert selected.omitted_indices == () and selected.dataset_id is None
    assert selected.to_grid().atoms == []
    assert math.copysign(1, selected.to_grid().values[0]) == -1
    assert selected.to_grid().values[1] == 1.2345678901234567e-120


def test_immutable_data_and_independent_native_arrays():
    doc = read_cube_datasets(cube())
    with pytest.raises(FrozenInstanceError):
        doc.shape = (9, 9, 9)
    with pytest.raises(TypeError):
        doc.datasets[0][0] = 99
    selected = doc.select(np.int64(1))
    with pytest.raises(FrozenInstanceError):
        selected.dataset_index = 0
    grid = selected.to_grid()
    grid.values[:] = 99
    assert doc.datasets[1] == (10, 20)
    np.testing.assert_array_equal(selected.to_grid().values, [10, 20])


@pytest.mark.parametrize(
    "index,error",
    [
        (-1, ValueError),
        (2, ValueError),
        (True, TypeError),
        (np.bool_(False), TypeError),
        (1.0, TypeError),
        ("1", TypeError),
    ],
)
def test_selection_requires_explicit_valid_positional_index(index, error):
    with pytest.raises(error):
        read_cube_datasets(cube()).select(index)


@pytest.mark.parametrize(
    "source",
    [
        "",
        "one comment\n",
        cube().replace("-1 .125", "-0 .125"),
        cube().replace("-1 .125", "-1_0 .125"),
        cube().replace("-1 .125", "-99999999999999999999999999 .125"),
        cube().replace("-1 .125 -.25 .375", "-1 0 0"),
        cube().replace("1 .5 .1 0", "-1 .5 .1 0"),
        cube().replace("1 0 .75 .2", "0 0 .75 .2"),
        cube().replace("2 .3 0 1", "2 .3 0"),
        cube().replace("2 7 9", "2 7 7"),
        cube().replace("2 7 9", "0"),
        cube().replace("2 7 9", "65"),
        cube().replace("2 7 9", "2 7 9223372036854775808"),
        cube(nval=2),
        cube(1, nval=0),
        cube(1, nval=65),
        cube(values="1 2 3"),
        cube(values="1 2 3 4 5"),
        cube(values="1 2 nan 4"),
        cube(values="1 2 1e999 4"),
        cube(values="1 2 1_0 4"),
        cube(values="1 2 1D2 4"),
        cube().replace(".125 -.25", "inf -.25"),
        cube().replace(".5 .1 0", "nan .1 0"),
        cube().replace("6 4.0", "6 1e999"),
        cube().replace("6 4.0", "0 0.0"),
        cube().replace(".1 .2 .3", ".1 .2"),
        cube().replace("synthetic", "synthe\x00tic"),
        cube().replace("synthetic", "s" * 16385),
        cube().replace("synthetic", "\ud800"),
    ],
)
def test_malformed_ambiguous_nonfinite_or_trailing_data_is_refused(source):
    with pytest.raises(ValueError):
        read_cube_datasets(source)


@pytest.mark.parametrize("value", [None, b"cube", 42, ["cube"]])
def test_input_must_be_text(value):
    with pytest.raises(TypeError):
        read_cube_datasets(value)


@pytest.mark.parametrize(
    "name,hard",
    [
        ("max_input_bytes", 8388608),
        ("max_atoms", 50000),
        ("max_grid_points", 262144),
        ("max_total_values", 262144),
    ],
)
def test_limits_may_only_be_lowered(name, hard):
    for value, error in [
        (0, ValueError),
        (hard + 1, ValueError),
        (True, TypeError),
        (1.0, TypeError),
    ]:
        with pytest.raises(error):
            read_cube_datasets(cube(), **{name: value})


def test_input_byte_limit_uses_utf8_and_exact_boundary():
    text = cube().replace("title", "日本語")
    size = len(text.encode())
    assert read_cube_datasets(text, max_input_bytes=size).dataset_count == 2
    with pytest.raises(ValueError, match="byte limit"):
        read_cube_datasets(text, max_input_bytes=size - 1)


def test_lowered_limits_and_declared_size_refusal_before_data():
    for text, kwargs in [
        (cube(-2), {"max_atoms": 1}),
        (cube(), {"max_grid_points": 1}),
        (cube(), {"max_total_values": 3}),
        (cube(shape=(2**62, 2**62, 2**62)), {}),
    ]:
        with pytest.raises(ValueError, match="exceeds limit"):
            read_cube_datasets(text, **kwargs)


def test_actual_dataset_and_total_scalar_boundaries():
    doc = read_cube_datasets(
        cube(
            ids=tuple(range(64)),
            shape=(1, 1, 1),
            values=" ".join(str(i) for i in range(64)),
        )
    )
    assert doc.dataset_count == 64 and doc.select(63).to_grid().get(0, 0, 0) == 63
    doc = read_cube_datasets(
        cube(
            natoms=0,
            shape=(64, 64, 64),
            values=("0 0 0 0 0 0\n" * (262144 // 6) + "0 0 0 0"),
        )
    )
    assert len(doc.datasets[0]) == 262144
    assert doc.select(0).to_grid().shape == (64, 64, 64)
    with pytest.raises(ValueError, match="total scalar"):
        read_cube_datasets(cube(natoms=1, nval=2, shape=(64, 64, 64), values=""))


def test_actual_atom_boundary_and_early_overflow():
    doc = read_cube_datasets(cube(natoms=50000, values="1 2"))
    assert len(doc.select(0).to_grid().atoms) == 50000
    with pytest.raises(ValueError, match="atom count"):
        read_cube_datasets(cube().replace("-1 .125", "50001 .125"))
