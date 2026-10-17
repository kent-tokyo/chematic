"""Explicit checked XYZ output profile; no bond/graph projection."""

import itertools
import math
import struct

import pytest
import chematic


def tokens(text):
    return [line.split() for line in text.splitlines()[2:]]


def test_requested_20_12f_format_and_original_order():
    coords = [
        [1.2345678901234, -2.3456789012345, -0.0],
        [0, 0, 0],
        [0.9572345678912, 0, 0],
    ]
    text = chematic.write_xyz_checked(["O", "H", "H"], coords, comment="water")
    assert text.splitlines()[:2] == ["3", "water"]
    assert text.splitlines()[2] == "O  " + " ".join(f"{v:20.12f}" for v in coords[0])
    rows = tokens(text)
    assert [r[0] for r in rows] == ["O", "H", "H"]
    assert rows[0][1:] == ["1.234567890123", "-2.345678901234", "-0.000000000000"]
    assert math.copysign(1, float(rows[0][3])) == -1
    assert coords[0][0] == 1.2345678901234


@pytest.mark.parametrize("notation", ["fixed", "scientific"])
@pytest.mark.parametrize("places", [0, 1, 6, 12, 16, 17])
def test_decimal_places_against_python_numeric_formatter(notation, places):
    coords = [[1.25, -2.75, -0.0]]
    row = tokens(
        chematic.write_xyz_checked(
            ["Ne"], coords, notation=notation, decimal_places=places, min_width=0
        )
    )[0]
    for token, value in zip(row[1:], coords[0]):
        expected = format(value, f'.{places}{"f" if notation=="fixed" else "e"}')
        if notation == "fixed":
            assert token == expected
        else:
            a, b = token.split("e")
            c, d = expected.split("e")
            assert a == c and int(b) == int(d)


def test_scientific_preserves_finite_extremes_and_signed_zero():
    values = [
        0.0,
        -0.0,
        5e-324,
        2.2250738585072014e-308,
        1.7976931348623157e308,
        -1.7976931348623157e308,
        1 / 3,
    ]
    for value in values:
        text = chematic.write_xyz_checked(
            ["C"], [[value] * 3], notation="scientific", decimal_places=16, min_width=0
        )
        for token in tokens(text)[0][1:]:
            assert struct.pack(">d", float(token)) == struct.pack(">d", value)
    for value in values[-3:-1]:
        with pytest.raises(ValueError, match="rounded coordinate"):
            chematic.write_xyz_checked(
                ["C"], [[value, 0, 0]], notation="scientific", decimal_places=0
            )


def test_large_fixed_fields_expand_without_fusing_tokens():
    text = chematic.write_xyz_checked(
        ["Og"], [[1e100, -1e100, 0]], min_width=1, decimal_places=17
    )
    assert len(tokens(text)[0]) == 4
    assert float(tokens(text)[0][1]) == 1e100


@pytest.mark.parametrize(
    "kwargs",
    [
        {"decimal_places": 18},
        {"decimal_places": -1},
        {"min_width": 65},
        {"min_width": -1},
        {"max_output_bytes": 0},
        {"max_output_bytes": 64 * 1024 * 1024 + 1},
        {"notation": "auto"},
        {"notation": "FIXED"},
        {"comment": "bad\nnew frame"},
        {"comment": "Properties=species:S:1:pos:R:3"},
        {"comment": "Lattice=cell"},
        {"comment": "café"},
        {"comment": "x" * 4097},
        {"comment": "\t"},
    ],
)
def test_invalid_options(kwargs):
    with pytest.raises(ValueError):
        chematic.write_xyz_checked(["O"], [[0, 0, 0]], **kwargs)


@pytest.mark.parametrize("field", ["decimal_places", "min_width", "max_output_bytes"])
@pytest.mark.parametrize("value", [True, 1.5, "2"])
def test_strict_integer_options(field, value):
    with pytest.raises((ValueError, TypeError)):
        chematic.write_xyz_checked(["O"], [[0, 0, 0]], **{field: value})


@pytest.mark.parametrize(
    "elements,positions",
    [
        (["O"], []),
        ([], [[0, 0, 0]]),
        (["O"], [[0, 0]]),
        (["O"], [[0, 0, 0, 0]]),
        (["D"], [[0, 0, 0]]),
        (["cl"], [[0, 0, 0]]),
        (["*"], [[0, 0, 0]]),
        ("O", [[0, 0, 0]]),
        (["O"], "123"),
        (["O"], ["123"]),
        (["O"], [[True, 0, 0]]),
    ],
)
def test_invalid_geometry(elements, positions):
    with pytest.raises((ValueError, TypeError)):
        chematic.write_xyz_checked(elements, positions)


@pytest.mark.parametrize("value", [math.nan, math.inf, -math.inf])
def test_nonfinite_coordinates(value):
    with pytest.raises(ValueError, match="finite"):
        chematic.write_xyz_checked(["O"], [[0, value, 0]])


def test_finite_output_budget_exact_endpoint_and_empty_geometry():
    text = chematic.write_xyz_checked(["O"], [[0, 0, 0]])
    assert (
        chematic.write_xyz_checked(["O"], [[0, 0, 0]], max_output_bytes=len(text))
        == text
    )
    with pytest.raises(ValueError, match="byte limit"):
        chematic.write_xyz_checked(["O"], [[0, 0, 0]], max_output_bytes=len(text) - 1)
    assert chematic.write_xyz_checked([], []) == "0\n\n"


def test_bounded_lazy_geometry_and_original_iterator_errors():
    text = chematic.write_xyz_checked(
        iter(["O", "H"]), (iter(row) for row in [[0, 0, 0], [1, 0, 0]])
    )
    assert len(tokens(text)) == 2
    with pytest.raises(ValueError, match="100000"):
        chematic.write_xyz_checked(itertools.repeat("H"), itertools.repeat([0, 0, 0]))

    def broken():
        yield 0.0
        yield 0.0
        yield 0.0
        raise RuntimeError("row failed")

    with pytest.raises(RuntimeError, match="row failed"):
        chematic.write_xyz_checked(["H"], [broken()])

    def trailing():
        yield [0, 0, 0]
        raise RuntimeError("outer failed")

    with pytest.raises(RuntimeError, match="outer failed"):
        chematic.write_xyz_checked(["H"], trailing())


def test_old_molecule_writer_default_is_unchanged():
    mol = chematic.from_smiles("O")
    assert "1.234568" in mol.to_xyz([[1.2345678901234, 0, 0]], "legacy")
