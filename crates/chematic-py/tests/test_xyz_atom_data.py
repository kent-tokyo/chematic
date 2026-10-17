"""Geometry/payload preservation and strict input boundaries for xyz+data."""

from dataclasses import FrozenInstanceError
import math

import pytest
from chematic import XyzAtomDataFrame, parse_xyz_atom_data

VECTOR = "2\r\nopaque comment\r\nCl 1.25 -0 .5 +1e2 -.2 3\r\nC -2 4 0 5 6 7"


def test_vector_source_atom_order_and_frozen_snapshot():
    frame = parse_xyz_atom_data(VECTOR.encode())
    assert isinstance(frame, XyzAtomDataFrame)
    assert frame.elements == ("Cl", "C")
    assert frame.positions == ((1.25, -0.0, 0.5), (-2.0, 4.0, 0.0))
    assert math.copysign(1, frame.positions[0][1]) == -1
    assert frame.values == ((100.0, -0.2, 3.0), (5.0, 6.0, 7.0))
    assert frame.source_text.encode() == VECTOR.encode()
    assert repr(frame) == "XyzAtomDataFrame(atoms=2, components=3)"
    with pytest.raises(FrozenInstanceError):
        frame.comment = "new"
    with pytest.raises(TypeError):
        frame.values[0][0] = 0
    again = parse_xyz_atom_data(frame.to_xyz())
    assert (again.elements, again.positions, again.values, again.comment) == (
        frame.elements,
        frame.positions,
        frame.values,
        frame.comment,
    )


@pytest.mark.parametrize("components", [1, 3, 9])
def test_explicit_width_no_guessing_and_all_columns(components):
    data = tuple(float(i) / 7 for i in range(components))
    text = "1\n\nO 0 1 2 " + " ".join(repr(x) for x in data) + "\n"
    frame = XyzAtomDataFrame(text, components=components)
    assert frame.values == (data,) and frame.components == components
    assert frame.comment == ""
    assert f"payload:R:{components}" in frame.to_extxyz("payload")
    assert parse_xyz_atom_data(frame.to_xyz(), components=components).values == (data,)
    if components != 3:
        with pytest.raises(ValueError):
            parse_xyz_atom_data(text)


@pytest.mark.parametrize("value", [True, 3.0, "3", None])
def test_component_type(value):
    with pytest.raises(TypeError):
        parse_xyz_atom_data(VECTOR, components=value)


@pytest.mark.parametrize("value", [0, 2, 4, 8, 10])
def test_component_value(value):
    with pytest.raises(ValueError):
        parse_xyz_atom_data(VECTOR, components=value)


@pytest.mark.parametrize("text", [None, 1, bytearray(b"1"), memoryview(b"1")])
def test_serialized_types(text):
    with pytest.raises(TypeError):
        parse_xyz_atom_data(text)


@pytest.mark.parametrize(
    "text",
    [
        b"\xff",
        "\ud800",
        "1\rbad\rC 0 0 0 1 2 3",
        "1\nline\x00\nC 0 0 0 1 2 3",
        "1\nline\u2028\nC 0 0 0 1 2 3",
    ],
)
def test_encoding_and_controls(text):
    with pytest.raises(ValueError):
        parse_xyz_atom_data(text)


@pytest.mark.parametrize(
    "count", ["", "-1", "+1", "0", "50001", "1.0", "0000001", "１"]
)
def test_invalid_counts(count):
    with pytest.raises(ValueError, match="atom count"):
        parse_xyz_atom_data(count + "\ncomment\nC 0 0 0 1 2 3\n")


@pytest.mark.parametrize(
    "symbol", ["D", "T", "X", "*", "6", "CL", "cl", "C13", "[13C]"]
)
def test_unresolved_element_aliases_are_refused(symbol):
    with pytest.raises(ValueError, match="element"):
        parse_xyz_atom_data(f"1\ncomment\n{symbol} 0 0 0 1 2 3\n")


@pytest.mark.parametrize(
    "token",
    [
        "NaN",
        "inf",
        "-Infinity",
        "1e309",
        "1e-9999",
        "1_0",
        "1D2",
        "0x1p0",
        "１",
        "true",
    ],
)
def test_invalid_numeric_data(token):
    with pytest.raises(ValueError):
        parse_xyz_atom_data(f"1\n\nC 0 0 0 {token} 2 3\n")
    with pytest.raises(ValueError):
        parse_xyz_atom_data(f"1\n\nC {token} 0 0 1 2 3\n")


@pytest.mark.parametrize(
    "text",
    [
        "1\n",
        "1\ncomment\n",
        "1\ncomment\nC 0 0 0 1 2",
        "1\ncomment\nC 0 0 0 1 2 3 4",
        VECTOR + "\n1\nsecond\nH 0 0 0 1 2 3",
        "2\ncomment\nC 0 0 0 1 2 3\n",
    ],
)
def test_truncation_extra_columns_and_multiple_frames(text):
    with pytest.raises(ValueError):
        parse_xyz_atom_data(text)


@pytest.mark.parametrize(
    "comment",
    ["Properties=species:S:1:pos:R:3:data:R:3", 'Lattice="1 0 0 0 1 0 0 0 1"'],
)
def test_existing_extxyz_schema_must_use_existing_reader(comment):
    with pytest.raises(ValueError, match="from_extxyz"):
        parse_xyz_atom_data("1\n" + comment + "\nC 0 0 0 1 2 3\n")


@pytest.mark.parametrize(
    "name",
    [
        "",
        "pos",
        "species",
        "numbers",
        "positions",
        "Z",
        "a:b",
        "a b",
        "_private",
        "a" * 33,
        "a\nrow",
        "数値",
    ],
)
def test_property_name_rejections(name):
    with pytest.raises(ValueError):
        parse_xyz_atom_data(VECTOR).to_extxyz(name)


def test_minimum_subnormal_maximum_float_and_negative_zero():
    value = float.fromhex("0x0.0000000000001p-1022")
    huge = float.fromhex("0x1.fffffffffffffp+1023")
    frame = parse_xyz_atom_data(f"1\n\nOg -0 {huge!r} {value!r} 0.0 -0e-9999 {value!r}")
    again = parse_xyz_atom_data(frame.to_xyz())
    assert again.positions == frame.positions and again.values == frame.values
    assert math.copysign(1, again.values[0][1]) == -1


def test_size_limits_and_trailing_whitespace():
    with pytest.raises(ValueError, match="16 MiB"):
        parse_xyz_atom_data(b"x" * ((16 << 20) + 1))
    with pytest.raises(ValueError, match="4096 bytes"):
        parse_xyz_atom_data("1\n" + "漢" * 1366 + "\nC 0 0 0 1 2 3\n")
    with pytest.raises(ValueError, match="4096 bytes"):
        parse_xyz_atom_data("1\n\nC " + "0 " * 2100)
    text = "50000\nmaximum count\n" + "C 0 0 0 1 2 3\n" * 50000
    frame = parse_xyz_atom_data(text)
    assert len(frame.elements) == 50000
    assert len(parse_xyz_atom_data(frame.to_xyz()).values) == 50000
    assert parse_xyz_atom_data(VECTOR + "\n\t\n  ").source_text.endswith("\n\t\n  ")


def test_extxyz_comment_escaping_and_property_type():
    frame = parse_xyz_atom_data('1\n"quoted" ' + chr(92) + "path\nC 0 0 0 1 2 3\n")
    header = frame.to_extxyz().splitlines()[1]
    assert "_JSON" in header and "comment=" in header
    with pytest.raises(TypeError):
        frame.to_extxyz(3)
