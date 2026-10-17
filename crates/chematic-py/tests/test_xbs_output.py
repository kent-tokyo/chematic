"""XBS output contracts: scene data, explicit styles and bounded inputs."""

from collections.abc import Sequence
from dataclasses import FrozenInstanceError
import math

import pytest
from chematic import XbsSpecies, XbsBondRule, write_xbs


def test_public_scene_preserves_atom_order_numbers_and_styles():
    xyz = [[1.2345678901234567, -0.0, 1e-200], [-2.5, 1e6, -1e6]]
    labels = ["Cl1", "C2"]
    styles = [XbsSpecies("C2", 0.7, 0), XbsSpecies("Cl1", 1.3, 1)]
    rules = [XbsBondRule("Cl1", "C2", 0.3, 1.7, 0.12, 0.5)]
    text = write_xbs(labels, xyz, species=styles, bonds=rules)
    rows = [line.split() for line in text.splitlines()]
    assert [row[:2] for row in rows[:2]] == [["atom", "Cl1"], ["atom", "C2"]]
    assert [[float(x) for x in row[2:]] for row in rows[:2]] == xyz
    assert math.copysign(1, float(rows[0][3])) == -1
    assert rows[2] == ["spec", "C2", "0.69999999999999996", "0"]
    assert [float(x) for x in rows[-1][3:]] == [0.3, 1.7, 0.12, 0.5]
    assert text.endswith("\n")
    assert text == write_xbs(labels, xyz, species=styles, bonds=rules)
    assert labels == ["Cl1", "C2"] and xyz[0][0] == 1.2345678901234567


def test_no_implicit_bonds_or_styles_and_snapshot_output():
    xyz = [[0, 0, 0]]
    text = write_xbs(("Site1",), xyz, species=(XbsSpecies("Site1", 1, 0.5),))
    xyz[0][0] = 9
    assert text == "atom Site1 0 0 0\nspec Site1 1 0.5\n"
    with pytest.raises(FrozenInstanceError):
        XbsSpecies("C", 1, 0.5).radius = 2


@pytest.mark.parametrize(
    "label",
    ["", "C C", "C\natom", "*", "C*", "C+", "1C", "Carbon123", "Ｃ", "C-a", "C\x00"],
)
def test_invalid_labels(label):
    with pytest.raises(ValueError):
        XbsSpecies(label, 1, 0.5)


@pytest.mark.parametrize("value", [True, "1", 1j, None])
def test_numeric_type_errors(value):
    with pytest.raises(TypeError):
        XbsSpecies("C", value, 0.5)
    with pytest.raises(TypeError):
        XbsBondRule("C", "C", 0, 2, 0.1, value)
    with pytest.raises(TypeError):
        write_xbs(["C"], [[value, 0, 0]], species=[XbsSpecies("C", 1, 0.5)])


@pytest.mark.parametrize("radius", [0, -1, math.nan, math.inf, 1e7, 10**1000])
def test_radius_ranges(radius):
    with pytest.raises(ValueError):
        XbsSpecies("C", radius, 0.5)
    with pytest.raises(ValueError):
        XbsBondRule("C", "C", 0, 2, radius, 0.5)


@pytest.mark.parametrize("gray", [-0.1, 1.1, math.nan, math.inf])
def test_gray_ranges(gray):
    with pytest.raises(ValueError):
        XbsSpecies("C", 1, gray)


@pytest.mark.parametrize(
    "lower,upper",
    [(1, 1), (2, 1), (-1, 2), (0, 0), (0, 1e7), (math.nan, 2), (0, math.inf)],
)
def test_rule_ranges(lower, upper):
    with pytest.raises(ValueError):
        XbsBondRule("C", "C", lower, upper, 0.1, 0.5)


@pytest.mark.parametrize(
    "labels,xyz,error",
    [
        ([], [], ValueError),
        (["C"], [], ValueError),
        (["C"], [[0, 0]], ValueError),
        (["C"], [[0, 0, 0, 0]], ValueError),
        (["C"], [[math.nan, 0, 0]], ValueError),
        (["C"], [[1e7, 0, 0]], ValueError),
        ("C", [[0, 0, 0]], TypeError),
        (["C"], ["000"], TypeError),
        ([b"C"], [[0, 0, 0]], TypeError),
        ((x for x in ["C"]), [[0, 0, 0]], TypeError),
    ],
)
def test_atom_validation(labels, xyz, error):
    with pytest.raises(error):
        write_xbs(labels, xyz, species=[XbsSpecies("C", 1, 0.5)])


def test_style_and_rule_references():
    a = XbsSpecies("C", 1, 0.5)
    b = XbsSpecies("H", 0.5, 1)

    def write(species, bonds=()):
        return write_xbs(["C"], [[0, 0, 0]], species=species, bonds=bonds)

    for styles in ([], [a, a], [b]):
        with pytest.raises(ValueError):
            write(styles)
    with pytest.raises(TypeError):
        write([("C", 1, 0.5)])
    with pytest.raises(TypeError):
        write([a], [("C", "C", 0, 2, 0.1, 0.5)])
    with pytest.raises(ValueError):
        write([a], [XbsBondRule("C", "H", 0, 2, 0.1, 0.5)])
    with pytest.raises(ValueError):
        write(
            [a, b],
            [
                XbsBondRule("C", "H", 0, 2, 0.1, 0.5),
                XbsBondRule("H", "C", 0, 3, 0.2, 0.6),
            ],
        )


def test_caps_before_copy_or_indexing():
    class Huge(Sequence):
        def __init__(self, n):
            self.n = n

        def __len__(self):
            return self.n

        def __getitem__(self, i):
            raise AssertionError("must not index over-cap input")

    a = XbsSpecies("C", 1, 0.5)
    for kwargs in (
        dict(labels=Huge(10001), positions=[]),
        dict(labels=["C"], positions=Huge(10001)),
        dict(labels=["C"], positions=[[0, 0, 0]], species=Huge(129)),
        dict(labels=["C"], positions=[[0, 0, 0]], bonds=Huge(1025)),
    ):
        with pytest.raises(ValueError):
            write_xbs(**({"species": [a]} | kwargs))
    result = write_xbs(["C"] * 10000, [[0, 0, 0]] * 10000, species=[a])
    assert result.count("atom C ") == 10000


def test_indexed_sequence_does_not_consume_unbounded_iterator():
    class IndexedList(list):
        def __iter__(self):
            raise AssertionError("only bounded indexed access is allowed")

    assert "spec C 1 0.5" in write_xbs(
        ["C"], [[0, 0, 0]], species=IndexedList([XbsSpecies("C", 1, 0.5)])
    )
