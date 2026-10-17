"""Original reaction-role fixtures, bounded coordinate extraction and refusal."""

import itertools
import xml.etree.ElementTree as ET

import chematic
import pytest

ARROW = [[100, 0], [200, 0]]
COORDS = [
    [[0, 0], [30, 15], [60, 0]],
    [[90, 0]],
    [[150, -45]],
    [[240, 0], [270, 15], [300, 0]],
]


def test_roles_and_object_identity():
    text = chematic.write_reaction_cdxml_checked("CCO.O>Cl>CC=O", COORDS, ARROW)
    root = ET.fromstring(text)
    fragments = root.findall(".//fragment")
    step = root.find(".//step")
    assert step.attrib["ReactionStepReactants"].split() == [
        f.attrib["id"] for f in fragments[:2]
    ]
    assert step.attrib["ReactionStepObjectsAboveArrow"] == fragments[2].attrib["id"]
    assert step.attrib["ReactionStepProducts"] == fragments[3].attrib["id"]
    assert step.attrib["ReactionStepArrows"] == root.find(".//arrow").attrib["id"]
    ids = [node.attrib["id"] for node in root.iter() if "id" in node.attrib]
    assert len(ids) == len(set(ids))
    assert chematic.from_cdxml(text).smiles == chematic.from_smiles("CCO").smiles


def test_duplicate_components_and_no_agents():
    text = chematic.write_reaction_cdxml_checked("O.O>>O.O", [[[0, 0]]] * 4, ARROW)
    root = ET.fromstring(text)
    assert len(root.findall(".//fragment")) == 4
    assert root.find(".//step").attrib["ReactionStepObjectsAboveArrow"] == ""


@pytest.mark.parametrize("reaction", ["C", "C>C", ">>", "C>>", ">>C", "C>>>C"])
def test_invalid_reaction(reaction):
    with pytest.raises(ValueError):
        chematic.write_reaction_cdxml_checked(reaction, [], ARROW)


@pytest.mark.parametrize(
    "coords", [None, 1, "text", [[None]], [["text"]], [[[0, 0, 0]], [[1, 1]]]]
)
def test_bad_coordinate_types_or_shapes(coords):
    with pytest.raises((TypeError, ValueError)):
        chematic.write_reaction_cdxml_checked("C>>O", coords, ARROW)


@pytest.mark.parametrize(
    "coords",
    [[], [[[0, 0]]], [[[0, 0]]] * 3, [[], [[1, 1]]], [[[0, 0], [1, 1]], [[1, 1]]]],
)
def test_exact_coordinate_counts(coords):
    with pytest.raises(ValueError):
        chematic.write_reaction_cdxml_checked("C>>O", coords, ARROW)


@pytest.mark.parametrize(
    "value", [float("nan"), float("inf"), -float("inf"), 1_000_001, -1_000_001]
)
def test_nonfinite_and_large_values(value):
    with pytest.raises(ValueError):
        chematic.write_reaction_cdxml_checked("C>>O", [[[value, 0]], [[0, 0]]], ARROW)
    with pytest.raises(ValueError):
        chematic.write_reaction_cdxml_checked(
            "C>>O", [[[0, 0]], [[0, 0]]], [[0, 0], [value, 1]]
        )


@pytest.mark.parametrize("arrow", [None, [], [[0, 0]], [[0, 0]] * 3, [[0, 0], [0, 0]]])
def test_invalid_arrow(arrow):
    with pytest.raises((TypeError, ValueError)):
        chematic.write_reaction_cdxml_checked("C>>O", [[[0, 0]], [[0, 0]]], arrow)


@pytest.mark.parametrize(
    "smiles",
    [
        "c1ccccc1",
        "[C@H](F)(Cl)Br",
        "[CH3:1]",
        "*",
        "CC=CC",
        "C/C=C/C",
        "[O]",
        "[N]",
        "[CH3]",
        "[H]",
        "[SiH4]",
    ],
)
def test_unsupported_graphs(smiles):
    n = chematic.from_smiles(smiles).atom_count
    with pytest.raises(ValueError):
        chematic.write_reaction_cdxml_checked(
            smiles + ">>O", [[[0, 0]] * n, [[0, 0]]], ARROW
        )


def test_bounded_lazy_coordinates():
    seen = []

    def outer():
        for i in itertools.count():
            seen.append(i)
            yield [[0, 0]]

    with pytest.raises(ValueError, match="too many component"):
        chematic.write_reaction_cdxml_checked("C>>O", outer(), ARROW)
    assert len(seen) == 3
    seen.clear()

    def row():
        for i in itertools.count():
            seen.append(i)
            yield [0, 0]

    with pytest.raises(ValueError, match="too many atom"):
        chematic.write_reaction_cdxml_checked("C>>O", [row(), [[0, 0]]], ARROW)
    assert len(seen) == 2
    chematic.write_reaction_cdxml_checked(
        "C>>O", (iter([[0, 0]]) for _ in range(2)), ARROW
    )


def test_original_source_error():
    def broken():
        yield [[0, 0]]
        raise RuntimeError("coordinate source failed")

    with pytest.raises(RuntimeError, match="coordinate source failed"):
        chematic.write_reaction_cdxml_checked("C>>O", broken(), ARROW)


@pytest.mark.parametrize(
    "reaction",
    [
        "C" * 257 + ">>O",
        ".".join(["C"] * 33) + ">>O",
        "C" * 16385 + ">>O",
        ".".join(["C" * 256] * 4) + ">O>O",
    ],
)
def test_limits_before_coordinate_consumption(reaction):
    def source():
        raise AssertionError("must reject before reading coordinates")
        yield []

    with pytest.raises(ValueError):
        chematic.write_reaction_cdxml_checked(reaction, source(), ARROW)
