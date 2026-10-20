"""Original static-scene input controls for the XBS dialect in OB #2576."""

from dataclasses import FrozenInstanceError
import math
import itertools
import struct

import pytest

from chematic import XbsScene, XbsSpecies, XbsBondRule, read_xbs, write_xbs


SCENE = (
    "* Original abstract two-site scene\r\n"
    "\t bonds Site1 H 0 2.5 0.125 .75\r\n"
    "atom Site1 +1.25 -0.0 4.9406564584124654e-324\r\n"
    "\r\n"
    "spec H .5 1\r\n"
    " atom H -2.5 0 0\r\n"
    "spec Site1 1.0 .25\r\n"
)


def test_source_order_lineage_forward_references_and_deep_immutability():
    result = read_xbs(SCENE.encode("ascii"))
    assert isinstance(result, XbsScene)
    assert result.source_text == SCENE
    assert result.to_xbs() == SCENE
    assert result.to_xbs().encode("ascii") == SCENE.encode("ascii")
    assert result.labels == ("Site1", "H")
    assert result.atom_lines == (3, 6)
    assert result.species_lines == (5, 7)
    assert result.bond_lines == (2,)
    assert tuple(s.label for s in result.species) == ("H", "Site1")
    assert result.bonds == (XbsBondRule("Site1", "H", 0, 2.5, 0.125, 0.75),)
    assert result.positions[0] == (1.25, -0.0, float.fromhex("0x0.0000000000001p-1022"))
    assert math.copysign(1, result.positions[0][1]) == -1
    with pytest.raises(FrozenInstanceError):
        result.labels = ("C",)
    with pytest.raises(TypeError):
        result.positions[0][0] = 2
    with pytest.raises(FrozenInstanceError):
        result.species[0].radius = 2
    assert len(repr(result)) < 100


def test_normalized_output_round_trips_binary64_and_scene_fields():
    result = read_xbs(SCENE)
    normalized = result.to_xbs(preserve_source=False)
    assert normalized == write_xbs(
        result.labels, result.positions, species=result.species, bonds=result.bonds
    )
    assert "*" not in normalized and "\r" not in normalized
    again = XbsScene(normalized)
    assert again.labels == result.labels
    assert again.species == result.species
    assert again.bonds == result.bonds
    assert b"".join(
        struct.pack("!d", v) for row in again.positions for v in row
    ) == b"".join(struct.pack("!d", v) for row in result.positions for v in row)


def test_abstract_labels_extra_unused_styles_and_no_chemical_inference():
    r = read_xbs("spec Unused 1 0\natom Site_2 0 0 0\nspec Site_2 1 0.5")
    assert r.labels == ("Site_2",)
    assert len(r.species) == 2
    assert r.bonds == ()
    assert not hasattr(r, "molecule")


@pytest.mark.parametrize(
    "source", [None, 12, bytearray(b"atom C 0 0 0"), memoryview(b"")]
)
def test_wrong_input_types(source):
    with pytest.raises(TypeError):
        read_xbs(source)


@pytest.mark.parametrize(
    "source", ["spec C 1 0\natom C 0 0 \u2212\u00bd", b"\xff", "\ud800"]
)
def test_only_ascii_is_accepted(source):
    with pytest.raises(ValueError, match="ASCII"):
        read_xbs(source)


@pytest.mark.parametrize(
    "record",
    [
        "frame new",
        "scale 2",
        "save output.bs",
        "cell 1 1 1",
        "# comment",
        "ATOM C 0 0 0",
        "2",
        "C 0 0 0 2",
        "spec C 1 red",
        "spec C 1 0 0 0",
        "bonds * C 0 2 0.1 0",
        "atom C+ 0 0 0",
        "atom Carbon123 0 0 0",
        "atom C 0 0",
        "atom C 0 0 0 extra",
        "spec C 1",
        "bonds C C 0 1 .1",
        "atom C 0 0 0 * inline comment",
        "atom C 0 0 0; spec C 1 0",
    ],
)
def test_unknown_and_unsupported_records_are_errors_at_the_physical_line(record):
    with pytest.raises(ValueError, match="XBS line 3:"):
        read_xbs("* heading\n\n" + record + "\n")


@pytest.mark.parametrize(
    "value",
    ["NaN", "inf", "1e9999", "1e-9999", "1_000", "0x1p0", "1D0", "true", "1000001"],
)
def test_invalid_numeric_literals_fail_without_partial_scene(value):
    with pytest.raises(ValueError, match="XBS line 2:"):
        read_xbs(f"spec C 1 0\natom C {value} 0 0\n")


@pytest.mark.parametrize(
    "record",
    [
        "spec C 0 0",
        "spec C -1 0",
        "spec C 1 1.1",
        "spec C 1 -1",
        "bonds C C 2 1 0.1 0",
        "bonds C C 0 2 0 0",
        "bonds C C -1 2 .1 0",
        "bonds C C 0 2 0.1 2",
    ],
)
def test_writer_style_and_rule_bounds_are_enforced(record):
    with pytest.raises(ValueError, match="XBS line 2:"):
        read_xbs("atom C 0 0 0\n" + record)


@pytest.mark.parametrize(
    "source,line,message",
    [
        ("atom C 0 0 0", 1, "missing style"),
        ("spec C 1 0\natom C 0 0 0\nbonds C H 0 2 .1 0", 3, "undefined species"),
        ("spec C 1 0\nspec C 2 1\natom C 0 0 0", 2, "duplicate species"),
        (
            "atom C 0 0 0\nspec C 1 0\nspec H 1 1\nbonds C H 0 2 .1 0\nbonds H C 0 3 .2 1",
            5,
            "duplicate undirected",
        ),
    ],
)
def test_reference_and_duplicate_errors_retain_source_lines(source, line, message):
    with pytest.raises(ValueError, match=f"XBS line {line}:.*{message}"):
        read_xbs(source)


@pytest.mark.parametrize("source", ["", " \n\t\n* comment only\n", "spec C 1 0\n"])
def test_empty_scene_is_refused(source):
    with pytest.raises(ValueError, match="at least one atom"):
        read_xbs(source)


@pytest.mark.parametrize("control", ["\x00", "\x0b", "\x0c", "\x1c", "\x7f", "\r"])
def test_controls_and_bare_cr_are_refused_even_in_comments(control):
    with pytest.raises(ValueError, match="XBS line 1:"):
        read_xbs("* comment" + control + "ending\natom C 0 0 0\nspec C 1 0\n")


def test_byte_line_and_record_limits():
    with pytest.raises(ValueError, match="4 MiB"):
        read_xbs(b"*" * ((4 << 20) + 1))
    with pytest.raises(ValueError, match="20000"):
        read_xbs("\n" * 20001)
    with pytest.raises(ValueError, match="4096"):
        read_xbs("*" * 4097)
    with pytest.raises(ValueError, match="atom count"):
        read_xbs("spec C 1 0\n" + "atom C 0 0 0\n" * 10001)
    with pytest.raises(ValueError, match="species count"):
        read_xbs("atom C 0 0 0\n" + "\n".join(f"spec S{i} 1 0" for i in range(129)))
    rules = [
        f"bonds S{i} S{j} 0 2 .1 0"
        for i, j in itertools.islice(
            itertools.combinations_with_replacement(range(46), 2), 1025
        )
    ]
    with pytest.raises(ValueError, match="bond-rule count"):
        read_xbs("atom C 0 0 0\n" + "\n".join(rules))
    styles = "\n".join(f"spec S{i} 1 0" for i in range(128))
    accepted = read_xbs("atom S0 0 0 0\n" + styles + "\n" + "\n".join(rules[:1024]))
    assert len(accepted.species) == 128
    assert len(accepted.bonds) == 1024


def test_supported_boundary_scene_and_line_counts():
    r = read_xbs("spec C 1 0\n" + "atom C 1000000 -1000000 -0\n" * 10000)
    assert len(r.positions) == 10000
    assert r.atom_lines[-1] == 10001
    padded = "*" * 4096 + "\n" + "\n" * 19997 + "atom C 0 0 0\nspec C 1 0"
    assert read_xbs(padded).atom_lines == (19999,)
    base = "atom C 0 0 0\nspec C 1 0\n"
    blocks, remainder = divmod((4 << 20) - len(base), 4096)
    tail = "*" * (remainder - 1) + "\n" if remainder else ""
    maximum = ("*" * 4095 + "\n") * blocks + tail + base
    assert len(maximum.encode("ascii")) == 4 << 20
    assert read_xbs(maximum).to_xbs() == maximum


def test_constructor_also_validates_and_source_policy_requires_bool():
    with pytest.raises(ValueError):
        XbsScene("frame bad")
    for invalid in (0, "yes", None):
        with pytest.raises(TypeError):
            read_xbs(SCENE).to_xbs(preserve_source=invalid)
