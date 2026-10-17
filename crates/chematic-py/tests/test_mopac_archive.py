"""Original fixtures for the bounded MOPAC archive reader."""

import pytest
import chematic

CART = "Summary\nFINAL GEOMETRY OBTAINED\nPM7 XYZ\nWater\n\nO 0 1 0 1 0 1\nH 1 1 0 1 0 1\nH 0 1 1 1 0 1\n"


def test_cartesian_archive():
    a = chematic.parse_mopac_archive(CART)
    assert a["source_text"] == CART
    assert a["title"] == "Water"
    assert a["summary_text"] == "Summary\n"
    assert a["atoms"][1] == {
        "element": "H",
        "label": None,
        "values": [1.0, 0.0, 0.0],
        "optimization_flags": ["1", "1", "1"],
        "references": None,
        "partial_charge": None,
    }


def test_internal_archive_with_symmetry():
    text = "Summary\nFINAL GEOMETRY OBTAINED CHARGE\nSYMMETRY\nOriginal H2\n\nH 0 0 0 0 0 0 0 0 0 0.1\nH .74 1 0 0 0 0 1 0 0 -0.1\n\n1 1 2\n"
    a = chematic.parse_mopac_archive(text)
    assert a["atoms"][1]["references"] == [1, 0, 0]
    assert a["atoms"][1]["values"] == [0.74, 0.0, 0.0]
    assert a["atoms"][1]["partial_charge"] == -0.1
    assert a["trailing_text"] == "\n1 1 2\n"


@pytest.mark.parametrize(
    "kwargs", [{"max_atoms": 2}, {"max_input_bytes": 10}, {"max_line_len": 4}]
)
def test_limits(kwargs):
    with pytest.raises(ValueError):
        chematic.parse_mopac_archive(CART, **kwargs)


@pytest.mark.parametrize(
    "text",
    [
        CART + CART,
        CART.replace("OBTAINED", "UNKNOWN"),
        CART.replace("H 1", "H NaN"),
        CART.replace("PM7 XYZ", "A0"),
    ],
)
def test_fail_closed(text):
    with pytest.raises(ValueError):
        chematic.parse_mopac_archive(text)
