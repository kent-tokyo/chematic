"""Periodic conversion preserves basis, geometry, provenance and constraints."""

from dataclasses import FrozenInstanceError

import numpy as np
import pytest
import chematic as ch


TMOL = """# original comment
$periodic 3
$lattice angs
 4 0.5 0
 0 5 0.75
 0.25 0 6
$coord frac
 -0.25 0.5 1.25 o f
 0.25 0 0.125 si
 0.5 0.25 -0.5 o
$end
"""


def test_bidirectional_conversion_retains_geometry_and_atom_mapping():
    original = ch.parse_turbomole(TMOL)
    result = ch.turbomole_to_poscar_checked(TMOL)
    assert result.source_text == TMOL
    assert (result.source_format, result.target_format) == ("turbomole", "poscar")
    assert result.source_indices == (0, 2, 1)
    assert result.text.splitlines()[5:9] == [
        "O Si",
        "2 1",
        "Selective dynamics",
        "Direct",
    ]
    structure = ch.PeriodicStructure.from_poscar(result.text)
    np.testing.assert_allclose(
        structure.lattice.matrix, original["lattice_angstrom"], atol=1e-14
    )
    np.testing.assert_allclose(
        structure.cartesian_positions(),
        np.array(original["positions_angstrom"])[[0, 2, 1]],
        atol=1e-13,
    )
    assert structure.fractional_positions()[0, 0] < 0
    assert structure.fractional_positions()[0, 2] > 1
    back = ch.poscar_to_turbomole_checked(result.text)
    assert back.source_text == result.text
    assert back.source_indices == (0, 1, 2)
    assert (back.source_format, back.target_format) == ("poscar", "turbomole")
    restored = ch.parse_turbomole(back.text)
    assert restored["elements"] == ["O", "O", "Si"]
    assert restored["frozen"] == [True, False, False]
    np.testing.assert_allclose(
        restored["lattice_angstrom"], original["lattice_angstrom"], atol=1e-14
    )
    np.testing.assert_allclose(
        restored["positions_angstrom"],
        np.array(original["positions_angstrom"])[[0, 2, 1]],
        atol=1e-13,
    )
    with pytest.raises(FrozenInstanceError):
        result.text = "changed"


@pytest.mark.parametrize("unit", ["", " bohr", " angs"])
def test_native_turbomole_units_and_cell(unit):
    text = (
        f"$periodic 3\n$cell{unit}\n4 5 6 80 90 100\n$coord{unit}\n1 2 3 si f\n$end\n"
    )
    reference = ch.parse_turbomole(text)
    result = ch.turbomole_to_poscar_checked(text)
    output = ch.PeriodicStructure.from_poscar(result.text)
    np.testing.assert_allclose(
        output.lattice.matrix, reference["lattice_angstrom"], atol=1e-14
    )
    np.testing.assert_allclose(
        output.cartesian_positions(), reference["positions_angstrom"], atol=1e-13
    )


@pytest.mark.parametrize(
    "scale,mode,position",
    [
        ("2", "Cartesian", "1 2 3"),
        ("-120", "Direct", "-0.25 0.5 1.25"),
        ("2 3 4", "Cartesian", "1 2 3"),
        ("1", "Direct", "0.5 0 0"),
    ],
)
def test_poscar_native_scale_rules(scale, mode, position):
    text = f"test\n{scale}\n4 0.5 0\n0 5 0.75\n0.25 0 6\nSi\n1\n{mode}\n{position}\n"
    original = ch.PeriodicStructure.from_poscar(text)
    result = ch.poscar_to_turbomole_checked(text)
    output = ch.parse_turbomole(result.text)
    np.testing.assert_allclose(
        output["lattice_angstrom"], original.lattice.matrix, atol=1e-13
    )
    np.testing.assert_allclose(
        output["positions_angstrom"], original.cartesian_positions(), atol=1e-13
    )
    assert output["frozen"] == [False]


def test_duplicate_poscar_species_groups_have_explicit_permutation():
    text = "test\n1\n4 0 0\n0 5 0\n0 0 6\nO Si O\n1 1 1\nSelective dynamics\nDirect\n0 0 0 F F F\n0.1 0.2 0.3 T T T\n0.4 0.5 0.6 T T T\n"
    result = ch.poscar_to_turbomole_checked(text)
    assert result.source_indices == (0, 2, 1)
    geometry = ch.parse_turbomole(result.text)
    assert geometry["elements"] == ["O", "O", "Si"]
    assert geometry["frozen"] == [True, False, False]
    np.testing.assert_allclose(
        geometry["positions_angstrom"], [[0, 0, 0], [1.6, 2.5, 3.6], [0.4, 1, 1.8]]
    )


@pytest.mark.parametrize(
    "bad",
    [
        "$coord\n0 0 0 si\n$end\n",
        TMOL.replace("$periodic 3", "$periodic 2"),
        TMOL.replace("$end", "$energy\n0\n$end"),
        TMOL.replace("4 0.5 0", "0 0 0"),
        TMOL.replace("-0.25", "nan"),
        TMOL.replace("$coord frac", "$coord file=coord"),
        TMOL.replace("$end", ""),
        TMOL.replace("$end", "$lattice angs\n4 0 0\n0 5 0\n0 0 6\n$end"),
    ],
)
def test_turbomole_unsupported_inputs_refuse(bad):
    with pytest.raises(ValueError):
        ch.turbomole_to_poscar_checked(bad)


def test_partial_constraints_and_restart_data_refuse():
    text = ch.turbomole_to_poscar_checked(TMOL).text
    with pytest.raises(ValueError, match="fully"):
        ch.poscar_to_turbomole_checked(text.replace("F F F", "F T F"))
    velocities = "\n" + "0 0 0\n" * 3
    with pytest.raises(ValueError, match="restart"):
        ch.poscar_to_turbomole_checked(text + velocities)


@pytest.mark.parametrize(
    "count", ["50001", "18446744073709551616", "9" * 100, "-1", "1.5", "０", "0"]
)
def test_count_preflight_before_native_allocation(count, monkeypatch):
    import chematic._periodic_geometry as implementation

    class Trap:
        @staticmethod
        def from_poscar(text):
            pytest.fail("native allocation must not start")

    monkeypatch.setattr(implementation, "PeriodicStructure", Trap)
    with pytest.raises(ValueError):
        ch.poscar_to_turbomole_checked(
            f"test\n1\n1 0 0\n0 1 0\n0 0 1\nSi\n{count}\nDirect\n0 0 0\n"
        )


@pytest.mark.parametrize(
    "convert", [ch.turbomole_to_poscar_checked, ch.poscar_to_turbomole_checked]
)
def test_text_and_byte_limits(convert):
    for value in (None, b"text", 1, []):
        with pytest.raises(TypeError):
            convert(value)
    for value in ("\ud800", "x" * 1025, "x" * ((16 << 20) + 1), "é" * (9 << 20)):
        with pytest.raises(ValueError):
            convert(value)


def test_crlf_and_comments_remain_in_source_sidecar():
    text = TMOL.replace("\n", "\r\n")
    assert ch.turbomole_to_poscar_checked(text).source_text == text
