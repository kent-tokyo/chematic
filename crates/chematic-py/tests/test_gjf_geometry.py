"""GJF geometry values and job metadata must be observable without silent row loss."""
import json
from pathlib import Path

import pytest
import chematic


def gjf(rows, charge_mult="-1 2"):
    return f"%chk=probe.chk\n# HF/STO-3G\n\ngeometry\n\n{charge_mult}\n{rows}\n\n"


def test_cartesian_values_and_optional_freeze_codes_are_separate():
    text = gjf("C1 -1 1.25 2.5 3.75\n8 0 -4.5 5.25 -6.75\nCl12 7.25 -8.5 9.75")
    g = chematic.parse_gjf_geometry(text)
    assert g["coords"] == [[1.25, 2.5, 3.75], [-4.5, 5.25, -6.75], [7.25, -8.5, 9.75]]
    assert g["freeze_codes"] == [-1, 0, None]
    assert (g["charge"], g["multiplicity"]) == (-1, 2)
    assert g["mol"].atom_count == 3
    assert chematic.from_gjf(text).atom_count == 3
    assert "C" in g["mol"].to_xyz(g["coords"])


def test_commas_comments_crlf_and_input_units():
    text = gjf("C, -1, 1.2345678901234567, -2.5e-20, 3.75 ! fixed\n! comment\nH 0 2.5 3.5 4.5", "-1,2 ! spin")
    for value in (text, text.replace("\n", "\r\n"), text.replace("# HF/STO-3G", "# HF/STO-3G Units=Bohr")):
        g = chematic.parse_gjf_geometry(value)
        assert g["coords"] == [[1.2345678901234567, -2.5e-20, 3.75], [2.5, 3.5, 4.5]]
        assert g["freeze_codes"] == [-1, 0]


@pytest.mark.parametrize("bad", [
    "H", "H 1", "H 1 2", "H -1 1 2 3 4", "H 0.5 1 2 3", "H x 1 2 3",
    "H 2147483648 1 2 3", "Variables:", "Constants:", "C 0 1 2 3 H",
    "C -1 1 2 NaN", "C -1 1 2 inf", "C -1 1 2 1e999", "Xx 1 2 3",
])
@pytest.mark.parametrize("position", ["first", "last"])
def test_malformed_rows_do_not_return_a_partial_molecule(bad, position):
    rows = f"{bad}\nC 0 0 0" if position == "first" else f"C 0 0 0\n{bad}"
    text = gjf(rows)
    with pytest.raises(ValueError):
        chematic.parse_gjf_geometry(text)
    with pytest.raises(ValueError):
        chematic.from_gjf(text)
    with pytest.raises(ValueError):
        chematic.convert_format(text, "gjf", "cjson")


@pytest.mark.parametrize("cm", ["0 0", "0 -1", "0 1 0 1", "0", "a 1"])
def test_invalid_or_fragment_charge_spin_is_refused(cm):
    with pytest.raises(ValueError):
        chematic.parse_gjf_geometry(gjf("C 0 0 0", cm))


def test_multiple_jobs_require_explicit_selection_outside_this_reader():
    text = gjf("C 0 0 0")
    with pytest.raises(ValueError, match="Link1"):
        chematic.parse_gjf_geometry(text + "--Link1--\n" + text)


def test_input_and_atom_limits_are_enforced_at_the_boundary():
    text = gjf("C 0 1 2 3\nH 4 5 6")
    assert len(chematic.parse_gjf_geometry(text, max_atoms=2, max_input_bytes=len(text))["coords"]) == 2
    for kwargs in ({"max_atoms": 1}, {"max_input_bytes": len(text)-1}, {"max_atoms": 0}):
        with pytest.raises(ValueError, match="exceeding limit"):
            chematic.parse_gjf_geometry(text, **kwargs)


_CASES = json.loads((Path(__file__).parents[3] / "validation/gjf-cartesian-2026-10-22/cases.json").read_text(encoding="utf-8"))["cases"]


@pytest.mark.parametrize("case", _CASES, ids=lambda case: case["id"])
@pytest.mark.parametrize("newline", ["LF", "CRLF"])
def test_original_cartesian_packet(case, newline):
    text = case["input"].replace("\n", "\r\n") if newline == "CRLF" else case["input"]
    expected = case["expected"]
    if expected is None:
        with pytest.raises(ValueError):
            chematic.parse_gjf_geometry(text)
        return
    g = chematic.parse_gjf_geometry(text)
    assert {key: g[key] for key in ("coords", "freeze_codes", "charge", "multiplicity")} == {
        key: expected[key] for key in ("coords", "freeze_codes", "charge", "multiplicity")
    }
    assert [row.split()[0] for row in g["mol"].to_xyz(g["coords"]).splitlines()[2:]] == expected["elements"]
