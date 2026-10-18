"""Original PDB attributes, blank optional fields and immutable source views."""

from dataclasses import FrozenInstanceError
import itertools
from pathlib import Path

import numpy as np
import pytest
import chematic as C


def atom(
    serial=1,
    *,
    name=" CA ",
    alt="A",
    residue="VAL",
    chain="A",
    seq=101,
    insertion=" ",
    xyz=(1.125, -2.5, 3.75),
    occupancy="  0.50",
    b_factor=" 23.45",
    kind="ATOM",
    element=" C",
    charge="  ",
):
    return (
        f"{kind:6}{serial:5d} {name}{alt}{residue:3} {chain}{seq:4d}{insertion}   "
        + "".join(f"{value:8.3f}" for value in xyz)
        + occupancy
        + b_factor
        + " " * 10
        + element
        + charge
    )


def test_full_record_fields_and_hetero_atoms():
    first = atom()
    second = atom(
        2,
        name=" O1 ",
        alt="B",
        residue="LIG",
        chain="Z",
        seq=-2,
        insertion="I",
        xyz=(8, 9, 10),
        occupancy="  0.25",
        b_factor=" 45.67",
        kind="HETATM",
        element=" O",
        charge="1-",
    )
    text = "HEADER    original\n" + first + "\n" + second + "\nEND\n"
    result = C.read_pdb_atom_records(text)
    assert result.source_text == text and result.model_num is None
    assert result.unhandled_records == ((1, "HEADER"),)
    a, b = result.atoms
    assert isinstance(a, C.PdbAtomRecord) and isinstance(result, C.PdbAtomRead)
    assert (a.source_index, a.line_number, a.serial) == (0, 2, 1)
    assert (a.atom_name, a.alt_loc, a.residue_name, a.chain_id) == (
        " CA ",
        "A",
        "VAL",
        "A",
    )
    assert a.position_angstrom == (1.125, -2.5, 3.75)
    assert (a.occupancy, a.b_factor) == (0.5, 23.45)
    assert (b.record_type, b.residue_number, b.insertion_code) == ("HETATM", -2, "I")
    assert (b.element_field, b.charge_field) == (" O", "1-")
    assert b.record_text == second


def test_blank_optional_fields_differ_from_explicit_zero():
    text = (
        atom(1, occupancy=" " * 6, b_factor=" " * 6)
        + "\n"
        + atom(2, occupancy="  0.00", b_factor="  0.00")
        + "\nEND\n"
    )
    a, b = C.read_pdb_atom_records(text).atoms
    assert a.occupancy is None and a.b_factor is None
    assert b.occupancy == 0.0 and b.b_factor == 0.0
    short = C.read_pdb_atom_records(atom()[:54] + "\n").atoms[0]
    assert short.occupancy is None and short.b_factor is None
    assert short.element_field == short.charge_field == ""


@pytest.mark.parametrize("ending", ["\n", "\r\n"])
@pytest.mark.parametrize("model", [1, 7])
def test_single_model_source_line_endings_and_indices(ending, model):
    text = ending.join(
        [f"MODEL     {model:4d}", atom(), "", atom(2, alt="B"), "ENDMDL", "END", ""]
    )
    result = C.read_pdb_atom_records(text)
    assert result.source_text == text and result.model_num == model
    assert [a.line_number for a in result.atoms] == [2, 4]
    assert [a.source_index for a in result.atoms] == [0, 1]
    assert result.atoms[0].record_text == atom()


def test_all_720_atom_orders_preserve_correspondence():
    rows = [
        atom(
            i + 1,
            seq=i - 2,
            alt="AB"[i % 2],
            occupancy=f"{(i+1)/10:6.2f}",
            b_factor=f"{i*11.13:6.2f}",
            xyz=(i * 1.25, -i * 0.125, i * 0.75),
        )
        for i in range(6)
    ]
    for order in itertools.permutations(range(6)):
        text = "\n".join(rows[i] for i in order) + "\nEND\n"
        result = C.read_pdb_atom_records(text)
        assert len(result.atoms) == 6
        for index, a in enumerate(result.atoms):
            original = order[index]
            assert a.source_index == index and a.line_number == index + 1
            assert a.serial == original + 1
            assert a.record_text == rows[original]
            assert a.occupancy == pytest.approx((original + 1) / 10)
            assert a.position_angstrom == (
                original * 1.25,
                -original * 0.125,
                original * 0.75,
            )


def test_deep_immutability_and_repeat_read_ownership():
    text = atom() + "\nEND\n"
    first = C.read_pdb_atom_records(text)
    second = C.read_pdb_atom_records(text)
    assert first == second and first is not second
    assert first.atoms[0] is not second.atoms[0]
    for obj, name, value in [
        (first, "source_text", "bad"),
        (first.atoms[0], "occupancy", 0.0),
    ]:
        with pytest.raises(FrozenInstanceError):
            setattr(obj, name, value)
    assert isinstance(first.atoms, tuple)
    assert isinstance(first.atoms[0].position_angstrom, tuple)
    assert isinstance(first.unhandled_records, tuple)


def test_existing_anisou_ter_and_conect_source_is_retained():
    path = (
        Path(__file__).resolve().parents[3]
        / "validation/openbabel_issues/pdb-selection/rich.pdb"
    )
    text = path.read_text()
    result = C.read_pdb_atom_records(text)
    assert result.source_text == text
    assert len(result.atoms) == 6
    assert result.atoms[0].b_factor == 23.45
    assert "ANISOU" in text and "CONECT" in text and "TER" in text


def test_read_preserves_ter_and_metadata_without_export():
    # No selection text is exported; all source metadata remains available.
    text = (
        "REMARK    preserve me\n"
        + atom(1)
        + "\n"
        + f"TER   {2:5d}      VAL A 101 "
        + "\n"
        + atom(3, chain="B")
        + "\nEND\n"
    )
    result = C.read_pdb_atom_records(text)
    assert result.source_text == text
    assert [a.serial for a in result.atoms] == [1, 3]
    assert result.unhandled_records == ((1, "REMARK"),)


@pytest.mark.parametrize(
    "field,value",
    [
        ("occupancy", "   NaN"),
        ("occupancy", "   inf"),
        ("occupancy", "  junk"),
        ("b_factor", "   NaN"),
        ("b_factor", "  junk"),
    ],
)
def test_invalid_optional_numbers_fail_instead_of_becoming_defaults(field, value):
    with pytest.raises(ValueError, match="occupancy or temperature"):
        C.read_pdb_atom_records(atom(**{field: value}) + "\n")


@pytest.mark.parametrize("cut", [55, 59, 61, 65])
def test_nonblank_partial_optional_field_is_refused(cut):
    # Keep at least one nonblank character in an incomplete six-column field.
    line = atom(occupancy="123456", b_factor="123456")[:cut]
    with pytest.raises(ValueError, match="incomplete"):
        C.read_pdb_atom_records(line + "\n")


@pytest.mark.parametrize("bad", [None, b"ATOM", 12, True, []])
def test_contents_are_str_only(bad):
    with pytest.raises(TypeError):
        C.read_pdb_atom_records(bad)


@pytest.mark.parametrize(
    "name,maximum",
    [
        ("max_input_bytes", 8 * 1024 * 1024),
        ("max_atoms", 50000),
        ("max_records", 250000),
        ("max_line_len", 1024),
    ],
)
@pytest.mark.parametrize(
    "value,error",
    [
        (0, ValueError),
        (-1, ValueError),
        (True, TypeError),
        (1.5, TypeError),
        ("3", TypeError),
    ],
)
def test_invalid_limits(name, maximum, value, error):
    with pytest.raises(error):
        C.read_pdb_atom_records("", **{name: value})
    with pytest.raises(ValueError):
        C.read_pdb_atom_records("", **{name: maximum + 1})


@pytest.mark.parametrize(
    "kwargs",
    [
        {"max_input_bytes": 20},
        {"max_atoms": 1},
        {"max_records": 1},
        {"max_line_len": 40},
    ],
)
def test_lower_limits_and_exact_boundaries(kwargs):
    text = atom() + "\n" + atom(2) + "\nEND\n"
    with pytest.raises(ValueError):
        C.read_pdb_atom_records(text, **kwargs)
    assert (
        len(
            C.read_pdb_atom_records(
                text,
                max_input_bytes=np.int64(len(text)),
                max_atoms=np.int64(2),
                max_records=np.int64(3),
                max_line_len=np.int64(80),
            ).atoms
        )
        == 2
    )


@pytest.mark.parametrize(
    "text",
    [
        "REMARK λ\n",
        "MODEL        1\n"
        + atom()
        + "\nENDMDL\nMODEL        2\n"
        + atom()
        + "\nENDMDL\n",
        atom() + "\n" + atom() + "\n",
        atom() + "\nCONECT    1   99\nEND\n",
        atom() + "\nEND\n" + atom(2) + "\n",
    ],
)
def test_native_frame_serial_reference_and_ascii_checks(text):
    with pytest.raises(ValueError):
        C.read_pdb_atom_records(text)


def test_empty_read_and_uninterpreted_element_charge_fields():
    assert C.read_pdb_atom_records("").atoms == ()
    a = C.read_pdb_atom_records(atom(element="  ", charge="  ") + "\n").atoms[0]
    assert a.element_field == a.charge_field == "  "
    assert a.alt_loc == "A"
    assert (
        C.read_pdb_atom_records(atom(alt=" ", chain=" ", insertion=" ") + "\n")
        .atoms[0]
        .chain_id
        == " "
    )
