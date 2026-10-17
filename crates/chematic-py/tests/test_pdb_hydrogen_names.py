"""Original coordinate-only naming fixtures; no hydrogen or bond perception."""

from dataclasses import FrozenInstanceError

import numpy as np
import pytest

import chematic as ch


def atom(serial, name=" H  ", element="H", *, kind="ATOM", chain="A",
         residue=1, insertion=" ", resname="LIG", alt=" "):
    return (f"{kind:6}{serial:5d} {name}{alt}{resname:>3} {chain}{residue:4d}"
            f"{insertion}   {serial / 10:8.3f}{-1.:8.3f}{2.:8.3f}"
            f"{1.:6.2f}{12.34:6.2f}          {element:>2}  \n")


def anisou(line):
    return "ANISOU" + line[6:28] + "".join(f"{n:7d}" for n in (100, 200, 300, -10, 20, -30)) + line[70:]


def example():
    first = atom(7, " O  ", "O", kind="HETATM")
    h = atom(12)
    return (first + h + anisou(h) + atom(20, "H1  ", "C") + atom(31)
            + atom(40, " D1 ", "D") + atom(51, "KEEP", "H", kind="HETATM")
            + "CONECT    7   12   31\nCONECT   12    7\nEND\n")


@pytest.mark.parametrize("crlf", [False, True])
@pytest.mark.parametrize("final_newline", [False, True])
def test_only_atom_and_associated_anisou_name_columns_change(crlf, final_newline):
    source = example()
    if not final_newline:
        source = source.rstrip("\n")
    if crlf:
        source = source.replace("\n", "\r\n")
    result = ch.name_pdb_hydrogens_checked(source)
    assert result.source_text == source and result.atom_count == 6
    assert result.renamed_count == 3
    assert [a.source_index for a in result.assignments] == [1, 3, 5]
    assert [a.serial for a in result.assignments] == [12, 31, 51]
    assert [a.new_name.strip() for a in result.assignments] == ["H2", "H3", "H4"]
    assert result.updated_anisou_lines == (3,)
    before, after = source.splitlines(keepends=True), result.pdb.splitlines(keepends=True)
    for index, (a, b) in enumerate(zip(before, after)):
        if index in (1, 2, 4, 6):
            assert a[:12] == b[:12] and a[16:] == b[16:]
        else:
            assert a == b
    assert len(result.pdb) == len(source)
    assert ch.name_pdb_hydrogens_checked(result.pdb).pdb == result.pdb


def test_selected_indices_reserve_other_hydrogens_and_use_source_order():
    source = atom(10, " H2 ") + atom(20) + atom(30, " H1 ", "C") + atom(40) + "END\n"
    result = ch.name_pdb_hydrogens_checked(source, [np.int64(3), 1])
    assert [a.source_index for a in result.assignments] == [1, 3]
    assert [a.new_name.strip() for a in result.assignments] == ["H3", "H4"]
    assert result.pdb.splitlines()[0] == source.splitlines()[0]


@pytest.mark.parametrize("change", [dict(chain="B"), dict(residue=-2),
                                    dict(insertion="I"), dict(resname="HOH")])
def test_residue_identity_restarts_names(change):
    result = ch.name_pdb_hydrogens_checked(atom(1) + atom(2, **change) + "END\n")
    assert [a.new_name.strip() for a in result.assignments] == ["H1", "H1"]


def test_ter_segment_restarts_names_and_is_preserved():
    source = atom(1) + "TER       2      LIG A   1 \n" + atom(3) + "END\n"
    result = ch.name_pdb_hydrogens_checked(source)
    assert [a.new_name.strip() for a in result.assignments] == ["H1", "H1"]
    assert result.pdb.splitlines()[1] == source.splitlines()[1]


def test_noncontiguous_same_residue_and_numeric_spelling_share_names():
    one = atom(1)
    other = atom(2, chain="B")
    last = atom(3)
    last = last[:22] + "0001" + last[26:]
    result = ch.name_pdb_hydrogens_checked(one + other + last + "END\n")
    assert [a.new_name.strip() for a in result.assignments] == ["H1", "H1", "H2"]


@pytest.mark.parametrize("count", [9, 10, 99, 100, 999])
def test_four_column_name_boundaries(count):
    source = "".join(atom(n + 1) for n in range(count)) + "END\n"
    result = ch.name_pdb_hydrogens_checked(source)
    field = result.assignments[-1].new_name
    assert len(field) == 4 and field.strip() == f"H{count}"
    assert field.startswith(" ") == (count < 100)


def test_exhaustion_and_retained_name_collisions_refuse():
    source = "".join(atom(n + 1) for n in range(1000)) + "END\n"
    with pytest.raises(ValueError, match="exhausted"):
        ch.name_pdb_hydrogens_checked(source)
    retained = "".join(atom(n, f"H{n}".ljust(4), "C") for n in range(1, 1000))
    with pytest.raises(ValueError, match="exhausted"):
        ch.name_pdb_hydrogens_checked(retained + atom(1001) + "END\n")


def test_no_hydrogens_empty_selection_and_unchanged_assignment():
    source = atom(1, " H1 ", "C") + "END\n"
    result = ch.name_pdb_hydrogens_checked(source)
    assert result.pdb == source and result.assignments == () and result.renamed_count == 0
    assert ch.name_pdb_hydrogens_checked(example(), []).pdb == example()
    result = ch.name_pdb_hydrogens_checked(atom(1, " H1 ") + "END\n")
    assert len(result.assignments) == 1 and result.renamed_count == 0


def test_isotope_element_records_are_not_inferred_as_h():
    source = atom(1, " D1 ", "D") + atom(2, " T1 ", "T") + atom(3) + "END\n"
    result = ch.name_pdb_hydrogens_checked(source)
    assert [a.source_index for a in result.assignments] == [2]
    assert result.pdb.splitlines()[:2] == source.splitlines()[:2]
    with pytest.raises(ValueError, match="explicit element H"):
        ch.name_pdb_hydrogens_checked(source, [0])


def test_selected_altloc_refuses_but_unselected_alternatives_are_preserved():
    source = atom(1, alt="A") + "END\n"
    with pytest.raises(ValueError, match="alternate"):
        ch.name_pdb_hydrogens_checked(source)
    assert ch.name_pdb_hydrogens_checked(source, []).pdb == source
    source = atom(1, " H1 ", "C", alt="A") + atom(2) + "END\n"
    assert ch.name_pdb_hydrogens_checked(source).assignments[0].new_name.strip() == "H2"


def test_single_explicit_model_and_mixed_line_endings_stay_exact():
    source = "MODEL        7\r\n" + atom(1) + "ENDMDL\r\nEND"
    result = ch.name_pdb_hydrogens_checked(source)
    assert result.pdb == source.replace(" H  ", " H1 ", 1)
    with pytest.raises(ValueError):
        ch.name_pdb_hydrogens_checked(source[:-3] + source)


@pytest.mark.parametrize("text", [None, b"", 1, False])
def test_text_type(text):
    with pytest.raises(TypeError):
        ch.name_pdb_hydrogens_checked(text)


@pytest.mark.parametrize("indices,error", [([0], ValueError), ([1, 1], ValueError),
    ([-1], ValueError), ([6], ValueError), ([True], TypeError), ([1.0], TypeError),
    ([np.bool_(True)], TypeError), ("1", TypeError), (b"1", TypeError),
    ({1}, TypeError), (iter([1]), TypeError), ([1] * 7, ValueError)])
def test_invalid_selection(indices, error):
    with pytest.raises(error):
        ch.name_pdb_hydrogens_checked(example(), indices)


@pytest.mark.parametrize("source", ["", "END\n", "HEADER original\n" + atom(1),
    "LINK  original\n" + atom(1), "SIGATM original\n" + atom(1),
    atom(1) + atom(1), atom(1).replace("   0.100", "     nan"),
    atom(1) + "CONECT    1   99\nEND\n", atom(1) + "END\n" + atom(2),
    atom(1, element=""), atom(1)[:54] + "\n", "REMARK λ\n" + atom(1)])
def test_invalid_or_unsupported_source_refuses(source):
    with pytest.raises(ValueError):
        ch.name_pdb_hydrogens_checked(source)


def test_bad_anisou_identity_refuses_and_reports_are_frozen():
    source = atom(1)
    with pytest.raises(ValueError):
        ch.name_pdb_hydrogens_checked(source + anisou(source).replace(" H  ", " X  "))
    report = ch.name_pdb_hydrogens_checked(source)
    with pytest.raises(FrozenInstanceError):
        report.pdb = ""
    with pytest.raises(FrozenInstanceError):
        report.assignments[0].new_name = "oops"


@pytest.mark.parametrize("limit,value", [("_MAX_BYTES", 10), ("_MAX_ATOMS", 1),
                                         ("_MAX_SUFFIX", 1)])
def test_resource_caps(limit, value, monkeypatch):
    import chematic._pdb_hydrogen_names as implementation

    monkeypatch.setattr(implementation, limit, value)
    with pytest.raises(ValueError):
        ch.name_pdb_hydrogens_checked(atom(1) + atom(2) + "END\n")
