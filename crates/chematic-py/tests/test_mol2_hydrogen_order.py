"""Source-aware MOL2 H ordering and complete supported ID remapping."""

from dataclasses import FrozenInstanceError
import json

import pytest

import chematic as ch
import chematic._mol2_hydrogen_order as implementation


def ethanol(ending="\n"):
    rows = [
        "@<TRIPOS>MOLECULE",
        "ethanol",
        "9 8 1 0 0",
        "SMALL",
        "USER_CHARGES",
        "",
        "@<TRIPOS>ATOM",
        " 10 C1   0.000000 0 0 C.3  1 RES  -0.123456 BACKBONE",
        " 20 C2   1.500000 0 0 C.3  1 RES   0.234567",
        " 30 O1   2.800000 0 0 O.3  1 RES  -0.345678",
        " 40 HA   -0.5  0.9 0 H 1 RES 0.012345",
        " 50 HB   -0.5 -0.9 0 H 1 RES 0.012345",
        " 60 HC   0 0 1 H 1 RES 0.012345",
        " 70 HD   1.5 1 0 H 1 RES 0.012345",
        " 80 HE   1.5 0 1 H 1 RES 0.012345",
        " 90 HF   3.2 0.9 0 H 1 RES 0.172842",
        "@<TRIPOS>BOND",
        "1 10 20 1",
        "2 20 30 1",
        "3 10 40 1",
        "4 10 50 1",
        "5 10 60 1",
        "6 20 70 1",
        "7 20 80 1",
        "8 30 90 1",
        "@<TRIPOS>SUBSTRUCTURE",
        "1 RES 20 GROUP 0 **** **** 0 ROOT",
        "@<TRIPOS>UNITY_ATOM_ATTR",
        "30 1",
        "charge 0",
    ]
    return ending.join(rows) + ending


def inspect(text):
    return json.loads(ch.parse_mol2_record_json(text))


@pytest.mark.parametrize("ending", ["\n", "\r\n"])
def test_h_rows_follow_parent_and_every_reference_is_remapped(ending):
    source = ethanol(ending)
    result = ch.order_mol2_hydrogens_checked(source)
    assert result.source_text == source
    assert result.new_to_old == (0, 3, 4, 5, 1, 6, 7, 2, 8)
    assert result.old_to_new == (0, 4, 7, 1, 2, 3, 5, 6, 8)
    assert result.atom_id_map == (
        (10, 1),
        (20, 5),
        (30, 8),
        (40, 2),
        (50, 3),
        (60, 4),
        (70, 6),
        (80, 7),
        (90, 9),
    )
    assert result.hydrogen_parents == (
        (40, 10),
        (50, 10),
        (60, 10),
        (70, 20),
        (80, 20),
        (90, 30),
    )
    assert "1 RES 5 GROUP 0 **** **** 0 ROOT" in result.mol2
    assert "@<TRIPOS>UNITY_ATOM_ATTR" + ending + "8 1" in result.mol2
    assert "1 1 5 1" in result.mol2
    assert "8 8 9 1" in result.mol2
    before, after = inspect(source), inspect(result.mol2)
    for new, old in enumerate(result.new_to_old):
        assert after["coordinates"][new] == before["coordinates"][old]
        a, b = dict(after["atoms"][new]), dict(before["atoms"][old])
        del a["source_id"], b["source_id"]
        assert a == b
        assert after["formal_charges"][new] == before["formal_charges"][old]
    assert result.mol2.count(ending) == source.count(ending)
    assert "-0.123456 BACKBONE" in result.mol2
    with pytest.raises(FrozenInstanceError):
        result.mol2 = "changed"


def test_idempotent_and_already_grouped_source_is_exact():
    once = ch.order_mol2_hydrogens_checked(ethanol())
    twice = ch.order_mol2_hydrogens_checked(once.mol2)
    assert twice.mol2 == once.mol2
    assert twice.new_to_old == twice.old_to_new == tuple(range(9))
    assert twice.atom_id_map == tuple((i, i) for i in range(1, 10))


def test_no_h_leaves_arbitrary_ids_and_text_exact():
    source = "@<TRIPOS>MOLECULE\nions\n2 0 0 0 0\nSMALL\nNO_CHARGES\n@<TRIPOS>ATOM\n91 Na 0 0 0 Na\n7 Cl 9 0 0 Cl\n"
    result = ch.order_mol2_hydrogens_checked(source)
    assert result.mol2 == source
    assert result.atom_id_map == ((91, 91), (7, 7))
    assert result.hydrogen_parents == ()


def test_formal_charges_follow_new_indices_and_partial_charges_stay_separate():
    source = (
        "@<TRIPOS>MOLECULE\nion\n6 4 0 0 0\nSMALL\nUSER_CHARGES\n"
        "@<TRIPOS>ATOM\n1 N 0 0 0 N.4 1 ION -0.7\n2 Cl 9 0 0 Cl 1 ION 0.7\n"
        "3 H1 1 1 1 H 1 ION 0\n4 H2 -1 -1 1 H 1 ION 0\n5 H3 -1 1 -1 H 1 ION 0\n6 H4 1 -1 -1 H 1 ION 0\n"
        "@<TRIPOS>BOND\n1 1 3 1\n2 1 4 1\n3 1 5 1\n4 1 6 1\n"
        "@<TRIPOS>UNITY_ATOM_ATTR\n1 1\ncharge 1\n2 1\ncharge -1\n"
    )
    result = ch.order_mol2_hydrogens_checked(source)
    assert inspect(result.mol2)["formal_charges"] == [1, 0, 0, 0, 0, -1]
    assert "1 N 0 0 0 N.4 1 ION -0.7" in result.mol2
    assert "6 Cl 9 0 0 Cl 1 ION 0.7" in result.mol2


@pytest.mark.parametrize("section", ["SET", "UNITY_BOND_ATTR", "CRYSIN", "UNKNOWN"])
def test_opaque_sections_refuse_even_without_reordering(section):
    source = ethanol() + "@<TRIPOS>" + section + "\n1 2 3\n"
    with pytest.raises(ValueError, match="unsupported.*section"):
        ch.order_mol2_hydrogens_checked(source)


@pytest.mark.parametrize(
    "change",
    [
        lambda s: s.replace("charge 0", "color red"),
        lambda s: s.replace("30 1\ncharge 0", "30 2\ncharge 0\nlabel abc"),
        lambda s: s.replace("1 RES 20 GROUP", "1 RES 999 GROUP"),
        lambda s: s.replace("1 RES 20 GROUP", "bad RES 20 GROUP"),
        lambda s: s.replace("1 RES 20 GROUP", "2 RES 20 GROUP"),
        lambda s: s.replace(
            "@<TRIPOS>SUBSTRUCTURE\n", "@<TRIPOS>SUBSTRUCTURE\n1 RES 10 GROUP\n"
        ),
        lambda s: s.replace("3 10 40 1", "3 50 40 1"),
        lambda s: s.replace("3 10 40 1", "3 10 40 2"),
        lambda s: s.replace("9 8 1 0 0", "9 7 1 0 0").replace("3 10 40 1\n", ""),
        lambda s: s.replace("8 30 90 1", "8 30 999 1"),
        lambda s: s.replace(" 90 HF", " 10 HF"),
        lambda s: s + s,
    ],
)
def test_unsupported_or_invalid_sources_refuse(change):
    with pytest.raises(ValueError):
        ch.order_mol2_hydrogens_checked(change(ethanol()))


@pytest.mark.parametrize(
    "source",
    ["", "not mol2", "\ud800", "a" * ((16 << 20) + 1)],
    ids=["empty", "invalid-mol2", "unpaired-surrogate", "oversize"],
)
def test_bad_text_or_byte_limits(source):
    with pytest.raises(ValueError):
        ch.order_mol2_hydrogens_checked(source)


@pytest.mark.parametrize("source", [None, b"mol2", 42])
def test_wrong_types(source):
    with pytest.raises(TypeError):
        ch.order_mol2_hydrogens_checked(source)


def test_comments_whitespace_and_missing_terminal_newline_outside_atoms():
    source = ethanol().replace(" 20 C2", "# between atom rows\n 20 C2").rstrip("\n")
    result = ch.order_mol2_hydrogens_checked(source)
    assert result.mol2.splitlines()[8] == "# between atom rows"
    assert not result.mol2.endswith("\n")
    assert result.mol2.count("# between atom rows") == 1


def test_inconsistent_atom_terminators_refuse():
    source = ethanol().replace(
        " 20 C2   1.500000 0 0 C.3  1 RES   0.234567\n",
        " 20 C2   1.500000 0 0 C.3  1 RES   0.234567\r\n",
    )
    with pytest.raises(ValueError, match="common LF or CRLF"):
        ch.order_mol2_hydrogens_checked(source)


def test_mapped_readback_disagreement_never_returns_a_result(monkeypatch):
    original = implementation.parse_mol2_record_json
    calls = []

    def changed(text):
        result = json.loads(original(text))
        calls.append(text)
        if len(calls) == 2:
            result["formal_charges"][0] = 7
        return json.dumps(result)

    monkeypatch.setattr(implementation, "parse_mol2_record_json", changed)
    with pytest.raises(RuntimeError, match="mapped readback"):
        ch.order_mol2_hydrogens_checked(ethanol())
    assert len(calls) == 2


def test_actual_atom_limit_and_rejection_before_native(monkeypatch):
    header = (
        "@<TRIPOS>MOLECULE\nhelium\n16384 0 0 0 0\nSMALL\nNO_CHARGES\n@<TRIPOS>ATOM\n"
    )
    source = header + "".join(f"{i} He{i} {i} 0 0 He\n" for i in range(1, 16385))
    result = ch.order_mol2_hydrogens_checked(source)
    assert result.mol2 == source
    assert len(result.new_to_old) == 16384

    def forbidden(text):
        pytest.fail("oversized input must not reach native parsing")

    monkeypatch.setattr(implementation, "parse_mol2_record_json", forbidden)
    with pytest.raises(ValueError, match="16384 atoms"):
        ch.order_mol2_hydrogens_checked(source + "16385 He 0 0 0 He\n")


def test_bond_and_physical_line_limits_precede_native(monkeypatch):
    def forbidden(text):
        pytest.fail("oversized input must not reach native parsing")

    monkeypatch.setattr(implementation, "parse_mol2_record_json", forbidden)
    head = "@<TRIPOS>MOLECULE\nexample\n1 32769\nSMALL\nNO_CHARGES\n@<TRIPOS>ATOM\n1 He 0 0 0 He\n@<TRIPOS>BOND\n"
    with pytest.raises(ValueError, match="32768 bonds"):
        ch.order_mol2_hydrogens_checked(head + "1 1 2 1\n" * 32769)
    with pytest.raises(ValueError, match="physical line limits"):
        ch.order_mol2_hydrogens_checked("#" + "x" * 16384 + "\n" + ethanol())
