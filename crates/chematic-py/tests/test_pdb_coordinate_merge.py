"""Synthetic source-preserving PDB assembly and reference-remapping contracts."""

from collections.abc import Mapping
from dataclasses import FrozenInstanceError
import itertools
import string

import pytest

from chematic import merge_pdb_coordinate_blocks as merge, select_pdb_atoms


def atom(
    serial=1,
    chain="A",
    kind="ATOM",
    name=" CA ",
    residue="ALA",
    number=7,
    alt=" ",
    insertion=" ",
    element="C",
    charge="  ",
    x=1.25,
):
    return f"{kind:<6}{serial:5} {name:4}{alt}{residue:3} {chain}{number:4}{insertion}   {x:8.3f}{-2.5:8.3f}{3.75:8.3f}{.5:6.2f}{23.45:6.2f}          {element:>2}{charge:2}"


def anisou(line):
    return (
        "ANISOU"
        + line[6:28]
        + "".join(f"{v:7}" for v in [100, 200, 300, -4, 5, -6])
        + line[70:]
    )


def ter(line, serial=None):
    return f"TER   {serial if serial is not None else int(line[6:11])+1:5}      {line[17:27]}"


def block(*lines):
    return "\n".join(lines) + "\nEND\n"


def atom_lines(text):
    return [line for line in text.splitlines() if line.startswith(("ATOM  ", "HETATM"))]


def test_protein_ligand_merge_keeps_record_kind_ca_element_and_source_boundaries():
    protein = block(atom(20), atom(35, name=" N  ", element="N"))
    ligand = block(atom(20, chain="L", kind="HETATM", residue="LIG", name=" C1 "))
    result = merge([protein, ligand])
    lines = result.pdb.splitlines()
    assert [s[:6].strip() for s in lines] == [
        "ATOM",
        "ATOM",
        "TER",
        "HETATM",
        "TER",
        "END",
    ]
    assert [int(s[6:11]) for s in lines[:-1]] == [1, 2, 3, 4, 5]
    assert result.inserted_ter_after_sources == (0, 1)
    assert result.source_texts == (protein, ligand)
    expected = atom_lines(protein) + atom_lines(ligand)
    for actual, original in zip(atom_lines(result.pdb), expected):
        assert actual[:6] + actual[11:] == original[:6] + original[11:]
    assert atom_lines(result.pdb)[0][12:16] == " CA "
    assert atom_lines(result.pdb)[0][76:78] == " C"
    assert [
        (a.source_index, a.source_atom_index, a.old_serial, a.new_serial, a.record_type)
        for a in result.atoms
    ] == [(0, 0, 20, 1, "ATOM"), (0, 1, 35, 2, "ATOM"), (1, 0, 20, 4, "HETATM")]
    assert select_pdb_atoms(result.pdb)["pdb"] == result.pdb


def test_anisou_ter_conect_references_follow_renumbering_and_chain_mapping():
    first = atom(40, alt="A", insertion="I", number=-2, charge="1+")
    second = atom(70, name=" CB ")
    source = block(
        first,
        anisou(first),
        second,
        ter(second),
        "CONECT   40   70   70",
        "CONECT   70   40",
    )
    ligand = block(atom(40, kind="HETATM", residue="LIG"))
    maps = [{"A": "P"}, {"A": "L"}]
    result = merge([source, ligand], chain_maps=maps)
    assert result.inserted_ter_after_sources == (1,)
    lines = result.pdb.splitlines()
    assert lines[1][:6] == "ANISOU" and int(lines[1][6:11]) == 1
    assert lines[1][21] == "P" and lines[1][28:] == anisou(first)[28:]
    assert lines[0][16:21] == first[16:21] and lines[0][22:] == first[22:]
    assert [s for s in lines if s.startswith("CONECT")] == [
        "CONECT    1    2    2",
        "CONECT    2    1",
    ]
    assert lines.index("CONECT    1    2    2") > max(
        i for i, s in enumerate(lines) if s.startswith("TER")
    )
    assert [(a.source_chain, a.output_chain) for a in result.atoms] == [
        ("A", "P"),
        ("A", "P"),
        ("A", "L"),
    ]
    maps[0]["A"] = "Z"
    assert result.atoms[0].output_chain == "P"


def test_multiple_source_chains_blank_chain_and_input_order_are_preserved():
    source = block(atom(1, chain=" "), atom(9, chain="B"))
    ligand = block(atom(1, chain="L", kind="HETATM"))
    result = merge((s for s in [source, ligand]), chain_maps=({" ": "A"}, {}))
    assert [a.output_chain for a in result.atoms] == ["A", "B", "L"]
    assert [line[21] for line in atom_lines(result.pdb)] == ["A", "B", "L"]


def test_crlf_unterminated_input_and_existing_final_ter():
    first = atom(7)
    source = first + "\r\n" + ter(first) + "\r\n"
    result = merge([source])
    assert result.source_texts == (source,)
    assert "\r" not in result.pdb and result.pdb.endswith("END\n")
    assert result.inserted_ter_after_sources == ()
    assert (
        len([line for line in result.pdb.splitlines() if line.startswith("TER")]) == 1
    )


def test_metadata_loss_requires_explicit_consent_and_retains_original_lines():
    source = "HEADER    original\nSEQRES    original\n" + block(atom())
    with pytest.raises(ValueError, match="consent"):
        merge([source])
    result = merge([source], allow_metadata_loss=True)
    assert result.source_texts == (source,)
    assert [
        (loss.source_index, loss.line, loss.record_name) for loss in result.losses
    ] == [(0, 1, "HEADER"), (0, 2, "SEQRES")]
    assert "HEADER" not in result.pdb
    with pytest.raises(FrozenInstanceError):
        result.pdb = "changed"
    with pytest.raises(FrozenInstanceError):
        result.atoms[0].new_serial = 3


@pytest.mark.parametrize(
    "maps,error",
    [
        (None, ValueError),
        ([{}, {}], ValueError),
        ([{"A": "Q"}, {"A": "Q"}], ValueError),
        ([{"Z": "Q"}, {"A": "L"}], ValueError),
        (
            [
                {},
            ],
            ValueError,
        ),
        ([{}, {}, {}], ValueError),
        (["A", {}], TypeError),
        ([{"A": ""}, {}], ValueError),
        ([{"A": "AA"}, {}], ValueError),
        ([{"A": "!"}, {}], ValueError),
        ([{"A": 1}, {}], TypeError),
        ([{1: "A"}, {}], TypeError),
    ],
)
def test_colliding_or_malformed_chain_maps_refuse(maps, error):
    with pytest.raises(error):
        merge([block(atom()), block(atom())], chain_maps=maps)


def test_noninjective_mapping_within_one_source_refuses():
    with pytest.raises(ValueError, match="collide"):
        merge(
            [block(atom(1, chain="A"), atom(2, chain="B"))],
            chain_maps=[{"A": "C", "B": "C"}],
        )


@pytest.mark.parametrize("flag", [1, "yes", None])
def test_metadata_consent_requires_bool(flag):
    with pytest.raises(TypeError):
        merge([block(atom())], allow_metadata_loss=flag)


@pytest.mark.parametrize(
    "sources,error",
    [
        (None, TypeError),
        ("PDB", TypeError),
        ([b"PDB"], TypeError),
        ([1], TypeError),
        ([], ValueError),
        ([""], ValueError),
        (["END\n"], ValueError),
        ([block(atom(chain="!"))], ValueError),
        (["REMARK λ\n"], ValueError),
    ],
)
def test_source_types_and_empty_or_unsupported_chains(sources, error):
    with pytest.raises(error):
        merge(sources)


@pytest.mark.parametrize(
    "source",
    [
        block(atom(), atom()),
        block(atom(x=float("nan"))),
        block(atom()).replace("  0.50", "   inf"),
        block(atom(), "CONECT    1   99"),
        block(atom(), "CONECT    1    1"),
        block(atom(), "ANISOU    9"),
        block(atom(), ter(atom(), 99)),
        "MODEL        1\n" + atom() + "\nENDMDL\nEND\n",
        block(atom()) + atom(),
        "REMARK" + "x" * 1025 + "\n" + block(atom()),
    ],
)
def test_native_validation_is_not_bypassed_by_metadata_consent(source):
    with pytest.raises(ValueError):
        merge([source], allow_metadata_loss=True)


def test_source_and_mapping_exceptions_propagate_before_native_calls(monkeypatch):
    import chematic._pdb_merge as impl

    def forbidden(*args, **kwargs):
        raise AssertionError("native validation must not start")

    monkeypatch.setattr(impl, "_select", forbidden)
    error = RuntimeError("source interrupted")

    def sources():
        yield block(atom())
        raise error

    class Broken(Mapping):
        def __len__(self):
            return 1

        def __iter__(self):
            return iter(["A"])

        def __getitem__(self, key):
            raise error

    for source, maps in [(sources(), None), ([block(atom())], [Broken()])]:
        with pytest.raises(RuntimeError) as caught:
            merge(source, chain_maps=maps)
        assert caught.value is error
    with pytest.raises(ValueError, match="32 entries"):
        merge([block(atom())], chain_maps=itertools.repeat({}))


def test_source_count_byte_record_and_atom_caps():
    source = block(atom())
    chains = string.ascii_letters[:32]
    result = merge(itertools.repeat(source, 32), chain_maps=({"A": c} for c in chains))
    assert len(result.atoms) == 32 and len(result.source_texts) == 32
    with pytest.raises(ValueError, match="32 sources"):
        merge(itertools.repeat(source))
    with pytest.raises(ValueError, match="8 MiB"):
        merge(["x" * (8 * 1024 * 1024 + 1)])
    with pytest.raises(ValueError, match="aggregate"):
        merge(["x" * (8 * 1024 * 1024)] * 3)
    with pytest.raises(ValueError, match="aggregate"):
        merge(["\n" * 250001])
    exact_records = atom() + "\n" + "\n" * 249998 + "END\n"
    assert len(merge([exact_records]).atoms) == 1
    large = block(*(atom(i + 1) for i in range(50000)))
    assert len(merge([large]).atoms) == 50000
    with pytest.raises(ValueError):
        merge([large, block(atom(chain="B"))])


def test_actual_five_column_serial_boundary_and_overflow():
    # 50000 atoms + 49999 TER: the final permitted serial is exactly 99999.
    lines = [atom(1), atom(2), ter(atom(2))]
    for serial in range(4, 99999, 2):
        row = atom(serial)
        lines.extend([row, ter(row)])
    result = merge([block(*lines)])
    assert len(result.atoms) == 50000
    assert (
        max(int(s[6:11]) for s in result.pdb.splitlines() if s.startswith("TER"))
        == 99999
    )
    # Every atom is a separate segment: the final source atom 99999 needs TER 100000.
    lines = []
    for serial in range(1, 99999, 2):
        row = atom(serial)
        lines.extend([row, ter(row)])
    lines.append(atom(99999))
    with pytest.raises(ValueError, match="five-column"):
        merge([block(*lines)])
