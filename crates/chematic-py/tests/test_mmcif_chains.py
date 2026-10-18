"""Explicit chain partitions, identity retention and bounded eager reads."""

from dataclasses import FrozenInstanceError
import itertools

import numpy as np
import pytest

import chematic as C

COLUMNS = """group_PDB id type_symbol label_atom_id label_alt_id label_comp_id
label_asym_id label_entity_id label_seq_id pdbx_PDB_ins_code Cartn_x Cartn_y
Cartn_z occupancy B_iso_or_equiv auth_seq_id auth_comp_id auth_asym_id
 auth_atom_id pdbx_PDB_model_num pdbx_formal_charge custom_payload""".split()

ROWS = [
    "ATOM 1 C CA . GLY Alpha 1 1 ? 0.125 -2.5 3.75 1 12.5 10 GLY A CA 1 ? opaque1".split(),
    "ATOM 2 N N . GLY Alpha 1 1 ? 1.25 -2 4 1 13.5 10 GLY A N 1 ? opaque2".split(),
    "HETATM 3 N N1 AA A1AJ1 Ligand 2 . AB 7.125 0 1 0.4 22.5 50 AUTHOR_LONG A AN1 1 1 opaque3".split(),
    "HETATM 4 N N1 BB A1AJ1 Ligand 2 . AB 7.25 0 1 0.6 23.5 50 AUTHOR_LONG A AN1 1 1 opaque4".split(),
    "ATOM 1 C CA . GLY Alpha 1 1 ? -1 8 2 1 12.5 10 GLY A CA 2 ? opaque5".split(),
    "ATOM 2 N N . GLY Alpha 1 1 ? -2 8 3 1 13.5 10 GLY A N 2 ? opaque6".split(),
    "HETATM 5 O O . HOH Solvent 3 . ? 15 20 25 1 30 70 HOH A O 2 ? opaque7".split(),
    "HETATM 6 O O . HOH Solvent 3 . ? 16 21 26 1 31 71 HOH A O 2 ? opaque8".split(),
]


def render(rows=ROWS):
    return (
        "data_original_chains\n# original synthetic input\n_entry.id original\nloop_\n"
        + "".join("_atom_site." + column + "\n" for column in COLUMNS)
        + "\n".join(" ".join(row) for row in rows)
        + "\n#\n"
    )


def test_model_and_label_chains_do_not_merge_author_aliases():
    text = render()
    result = C.read_mmcif_chains(text)
    assert result.source_text == text
    assert [(c.model_num, c.label_asym_id) for c in result.chains] == [
        (1, "Alpha"),
        (1, "Ligand"),
        (2, "Alpha"),
        (2, "Solvent"),
    ]
    assert [[a.source_index for a in c.atoms] for c in result.chains] == [
        [0, 1],
        [2, 3],
        [4, 5],
        [6, 7],
    ]
    assert all(a.residue.auth_asym_id == "A" for c in result.chains for a in c.atoms)
    assert "_atom_site.custom_payload" in result.unhandled_columns
    assert "opaque8" in result.source_text
    assert all(isinstance(c, C.MmcifChain) for c in result.chains)


def test_long_aliases_alternates_insertion_charge_and_numeric_fields():
    result = C.read_mmcif_chains(render())
    atoms = result.chains[1].atoms
    assert [a.alt_loc for a in atoms] == ["AA", "BB"]
    assert atoms[0].label_atom_id == "N1"
    assert atoms[0].auth_atom_id == "AN1"
    assert atoms[0].position_angstrom == (7.125, 0.0, 1.0)
    assert (atoms[0].occupancy, atoms[0].b_iso, atoms[0].formal_charge) == (
        0.4,
        22.5,
        1,
    )
    residue = atoms[0].residue
    assert residue is atoms[1].residue
    assert isinstance(residue, C.MmcifChainResidue)
    assert (residue.label_comp_id, residue.auth_comp_id) == ("A1AJ1", "AUTHOR_LONG")
    assert residue.label_seq_id is None
    assert (residue.auth_seq_id, residue.insertion_code, residue.entity_id) == (
        "50",
        "AB",
        "2",
    )
    waters = result.chains[3].atoms
    assert [a.residue.auth_seq_id for a in waters] == ["70", "71"]
    assert waters[0].residue != waters[1].residue


def test_all_720_interleavings_preserve_chain_and_source_correspondence():
    for order in itertools.permutations(range(6)):
        rows = [ROWS[i] for i in order]
        result = C.read_mmcif_chains(render(rows))
        expected_keys = list(dict.fromkeys((int(r[19]), r[6]) for r in rows))
        assert [(c.model_num, c.label_asym_id) for c in result.chains] == expected_keys
        seen = []
        for chain in result.chains:
            indices = [
                i
                for i, r in enumerate(rows)
                if (int(r[19]), r[6]) == (chain.model_num, chain.label_asym_id)
            ]
            assert [a.source_index for a in chain.atoms] == indices
            for atom in chain.atoms:
                row = rows[atom.source_index]
                assert atom.position_angstrom == tuple(float(x) for x in row[10:13])
                assert atom.serial == int(row[1])
                assert atom.residue.label_comp_id == row[5]
                seen.append(atom.source_index)
        assert sorted(seen) == list(range(6))


def test_outputs_are_deeply_immutable_and_repeat_reads_are_independent():
    text = render()
    first = C.read_mmcif_chains(text)
    second = C.read_mmcif_chains(text)
    assert first == second and first is not second
    for obj, name, replacement in [
        (first, "source_text", "broken"),
        (first.chains[0], "atoms", ()),
        (first.chains[0].atoms[0], "serial", 99),
        (first.chains[0].atoms[0].residue, "label_asym_id", "X"),
    ]:
        with pytest.raises(FrozenInstanceError):
            setattr(obj, name, replacement)
    assert isinstance(first.chains, tuple)
    assert isinstance(first.chains[0].atoms, tuple)
    assert isinstance(first.chains[0].atoms[0].position_angstrom, tuple)
    assert C.read_mmcif_chains(text) == second


@pytest.mark.parametrize("label", ["?", ".", "''", "' '"])
def test_missing_label_never_falls_back_to_author_alias(label):
    rows = [ROWS[0].copy()]
    rows[0][6] = label
    with pytest.raises(ValueError, match="requires label_asym_id"):
        C.read_mmcif_chains(render(rows))


@pytest.mark.parametrize(
    "column,value",
    [(19, "0"), (19, "4294967297"), (10, "NaN"), (10, "inf"), (2, "Xx"), (15, "12A")],
)
def test_native_invalid_fields_fail_explicitly(column, value):
    rows = [ROWS[0].copy()]
    rows[0][column] = value
    with pytest.raises(ValueError):
        C.read_mmcif_chains(render(rows))


def test_duplicate_atoms_and_incomplete_rows_are_refused():
    with pytest.raises(ValueError, match="repeats"):
        C.read_mmcif_chains(render([ROWS[0], ROWS[0]]))
    with pytest.raises(ValueError, match="incomplete"):
        C.read_mmcif_chains(render([ROWS[0][:-1]]))


@pytest.mark.parametrize("value", [None, b"data_x", 12, [], True])
def test_input_is_contents_str_only(value):
    with pytest.raises(TypeError):
        C.read_mmcif_chains(value)


@pytest.mark.parametrize(
    "name,maximum",
    [
        ("max_input_bytes", 8 * 1024 * 1024),
        ("max_atoms", 50000),
        ("max_chains", 4096),
        ("max_line_len", 8192),
    ],
)
@pytest.mark.parametrize(
    "value,exception",
    [
        (0, ValueError),
        (-1, ValueError),
        (True, TypeError),
        (1.0, TypeError),
        ("2", TypeError),
        (None, TypeError),
    ],
)
def test_invalid_limits(name, maximum, value, exception):
    with pytest.raises(exception):
        C.read_mmcif_chains(render(), **{name: value})
    with pytest.raises(ValueError):
        C.read_mmcif_chains(render(), **{name: maximum + 1})


@pytest.mark.parametrize(
    "kwargs",
    [
        {"max_input_bytes": 32},
        {"max_atoms": 7},
        {"max_chains": 3},
        {"max_line_len": 16},
    ],
)
def test_exceeded_limits_refuse_entire_read(kwargs):
    with pytest.raises(ValueError):
        C.read_mmcif_chains(render(), **kwargs)
    assert len(C.read_mmcif_chains(render()).chains) == 4


def test_utf8_byte_and_exact_count_boundaries_with_numpy_integers():
    text = render().replace("# original synthetic input", "# 元データ")
    result = C.read_mmcif_chains(
        text,
        max_input_bytes=np.int64(len(text.encode())),
        max_atoms=np.int64(8),
        max_chains=np.int64(4),
    )
    assert result.source_text == text
    with pytest.raises(ValueError, match="byte limit"):
        C.read_mmcif_chains(text, max_input_bytes=len(text.encode()) - 1)


def test_missing_optional_author_aliases_remain_none():
    row = ROWS[0].copy()
    for i in [15, 16, 17, 18]:
        row[i] = "?"
    result = C.read_mmcif_chains(render([row])).chains[0].atoms[0]
    assert result.auth_atom_id is None
    assert result.residue.auth_asym_id is None
    assert result.residue.auth_comp_id is None
    assert result.residue.auth_seq_id is None
    assert result.residue.label_seq_id == 1


def test_native_first_loop_scope_is_explicit_and_source_keeps_other_blocks():
    text = render(ROWS[:1]) + render(ROWS[2:3]).replace(
        "data_original_chains", "data_second"
    )
    result = C.read_mmcif_chains(text)
    assert len(result.chains) == 1
    assert len(result.chains[0].atoms) == 1
    assert result.source_text == text
