"""Independent synthetic PDB model splitting contracts; no external fixtures."""

from dataclasses import FrozenInstanceError
from pathlib import Path

import numpy as np
import pytest

import chematic

RICH = (
    Path(__file__).resolve().parents[3]
    / "validation/openbabel_issues/pdb-selection/rich.pdb"
).read_text(encoding="utf-8")
ROWS = (
    "\n".join(
        line
        for line in RICH.splitlines()
        if not line.startswith("CONECT") and line != "END"
    )
    + "\n"
)
CONECT = "CONECT    1    2    5    5\n"


def models(count, start=1, step=1, rows=ROWS):
    return (
        "".join(f"MODEL     {start+i*step:4d}\n{rows}ENDMDL\n" for i in range(count))
        + "END\n"
    )


@pytest.mark.parametrize(
    "count,width", [(1, 1), (9, 1), (10, 2), (99, 2), (100, 3), (1000, 4), (9999, 4)]
)
def test_padding_and_lexical_sort(count, width):
    r = chematic.split_pdb_models(models(count), prefix="nmr")
    assert len(r.files) == count
    assert [f.filename for f in r.files] == sorted(f.filename for f in r.files)
    assert r.files[0].filename == f"nmr{1:0{width}d}.pdb"
    assert r.files[-1].filename == f"nmr{count:0{width}d}.pdb"
    assert r.files[-1].model_num == count
    assert all(f.atom_count == 6 for f in r.files)
    assert all(f.pdb.splitlines()[1:-2] == ROWS.splitlines() for f in r.files)
    assert not r.losses


@pytest.mark.parametrize("crlf", [False, True])
@pytest.mark.parametrize("final_end", [False, True])
def test_original_identity_envelopes_shared_conect_and_source(crlf, final_end):
    text = models(3, start=7, step=12).removesuffix("END\n") + CONECT
    if final_end:
        text += "END\n"
    if crlf:
        text = text.replace("\n", "\r\n")
    result = chematic.split_pdb_models(text)
    assert result.source_text == text
    assert [f.model_num for f in result.files] == [7, 19, 31]
    for f in result.files:
        assert "\r" not in f.pdb
        assert f.pdb.endswith("ENDMDL\n" + CONECT + "END\n")
        checked = chematic.select_pdb_atoms(f.pdb)
        assert checked["pdb"] == f.pdb
        assert not checked["retargeted_ter_lines"]
        assert not checked["filtered_conect_lines"]
        assert [a["serial"] for a in checked["atoms"]] == [1, 2, 5, 6, 9, 10]
        assert checked["atoms"][0]["alt_loc"] == "A"
        assert checked["atoms"][2]["residue_number"] == -2
        assert checked["atoms"][2]["insertion_code"] == "I"


def test_implicit_model_and_no_input_newline():
    r = chematic.split_pdb_models(ROWS.rstrip("\n"), prefix="run-2026")
    assert r.files == (chematic.PdbModelFile("run-20261.pdb", None, ROWS + "END\n", 6),)
    with pytest.raises(FrozenInstanceError):
        r.files[0].filename = "oops"
    with pytest.raises(FrozenInstanceError):
        r.source_text = "oops"


def test_metadata_requires_consent_and_loss_lines_reference_original():
    text = "HEADER    original\nNUMMDL       2\n" + models(2).replace(
        "ENDMDL\n", "REMARK    original inside model\nENDMDL\n"
    )
    text = text.replace("END\n", "MASTER    original count\nEND\n")
    with pytest.raises(ValueError, match="consent"):
        chematic.split_pdb_models(text)
    result = chematic.split_pdb_models(text, allow_metadata_loss=True)
    assert result.source_text == text
    assert [loss.record_name for loss in result.losses] == [
        "HEADER",
        "NUMMDL",
        "REMARK",
        "REMARK",
        "MASTER",
    ]
    for loss in result.losses:
        assert text.splitlines()[loss.line - 1].startswith(loss.record_name)
    assert all("REMARK" not in f.pdb and "MASTER" not in f.pdb for f in result.files)


@pytest.mark.parametrize(
    "text,match",
    [
        ("", "at least one"),
        ("\n\n", "at least one"),
        ("ENDMDL\n", "unmatched"),
        ("MODEL        1\nENDMDL\n", "no ATOM"),
        ("MODEL        1\n" + ROWS, "lacks ENDMDL"),
        ("MODEL        1\nMODEL        2\n" + ROWS + "ENDMDL\n", "nested"),
        (models(2, start=1, step=0), "duplicate MODEL"),
        (models(1) + ROWS, "after END"),
        (models(1).replace("END\n", ROWS + "END\n"), "outside model"),
        (ROWS + models(1), "mixed implicit"),
        (models(1).replace("ENDMDL", "CONECT    1    2\nENDMDL"), "CONECT must follow"),
        (models(1).replace("END\n", CONECT + models(1, start=2)), "nested MODEL"),
        (models(1).replace("ENDMDL", "ENDMDL x"), "malformed ENDMDL"),
        (models(1).replace("END\n", "END   x\n"), "malformed END"),
        (models(1).replace("MODEL        1", "MODEL        0"), "MODEL number"),
        (models(1).replace("MODEL        1", "MODEL       -1"), "MODEL number"),
        (models(1).replace("MODEL        1", "MODEL        1x"), "MODEL number"),
        (models(1).replace("MODEL        1", "MODEL x      1"), "MODEL number"),
        (models(1).replace("MODEL        1", "MODEL     10000"), "MODEL number"),
        (models(1).replace("MODEL        1", "MODEL"), "MODEL number"),
        (models(1).replace("END\n", "CONECT    1  999\nEND\n"), "CONECT"),
        (models(1).replace("   0.000", "     nan", 1), "invalid PDB x"),
        (models(1).replace("ATOM      2", "ATOM      1"), "duplicate"),
        (models(1).replace("ANISOU    2", "ANISOU   99"), "ANISOU"),
        (models(1).replace("TER       3", "TER      99"), "TER"),
        (models(1) + "\t", "printable ASCII"),
        (models(1) + "REMARK café\n", "ASCII"),
        (models(1) + "X" * 1025, "1024"),
    ],
)
def test_bad_input_fails_even_with_metadata_loss_consent(text, match):
    with pytest.raises(ValueError, match=match):
        chematic.split_pdb_models(text, allow_metadata_loss=True)


@pytest.mark.parametrize(
    "prefix",
    ["", "../split", "/tmp/file", "a/b", "a\\b", ".hidden", "café", "x\x00", "x" * 65],
)
def test_filename_stems(prefix):
    with pytest.raises(ValueError, match="prefix"):
        chematic.split_pdb_models(models(1), prefix=prefix)


@pytest.mark.parametrize(
    "name,number",
    [
        ("max_models", 0),
        ("max_models", 10000),
        ("max_output_bytes", 0),
        ("max_output_bytes", 64 * 1024 * 1024 + 1),
    ],
)
def test_limit_ranges(name, number):
    with pytest.raises(ValueError, match=name):
        chematic.split_pdb_models(models(1), **{name: number})


@pytest.mark.parametrize("name", ["max_models", "max_output_bytes"])
@pytest.mark.parametrize("number", [True, 1.5, "2"])
def test_limit_types(name, number):
    with pytest.raises(TypeError):
        chematic.split_pdb_models(models(1), **{name: number})


def test_output_expansion_budget_and_model_cap():
    text = models(2).replace("END\n", CONECT + "END\n")
    expected = sum(len(f.pdb) for f in chematic.split_pdb_models(text).files)
    assert (
        sum(
            len(f.pdb)
            for f in chematic.split_pdb_models(
                text, max_output_bytes=np.int64(expected)
            ).files
        )
        == expected
    )
    with pytest.raises(ValueError, match="aggregate output"):
        chematic.split_pdb_models(text, max_output_bytes=expected - 1)
    with pytest.raises(ValueError, match="model count"):
        chematic.split_pdb_models(text, max_models=1)
    assert len(chematic.split_pdb_models(text, max_models=np.int32(2)).files) == 2


def test_all_models_validate_before_result_or_io(tmp_path):
    good = models(1).removesuffix("END\n")
    bad = models(1, start=2).replace("   0.000", "     nan", 1)
    with pytest.raises(ValueError, match="output 2"):
        chematic.split_pdb_models(good + bad)
    assert not list(tmp_path.iterdir())


def test_shared_conect_must_resolve_in_every_model():
    short = (
        "\n".join(
            line
            for line in ROWS.splitlines()
            if line.startswith(("ATOM      1", "ANISOU    1"))
        )
        + "\n"
    )
    text = models(1).removesuffix("END\n") + models(1, start=2, rows=short).replace(
        "END\n", CONECT + "END\n"
    )
    with pytest.raises(ValueError, match="output 2"):
        chematic.split_pdb_models(text)


@pytest.mark.parametrize("text", [b"END\n", None, 1])
def test_text_type(text):
    with pytest.raises(TypeError):
        chematic.split_pdb_models(text)


def test_loss_flag_type():
    with pytest.raises(TypeError, match="bool"):
        chematic.split_pdb_models(models(1), allow_metadata_loss=1)


def test_input_record_bound_before_blank_line_materialization():
    with pytest.raises(ValueError, match="record limit"):
        chematic.split_pdb_models("\n" * 1_000_001)


def test_input_size_bound():
    with pytest.raises(ValueError, match="input byte limit"):
        chematic.split_pdb_models(" " * (64 * 1024 * 1024 + 1))
