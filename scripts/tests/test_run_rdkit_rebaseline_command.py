import pytest

from scripts.run_rdkit_rebaseline_command import expand, find_command, parse_values


def test_parse_values_and_expand_keep_each_argv_token_separate():
    values = parse_values(["PYTHON=/tmp/venv/bin/python", "OUTPUT_DIR=/tmp/out dir"])
    assert expand(["${PYTHON}", "--output", "${OUTPUT_DIR}/result.json"], values) == [
        "/tmp/venv/bin/python",
        "--output",
        "/tmp/out dir/result.json",
    ]


def test_expand_rejects_missing_placeholders():
    with pytest.raises(ValueError, match="OUTPUT_DIR"):
        expand(["${PYTHON}", "${OUTPUT_DIR}/result.json"], {"PYTHON": "python3"})


def test_command_catalog_exposes_python_chemistry_lane():
    lane, command = find_command(
        "python-2026.03.6-distributed", "python-chemistry-all"
    )
    assert lane["rdkit_version"] == "2026.03.6"
    assert set(command["covers"]) == {
        "smiles_parse_write",
        "cip",
        "smarts",
        "morgan",
    }
