from __future__ import annotations

import importlib.util
from pathlib import Path

import pytest


elixcee = pytest.importorskip("elixcee")

ROOT = Path(__file__).resolve().parents[3]
EXAMPLE_PATH = ROOT / "examples" / "xlsx_batch_analysis.py"
SPEC = importlib.util.spec_from_file_location("xlsx_batch_analysis", EXAMPLE_PATH)
assert SPEC is not None and SPEC.loader is not None
xlsx_batch_analysis = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(xlsx_batch_analysis)


def make_input(path: Path) -> None:
    with elixcee.create_stream_bounded(str(path), max_rows=10, max_columns=3) as writer:
        writer.append(["ID", "SMILES", "note"])
        writer.append(["aspirin", "CC(=O)Oc1ccccc1C(=O)O", "ok"])
        writer.append(["bad", "C1CC", "invalid"])
        writer.append(["blank", "", "missing"])


def test_xlsx_batch_analysis_round_trip(tmp_path: Path) -> None:
    input_path = tmp_path / "input.xlsx"
    output_path = tmp_path / "output.xlsx"
    make_input(input_path)

    assert xlsx_batch_analysis.main([str(input_path), str(output_path)]) == 0

    with elixcee.open_stream(str(output_path), include_row_numbers=True) as reader:
        rows = list(reader)
    assert rows[0][1] == list(xlsx_batch_analysis.OUTPUT_COLUMNS)
    assert rows[1][1][0] == 2
    assert rows[1][1][13] == "ok"
    assert rows[2][1][13] == "error"
    assert rows[3][1][13] == "error"


def test_xlsx_batch_analysis_does_not_publish_truncated_output(tmp_path: Path) -> None:
    input_path = tmp_path / "input.xlsx"
    output_path = tmp_path / "output.xlsx"
    make_input(input_path)

    assert (
        xlsx_batch_analysis.main([str(input_path), str(output_path), "--max-rows", "2"])
        == 1
    )
    assert not output_path.exists()
