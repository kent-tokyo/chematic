"""Accounting checks for the independent RDKit C++ old/new packet."""

import importlib.util
import gzip
import json
from pathlib import Path

import pytest


SCRIPT = Path(__file__).resolve().parents[1] / "summarize_rdkit_native_rebaseline.py"
SPEC = importlib.util.spec_from_file_location("native_summary", SCRIPT)
assert SPEC is not None and SPEC.loader is not None
MODULE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(MODULE)


def row(index: int = 0) -> dict:
    return {
        "input_index": index,
        "smiles": "CC",
        "status": "ok",
        "canonical": "CC",
        "cip_atoms": {},
        "cip_bonds": {},
        "morgan_on_bits": [2, 5],
        "smarts": [[[0]], []],
    }


def test_counts_each_changed_cell_without_losing_unchanged_rows():
    old = [row(0), row(1)]
    new = [row(0), row(1)]
    new[1]["canonical"] = "C-C"
    new[1]["smarts"] = [[], []]
    counts, deltas = MODULE.summarize(old, new, ["C", "N"])
    assert counts == {
        "rows": 2,
        "smarts_cells_compared": 4,
        "changed_rows": 1,
        "canonical_changed_rows": 1,
        "smarts_changed_cells": 1,
        "noncomparable_rows": 0,
        "parse_status_changed_rows": 0,
        "cip_atoms_changed_rows": 0,
        "cip_bonds_changed_rows": 0,
        "morgan_on_bits_changed_rows": 0,
    }
    assert deltas[0]["input_index"] == 1
    assert [item["operation"] for item in deltas[0]["differences"]] == [
        "canonical",
        "smarts",
    ]


def test_incomplete_or_reordered_rows_fail(tmp_path: Path):
    path = tmp_path / "rows.jsonl"
    path.write_text(json.dumps(row(1)) + "\n", encoding="utf-8")
    with pytest.raises(ValueError, match="reordered"):
        MODULE.read_rows(path, ["CC"], 2)
    path.write_text(json.dumps(row(0)) + "\n", encoding="utf-8")
    with pytest.raises(ValueError, match="expected 2 rows"):
        MODULE.read_rows(path, ["CC", "CC"], 2)


def test_invalid_bit_set_fails(tmp_path: Path):
    path = tmp_path / "rows.jsonl"
    bad = row()
    bad["morgan_on_bits"] = [5, 2]
    path.write_text(json.dumps(bad) + "\n", encoding="utf-8")
    with pytest.raises(ValueError, match="invalid Morgan"):
        MODULE.read_rows(path, ["CC"], 2)


def test_new_parse_failure_is_a_reported_oracle_change(tmp_path: Path):
    path = tmp_path / "rows.jsonl"
    path.write_text(
        json.dumps({"input_index": 0, "smiles": "CC", "status": "parse_failure"})
        + "\n",
        encoding="utf-8",
    )
    new = MODULE.read_rows(path, ["CC"], 2)
    counts, deltas = MODULE.summarize([row()], new, ["C", "N"])
    assert counts["rows"] == 1
    assert counts["noncomparable_rows"] == 1
    assert counts["parse_status_changed_rows"] == 1
    assert counts["smarts_cells_compared"] == 0
    assert deltas[0]["differences"][0]["operation"] == "smiles_parse"


def test_old_native_oracle_must_match_python_wheel_oracle(tmp_path: Path):
    native = row()
    baseline = {
        "input_index": 0,
        "smiles": "CC",
        "status": "completed",
        "smiles_parse_write": {"rdkit_canonical": "CC"},
        "cip": {"rdkit_atoms": {}, "rdkit_bonds": {}},
        "morgan": {"rdkit_sha256": MODULE.morgan_binary_sha256([2, 5])},
    }
    path = tmp_path / "python.jsonl.gz"
    with gzip.open(path, "wt", encoding="utf-8") as handle:
        handle.write(json.dumps(baseline) + "\n")
    assert MODULE.verify_python_baseline(path, [native]) == 1
    native["canonical"] = "C-C"
    with pytest.raises(ValueError, match="canonical differs"):
        MODULE.verify_python_baseline(path, [native])
