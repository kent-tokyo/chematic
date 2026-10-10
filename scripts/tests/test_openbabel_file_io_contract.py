"""Regression tests for the bounded Open Babel file-I/O contract."""

from __future__ import annotations

import copy
import json
from pathlib import Path

from scripts.check_openbabel_file_io_contract import validate_contract


ROOT = Path(__file__).resolve().parents[2]
CONTRACT = json.loads(
    (ROOT / "validation" / "openbabel_file_io_contract_v1.json").read_text(encoding="utf-8")
)


def test_checked_in_contract_is_valid():
    assert validate_contract(CONTRACT) == []


def test_contract_rejects_a_missing_tier_a_format():
    candidate = copy.deepcopy(CONTRACT)
    candidate["formats"] = candidate["formats"][:-1]
    assert any("format ids must be exactly" in error for error in validate_contract(candidate))


def test_contract_rejects_a_non_equivalent_speed_claim():
    candidate = copy.deepcopy(CONTRACT)
    candidate["win_criteria"]["performance"]["require_equivalent_output"] = False
    assert "performance wins must require equivalent output" in validate_contract(candidate)
