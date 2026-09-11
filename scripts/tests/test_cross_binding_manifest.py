"""Regression tests for the shared cross-binding manifest loader."""

import json
from pathlib import Path

import pytest

from scripts.check_cross_binding_manifest import (
    DuplicateJsonKeyError,
    _reject_duplicate_keys,
)


def test_current_manifest_has_no_duplicate_keys():
    path = Path(__file__).parents[2] / "validation" / "cross_binding_contract.json"
    document = json.loads(
        path.read_text(encoding="utf-8"), object_pairs_hook=_reject_duplicate_keys
    )
    assert document["operation_manifest"]["schema_version"] == 1


def test_manifest_loader_rejects_duplicate_keys():
    with pytest.raises(DuplicateJsonKeyError, match="duplicate JSON key: contract"):
        json.loads(
            '{"contract": "first", "contract": "shadowed"}',
            object_pairs_hook=_reject_duplicate_keys,
        )
