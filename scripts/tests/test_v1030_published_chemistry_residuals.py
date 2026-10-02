from pathlib import Path

import pytest

from scripts import check_v1030_published_chemistry_residuals as bridge


def test_published_residual_classification_covers_every_unresolved_row():
    result = bridge.check()
    assert result["published_rows"] == 10_000
    assert result["unresolved_reclassified"] == 69
    assert result["smarts_cells"] == 200
    assert result["cip_typed_abstentions"] == {
        "lone_pair_center": 1,
        "phosphorus_oracle_unstable": 4,
    }


def test_published_archive_mutation_fails_closed(tmp_path: Path, monkeypatch):
    changed = tmp_path / "changed.jsonl.gz"
    changed.write_bytes(bridge.PUBLISHED.read_bytes() + b"unexpected")
    monkeypatch.setattr(bridge, "PUBLISHED", changed)
    with pytest.raises(ValueError, match="archive changed"):
        bridge.check()


def test_adjudication_mutation_fails_closed(tmp_path: Path, monkeypatch):
    changed = tmp_path / "classification.json"
    changed.write_bytes(bridge.CLASSIFICATION.read_bytes() + b" ")
    monkeypatch.setattr(bridge, "CLASSIFICATION", changed)
    with pytest.raises(ValueError, match="classification changed"):
        bridge.check()
