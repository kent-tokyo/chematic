"""Resume and corpus-shard helpers of the published-wheel gate driver."""

import importlib.util
import json
from pathlib import Path

SCRIPT = Path(__file__).resolve().parents[1] / "published_wheel_chemistry_gates.py"
SPEC = importlib.util.spec_from_file_location("published_wheel_chemistry_gates", SCRIPT)
assert SPEC and SPEC.loader
MODULE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(MODULE)


def shard(tmp_path: Path, name: str, reactants: int, counts: dict, per_rule: dict, rows: str):
    summary = tmp_path / f"{name}.json"
    summary.write_text(json.dumps({
        "schema": "biotransformer-rule-corpus/v1",
        "rule_count": 3,
        "reactants": {"file": "r.smi", "sha256": "x", "every": 1, "limit": None, "count": reactants},
        "elapsed_seconds": 1.5,
        "counts": counts,
        "per_rule": per_rule,
    }))
    row_file = tmp_path / f"{name}-rows.jsonl"
    row_file.write_text(rows)
    return summary, row_file


def test_reactant_shards_add_up(tmp_path):
    a = shard(tmp_path, "a", 2,
              {"rules_run": 2, "rdkit_parse_error_rules": 1, "implicit:exact": 3, "explicit_h:differ": 1},
              {"R1": {"implicit:exact": 2, "explicit_h:differ": 1}, "R2": {"implicit:exact": 1},
               "R3": {"rdkit_parse_error": True}},
              '{"rule": "R1"}\n')
    b = shard(tmp_path, "b", 3,
              {"rules_run": 2, "rdkit_parse_error_rules": 1, "implicit:exact": 4},
              {"R1": {"implicit:exact": 3}, "R2": {"implicit:exact": 1},
               "R3": {"rdkit_parse_error": True}},
              "")
    out, rows = tmp_path / "merged.json", tmp_path / "merged-rows.jsonl"
    MODULE.merge_corpus_shards([a, b], out, rows)
    merged = json.loads(out.read_text())
    assert merged["counts"] == {"explicit_h:differ": 1, "implicit:exact": 7,
                                "rdkit_parse_error_rules": 1, "rules_run": 2}
    assert merged["per_rule"]["R1"] == {"implicit:exact": 5, "explicit_h:differ": 1}
    assert merged["per_rule"]["R3"] == {"rdkit_parse_error": True}
    assert merged["reactants"]["count"] == 5
    assert merged["shards"] == 2
    assert rows.read_text() == '{"rule": "R1"}\n'


def test_resume_only_skips_complete_json(tmp_path):
    path = tmp_path / "step.json"
    assert not MODULE.finished(path)
    path.write_text('{"counts": {"exact": 1')  # cut off mid-write
    assert not MODULE.finished(path)
    path.write_text('{"counts": {"exact": 1}}')
    assert MODULE.finished(path)
