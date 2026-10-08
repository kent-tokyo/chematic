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


CHEMATIC_SIDE = [
    "reaction_python_checked_candidates.py",
    "chematic_chemistry_dump.py",
    "check_python_smarts_parity_310k.py",
    "biotransformer_chematic_responses.py",
    "biotransformer_rules.py",
    "chematic_reaction_worker.py",
    "fetch_biotransformer_rules.py",
]


def test_chematic_side_scripts_import_no_rdkit_at_module_level():
    """The macOS and Windows jobs run these without RDKit installed."""
    import ast

    scripts = SCRIPT.parent
    for name in CHEMATIC_SIDE:
        tree = ast.parse((scripts / name).read_text(encoding="utf-8"))
        for node in tree.body:
            names = []
            if isinstance(node, ast.Import):
                names = [a.name for a in node.names]
            elif isinstance(node, ast.ImportFrom):
                names = [node.module or ""]
            assert not any(n.split(".")[0] == "rdkit" for n in names), name


def test_replay_answers_recorded_requests(tmp_path):
    import gzip
    import sys

    sys.path.insert(0, str(SCRIPT.parent))
    spec = importlib.util.spec_from_file_location("bt_corpus", SCRIPT.parent / "biotransformer_rule_corpus.py")
    corpus = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(corpus)

    rules = tmp_path / "rules.json"
    rules.write_text('{"reactions": {"one": {"btmrID": "R1", "smirks": "[C:1]>>[C:1]O"},'
                     ' "two": {"btmrID": "R2", "smirks": "[N:1]>>[N:1]C"}}}', encoding="utf-8")
    requests = tmp_path / "requests.tsv"
    requests.write_text("# header\nCC\t[H]C([H])([H])C([H])([H])[H]\n", encoding="utf-8")
    responses = tmp_path / "responses.jsonl.gz"
    with gzip.open(responses, "wt", encoding="utf-8") as f:
        f.write(json.dumps({"version": "0.0", "file": "x"}) + "\n")
        f.write(json.dumps([0, 0, "ok", None, [[{"smiles": "CCO"}]]]) + "\n")
        f.write(json.dumps([1, 1, "no_match", None, []]) + "\n")
    replay = corpus.Replay([responses], [rules], requests)
    assert replay.version == "0.0"
    got = replay.run({"smirks": "[C:1]>>[C:1]O", "smiles": "CC"})
    assert got == {"status": "ok", "detail": None, "products": [[{"smiles": "CCO"}]]}
    got = replay.run({"smirks": "[N:1]>>[N:1]C", "smiles": "[H]C([H])([H])C([H])([H])[H]"})
    assert got["status"] == "no_match"
    try:
        replay.run({"smirks": "[C:1]>>[C:1]O", "smiles": "CCC"})
    except SystemExit:
        pass
    else:
        raise AssertionError("an unrecorded request must stop the run")
