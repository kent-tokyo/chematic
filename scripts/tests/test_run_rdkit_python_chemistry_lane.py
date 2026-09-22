import json
from pathlib import Path

import pytest

from scripts.run_rdkit_python_chemistry_lane import (
    difference,
    example_record,
    load_queries,
    load_smiles,
    normalized_match_sets,
)


def test_loaders_keep_exposed_rows_and_queries_deterministic(tmp_path: Path):
    corpus = tmp_path / "rows.smi"
    corpus.write_text("CCO first\n# comment\nCCN second\n", encoding="utf-8")
    assert load_smiles(corpus, None) == ["CCO", "CCN"]
    assert load_smiles(corpus, 1) == ["CCO"]

    queries = tmp_path / "queries.json"
    queries.write_text(json.dumps({"queries": ["[OH]", "[#7]"]}), encoding="utf-8")
    assert load_queries(queries) == ["[OH]", "[#7]"]


def test_query_loader_rejects_duplicates(tmp_path: Path):
    path = tmp_path / "queries.json"
    path.write_text(json.dumps({"queries": ["C", "C"]}), encoding="utf-8")
    with pytest.raises(ValueError, match="duplicate"):
        load_queries(path)


def test_match_sets_compare_target_atoms_not_embedding_order():
    assert normalized_match_sets([(2, 1), (1, 2), (4,)]) == [[1, 2], [4]]


def test_difference_defaults_to_unresolved_but_accepts_adjudication():
    assert difference("cip", "labels differ")["classification"] == "unresolved"
    assert (
        difference("smiles_parse_write", "spelling only", "contract_difference")[
            "classification"
        ]
        == "contract_difference"
    )


def test_summary_examples_do_not_duplicate_full_operation_payloads():
    row = {
        "input_index": 3,
        "smiles": "CCO",
        "status": "completed",
        "differences": [{"operation": "smarts"}],
        "smarts": {"large": [1, 2, 3]},
    }
    assert example_record(row) == {
        "input_index": 3,
        "smiles": "CCO",
        "status": "completed",
        "differences": [{"operation": "smarts"}],
    }
