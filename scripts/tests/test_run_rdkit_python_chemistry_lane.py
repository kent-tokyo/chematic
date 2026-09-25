import json
from pathlib import Path

import pytest

from scripts.run_rdkit_python_chemistry_lane import (
    bond_endpoint_key,
    chematic_cip,
    count_difference_classes,
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


def test_difference_class_counts_count_each_operation_result():
    from collections import Counter

    counts: Counter[str] = Counter()
    count_difference_classes(
        [
            difference("smiles_parse_write", "spelling", "contract_difference"),
            difference("cip", "labels"),
            difference("smarts", "query cells"),
        ],
        counts,
    )
    assert counts == {
        "difference_class_contract_difference": 1,
        "difference_class_unresolved": 2,
    }


def test_bond_endpoint_key_is_independent_of_bond_index_and_direction():
    assert bond_endpoint_key(7, 2) == "2-7"
    assert bond_endpoint_key(2, 7) == "2-7"


def test_chematic_cip_reads_ez_bond_endpoints():
    # E/Z entries are keyed by the double bond's first atom (#634); the
    # explicit ``bond_atoms`` field wins when present.
    class FakeMolecule:
        bond_table = [
            (0, 1, "SINGLE", False),
            (4, 2, "DOUBLE", False),
            (5, 6, "DOUBLE", False),
        ]

        def cip_stereo(self, mode):
            assert mode == "accurate"
            return [
                {"atom_idx": 3, "descriptor": "R"},
                {"atom_idx": 4, "descriptor": "E"},
                {"atom_idx": 5, "descriptor": "Z", "bond_idx": 2, "bond_atoms": (5, 6)},
            ]

        def cip_stereo_unresolved(self):
            return []

    atoms, bonds, unresolved = chematic_cip(FakeMolecule())
    assert atoms == {3: "R"}
    assert bonds == {"2-4": "E", "5-6": "Z"}
    assert unresolved == {}


def test_chematic_cip_does_not_read_ez_atom_as_bond_index():
    # Historical comparator bug: atom 1 used to be read as bond 1 = (4, 2).
    class FakeMolecule:
        bond_table = [
            (0, 1, "SINGLE", False),
            (4, 2, "DOUBLE", False),
            (1, 7, "DOUBLE", False),
        ]

        def cip_stereo(self, mode):
            return [{"atom_idx": 1, "descriptor": "E"}]

        def cip_stereo_unresolved(self):
            return []

    _atoms, bonds, _unresolved = chematic_cip(FakeMolecule())
    assert bonds == {"1-7": "E"}


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
