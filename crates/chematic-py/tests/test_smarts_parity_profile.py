"""The opt-in SMARTS profile never turns a refusal into an empty match set."""

from pathlib import Path

import pytest

import chematic


def test_completed_match_and_no_match_are_distinct_from_refusal():
    mol = chematic.from_smiles("CCO")
    assert mol.find_matches_rdkit_parity("[#8]") == {
        "status": "ok", "reason": None, "matches": [[2]],
    }
    assert mol.find_matches_rdkit_parity("[#9]") == {
        "status": "ok", "reason": None, "matches": [],
    }


def test_ambiguous_ring_count_is_typed_unsupported():
    corpus = Path(__file__).resolve().parents[3] / "validation/benchmark_corpora/rdkit-js-browser-10k-v1.smi"
    smiles = corpus.read_text(encoding="utf-8").splitlines()[23]
    mol = chematic.from_smiles(smiles)
    result = mol.find_matches_rdkit_parity("[R3]")
    assert result == {
        "status": "typed_unsupported",
        "reason": "ring_model_ambiguous",
        "matches": None,
    }


def test_invalid_smarts_still_raises():
    with pytest.raises(ValueError, match="invalid SMARTS"):
        chematic.from_smiles("CCO").find_matches_rdkit_parity("[")
