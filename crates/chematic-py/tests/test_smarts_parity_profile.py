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


def test_order_dependent_ring_count_follows_rdkit_ring_order():
    # RDKit 2026.03's symmetrized SSSR depends on atom order for this
    # macrocycle; the default profile counts over the exact port of RDKit's
    # ring order (RDKit 2026.03.1 and 2026.03.6: atoms 3 and 4 are in three
    # rings) instead of refusing.
    corpus = Path(__file__).resolve().parents[3] / "validation/benchmark_corpora/rdkit-js-browser-10k-v1.smi"
    smiles = corpus.read_text(encoding="utf-8").splitlines()[23]
    mol = chematic.from_smiles(smiles)
    result = mol.find_matches_rdkit_parity("[R3]")
    assert result == {"status": "ok", "reason": None, "matches": [[3], [4]]}


def test_invalid_smarts_still_raises():
    with pytest.raises(ValueError, match="invalid SMARTS"):
        chematic.from_smiles("CCO").find_matches_rdkit_parity("[")


def test_rdkit_2026_09_1_profile_counts_relevant_cycles():
    # RDKit 2026.09.1 (native C++ build): no atom of this macrocycle is in
    # three rings, eight are in exactly one.
    corpus = Path(__file__).resolve().parents[3] / "validation/benchmark_corpora/rdkit-js-browser-10k-v1.smi"
    mol = chematic.from_smiles(corpus.read_text(encoding="utf-8").splitlines()[23])
    assert mol.find_matches_rdkit_parity("[R3]", profile="2026.09.1") == {
        "status": "ok", "reason": None, "matches": [],
    }
    r1 = mol.find_matches_rdkit_parity("[R1]", profile="2026.09.1")
    assert r1["matches"] == [[0], [1], [2], [5], [29], [30], [31], [32]]


def test_unknown_profile_raises():
    with pytest.raises(ValueError, match="unknown profile"):
        chematic.from_smiles("CCO").find_matches_rdkit_parity("[R1]", profile="2027.03.1")
